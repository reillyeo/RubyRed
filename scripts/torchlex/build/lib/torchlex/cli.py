import argparse
import glob
import os
import shutil
import sys
import time
from multiprocessing import Pool

from torchlex.config import load_config
from torchlex.demultiplexer import Demultiplexer
from torchlex.fastq_reader import FastqFile


def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="torchlex",
        description=(
            "Demultiplex Oxford Nanopore FASTQ files by barcode.\n\n"
            "Reads are assigned to one of up to 96 samples using Levenshtein-distance\n"
            "matching against spacer/linker anchor sequences flanking each barcode."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", nargs="+", required=True, metavar="PATH",
                   help="Input FASTQ/FASTQ.GZ files or a directory containing them")
    p.add_argument("-o", "--output", required=True,
                   help="Output directory")
    p.add_argument("-t", "--threads", type=int, default=1,
                   help="Number of parallel worker processes (default: 1)")
    p.add_argument("-c", "--config",
                   help="Path to a JSON config file (overrides built-in defaults)")
    p.add_argument("--barcodes",
                   help="Path to a custom barcodes FASTA (overrides bundled barcodes_custom2.fasta)")
    p.add_argument("--strict-filenames", action="store_true",
                   help="Enforce ONT filename format {flow_cell_id}_pass_{runid}_{idx}.fastq")
    p.add_argument("--min-qual", type=float, metavar="Q",
                   help="Minimum average Phred quality score (default: 8)")
    p.add_argument("--min-len", type=int, metavar="N",
                   help="Minimum read length after cropping (default: 900)")
    p.add_argument("--max-len", type=int, metavar="N",
                   help="Maximum read length (default: 1500)")
    p.add_argument("--max-dist", type=int, metavar="D",
                   help="Max Levenshtein distance for barcode assignment (default: 5, exclusive)")
    p.add_argument("--min-reads", type=int, metavar="N", default=1000,
                   help="Barcodes with fewer than N reads are moved to failed_barcodes/ (default: 1000)")
    return p.parse_args(argv)


def _resolve_inputs(paths):
    resolved = []
    for p in paths:
        if os.path.isdir(p):
            resolved.extend(sorted(glob.glob(os.path.join(p, "*.fastq"))))
            resolved.extend(sorted(glob.glob(os.path.join(p, "*.fastq.gz"))))
        elif os.path.isfile(p):
            resolved.append(p)
        else:
            matched = sorted(glob.glob(p))
            if not matched:
                print(f"warning: no files match '{p}'", file=sys.stderr)
            resolved.extend(matched)
    return resolved


def _chunk(lst, n):
    """Split lst into n roughly equal chunks."""
    k, m = divmod(len(lst), n)
    chunks, i = [], 0
    for j in range(n):
        size = k + (1 if j < m else 0)
        chunks.append(lst[i:i + size])
        i += size
    return [c for c in chunks if c]


# Top-level function required for multiprocessing pickling
def _worker(args):
    file_paths, config, strict_filenames, tmp_dir = args
    demux = Demultiplexer(config)
    n_reads = n_classified = n_bases = 0
    for path in file_paths:
        try:
            fq = FastqFile(path, strict_filenames=strict_filenames)
        except Exception as exc:
            print(f"[worker] skip {os.path.basename(path)}: {exc}", file=sys.stderr)
            continue
        classified = demux.demultiplex_file(fq, tmp_dir)
        n_reads += len(fq)
        n_classified += len(classified)
        n_bases += sum(len(r.seq) for r in fq.rows)
    return n_reads, n_classified, n_bases


def _merge(worker_dirs, output_dir):
    """Concatenate per-barcode files from each worker temp dir into output_dir."""
    barcode_files = set()
    for d in worker_dirs:
        for f in glob.glob(os.path.join(d, "BRK*.fastq")):
            barcode_files.add(os.path.basename(f))

    os.makedirs(output_dir, exist_ok=True)
    for name in sorted(barcode_files):
        with open(os.path.join(output_dir, name), "wb") as out:
            for d in worker_dirs:
                src = os.path.join(d, name)
                if os.path.exists(src):
                    with open(src, "rb") as src_f:
                        shutil.copyfileobj(src_f, out)

    for d in worker_dirs:
        shutil.rmtree(d, ignore_errors=True)


def _triage(output_dir, min_reads):
    """Move barcode FASTQ files below min_reads into failed_barcodes/."""
    failed_dir = os.path.join(output_dir, "failed_barcodes")
    passed, failed = [], []
    for f in sorted(glob.glob(os.path.join(output_dir, "BRK*.fastq"))):
        n = sum(1 for _ in open(f)) // 4
        if n < min_reads:
            os.makedirs(failed_dir, exist_ok=True)
            shutil.move(f, os.path.join(failed_dir, os.path.basename(f)))
            failed.append((os.path.basename(f), n))
        else:
            passed.append((os.path.basename(f), n))
    return passed, failed


def main(argv=None):
    args = _parse_args(argv)
    config = load_config(args.config)

    if args.barcodes:
        config["amplicones"]["barcodes_path"] = args.barcodes
    if args.min_qual is not None:
        config["sequence"]["minimum_average_phred_quality"] = args.min_qual
    if args.min_len is not None:
        config["sequence"]["min_sequence_length"] = args.min_len
    if args.max_len is not None:
        config["sequence"]["max_sequence_length"] = args.max_len
    if args.max_dist is not None:
        config["max_barcode_distance"] = args.max_dist

    fastq_paths = _resolve_inputs(args.input)
    if not fastq_paths:
        print("error: no input FASTQ files found", file=sys.stderr)
        sys.exit(1)

    n_workers = min(args.threads, len(fastq_paths))
    chunks = _chunk(fastq_paths, n_workers)
    tmp_dirs = [
        os.path.join(args.output, f".tmp_worker_{i}")
        for i in range(len(chunks))
    ]
    worker_args = [
        (chunk, config, args.strict_filenames, tmp_dir)
        for chunk, tmp_dir in zip(chunks, tmp_dirs)
    ]

    print(
        f"[torchlex] {len(fastq_paths)} files · {n_workers} worker(s) · {args.output}",
        file=sys.stderr,
    )

    t_start = time.perf_counter()

    # --- demultiplex ---
    if n_workers == 1:
        results = [_worker(worker_args[0])]
    else:
        with Pool(processes=n_workers) as pool:
            results = pool.map(_worker, worker_args)

    total_reads = total_classified = total_bases = 0
    for n_reads, n_classified, n_bases in results:
        total_reads += n_reads
        total_classified += n_classified
        total_bases += n_bases

    pct = 100.0 * total_classified / total_reads if total_reads else 0.0
    print(
        f"[torchlex] {total_classified:,}/{total_reads:,} reads classified ({pct:.1f}%)",
        file=sys.stderr,
    )

    t_demux = time.perf_counter()

    # --- merge worker outputs ---
    print("[torchlex] merging worker outputs ...", file=sys.stderr)
    _merge(tmp_dirs, args.output)

    t_merge = time.perf_counter()

    # --- triage by read count ---
    print(f"[torchlex] triaging barcodes (min reads: {args.min_reads}) ...", file=sys.stderr)
    passed, failed = _triage(args.output, args.min_reads)
    print(f"  passed : {len(passed)} barcodes → {args.output}/", file=sys.stderr)
    if failed:
        print(
            f"  failed : {len(failed)} barcodes → {args.output}/failed_barcodes/",
            file=sys.stderr,
        )
        for name, n in failed:
            print(f"           {name}  ({n} reads)", file=sys.stderr)

    # --- timing report ---
    t_end = time.perf_counter()
    elapsed = t_end - t_start
    elapsed_min = elapsed / 60.0
    gbp = total_bases / 1e9
    rate = elapsed_min / gbp if gbp > 0 else float("inf")

    print(f"\n[torchlex] timing", file=sys.stderr)
    print(f"  demux   : {(t_demux - t_start):.1f}s", file=sys.stderr)
    print(f"  merge   : {(t_merge - t_demux):.1f}s", file=sys.stderr)
    print(f"  triage  : {(t_end - t_merge):.1f}s", file=sys.stderr)
    print(f"  total   : {elapsed_min:.2f} min ({elapsed:.1f}s)", file=sys.stderr)
    print(f"\n[torchlex] throughput", file=sys.stderr)
    print(f"  input   : {gbp:.3f} Gbp across {len(fastq_paths)} files", file=sys.stderr)
    print(f"  speed   : {rate:.2f} min/Gbp  ({n_workers} thread(s))", file=sys.stderr)


if __name__ == "__main__":
    main()
