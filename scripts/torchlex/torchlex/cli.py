import argparse
import glob
import os
import shutil
import sys
import time
from multiprocessing import Pool

from tqdm import tqdm

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


# -----------------------------------------------------------------
# Per-process worker state (initialised once per worker process)
# -----------------------------------------------------------------
_state = {}


def _init_worker(config, output_dir, strict_filenames):
    _state["demux"] = Demultiplexer(config)
    _state["tmp_dir"] = os.path.join(output_dir, f".tmp_{os.getpid()}")
    _state["strict"] = strict_filenames


def _worker_file(path):
    """Process a single FASTQ file; called once per file by the pool."""
    try:
        fq = FastqFile(path, strict_filenames=_state["strict"])
    except Exception as exc:
        tqdm.write(f"[torchlex] skip {os.path.basename(path)}: {exc}")
        return 0, 0, 0
    classified = _state["demux"].demultiplex_file(fq, _state["tmp_dir"])
    n_bases = sum(len(r.seq) for r in fq.rows)
    return len(fq), len(classified), n_bases


# -----------------------------------------------------------------
# Merge & triage
# -----------------------------------------------------------------

def _merge(output_dir):
    """Concatenate per-worker-PID barcode files into the final output dir."""
    worker_dirs = glob.glob(os.path.join(output_dir, ".tmp_*"))
    barcode_files = set()
    for d in worker_dirs:
        for f in glob.glob(os.path.join(d, "BRK*.fastq")):
            barcode_files.add(os.path.basename(f))

    for name in tqdm(sorted(barcode_files), desc="merge   ", unit="barcode"):
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
    files = sorted(glob.glob(os.path.join(output_dir, "BRK*.fastq")))
    for f in tqdm(files, desc="triage  ", unit="barcode"):
        n = sum(1 for _ in open(f)) // 4
        if n < min_reads:
            os.makedirs(failed_dir, exist_ok=True)
            shutil.move(f, os.path.join(failed_dir, os.path.basename(f)))
            failed.append((os.path.basename(f), n))
        else:
            passed.append((os.path.basename(f), n))
    return passed, failed


# -----------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------

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
    os.makedirs(args.output, exist_ok=True)

    print(
        f"[torchlex] {len(fastq_paths)} files · {n_workers} worker(s) · {args.output}",
    )

    t_start = time.perf_counter()

    # --- demultiplex ---
    total_reads = total_classified = total_bases = 0

    init_args = (config, args.output, args.strict_filenames)

    with tqdm(total=len(fastq_paths), desc="demux   ", unit="file") as pbar:
        if n_workers == 1:
            _init_worker(*init_args)
            for path in fastq_paths:
                n_reads, n_classified, n_bases = _worker_file(path)
                total_reads += n_reads
                total_classified += n_classified
                total_bases += n_bases
                pct = 100.0 * total_classified / total_reads if total_reads else 0.0
                pbar.set_postfix(classified=f"{pct:.1f}%", reads=f"{total_reads/1e6:.2f}M")
                pbar.update(1)
        else:
            with Pool(processes=n_workers, initializer=_init_worker, initargs=init_args) as pool:
                for n_reads, n_classified, n_bases in pool.imap_unordered(_worker_file, fastq_paths):
                    total_reads += n_reads
                    total_classified += n_classified
                    total_bases += n_bases
                    pct = 100.0 * total_classified / total_reads if total_reads else 0.0
                    pbar.set_postfix(classified=f"{pct:.1f}%", reads=f"{total_reads/1e6:.2f}M")
                    pbar.update(1)

    pct = 100.0 * total_classified / total_reads if total_reads else 0.0
    print(f"[torchlex] {total_classified:,}/{total_reads:,} reads classified ({pct:.1f}%)")

    t_demux = time.perf_counter()

    # --- merge ---
    _merge(args.output)
    t_merge = time.perf_counter()

    # --- triage ---
    passed, failed = _triage(args.output, args.min_reads)
    print(f"[torchlex] passed: {len(passed)} barcodes  failed: {len(failed)} barcodes → failed_barcodes/")
    if failed:
        for name, n in failed:
            print(f"           {name}  ({n} reads)")

    # --- report ---
    t_end = time.perf_counter()
    elapsed_min = (t_end - t_start) / 60.0
    gbp = total_bases / 1e9
    rate = elapsed_min / gbp if gbp > 0 else float("inf")

    print(f"\n[torchlex] timing")
    print(f"  demux   : {t_demux - t_start:.1f}s")
    print(f"  merge   : {t_merge - t_demux:.1f}s")
    print(f"  triage  : {t_end - t_merge:.1f}s")
    print(f"  total   : {elapsed_min:.2f} min ({t_end - t_start:.1f}s)")
    print(f"\n[torchlex] throughput")
    print(f"  input   : {gbp:.3f} Gbp across {len(fastq_paths)} files")
    print(f"  speed   : {rate:.2f} min/Gbp  ({n_workers} thread(s))")


if __name__ == "__main__":
    main()
