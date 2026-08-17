import os

import numpy as np
from Bio.Seq import Seq
from polyleven import levenshtein

from torchlex.fastq_reader import (
    ClassifiedFastqFile,
    ClassifiedFastqRow,
    FastqFile,
    FastqRow,
)
from torchlex.utils import ave_qual, fasta_reader_for_barcodes


class Barcodes:
    """Load and index the 96-barcode reference FASTA.

    Each FASTA entry has 3 sequence lines: spacer (15 bp), barcode (24 bp), linker (15 bp).
    Entries suffixed _1st represent forward-strand reads; _2nd entries carry the same barcode
    but spacers/linkers are stored as-written so this class can reverse-complement them when
    loading (matching the orientation seen at the 3' end of a read).
    """

    def __init__(self, barcodes_path, barcode_length=24, spacer_length=15, linker_length=15):
        self.barcodes_path = barcodes_path
        self.barcode_length = barcode_length
        self.spacer_length = spacer_length
        self.linker_length = linker_length

        self._barcodes_fwd = {}
        self._barcodes_rev = {}
        self._spacers_fwd = set()
        self._spacers_rev = set()
        self._linkers_fwd = set()
        self._linkers_rev = set()

        self._load()

    def _load(self):
        with open(self.barcodes_path, "r") as f:
            for title, seq in fasta_reader_for_barcodes(f):
                if title.endswith("_1st"):
                    name = title[:-4]
                    self._barcodes_fwd[name] = seq[1]
                    self._spacers_fwd.add(seq[0])
                    self._linkers_fwd.add(seq[2])
                elif title.endswith("_2nd"):
                    name = title[:-4]
                    self._barcodes_rev[name] = seq[1]
                    self._spacers_rev.add(str(Seq(seq[0]).reverse_complement()))
                    self._linkers_rev.add(str(Seq(seq[2]).reverse_complement()))

    def get_barcodes(self, rev_complemented=False):
        return self._barcodes_rev if rev_complemented else self._barcodes_fwd

    def get_spacers(self, rev_complemented=False):
        return self._spacers_rev if rev_complemented else self._spacers_fwd

    def get_linkers(self, rev_complemented=False):
        return self._linkers_rev if rev_complemented else self._linkers_fwd


class Demultiplexer:
    """Assign each FASTQ read to one of up to 96 barcoded samples.

    Algorithm:
    1. Filter reads by average Phred quality and length.
    2. Extract a window from each end of the read.
    3. In each window, slide over the sequence to find the closest spacer and linker anchors.
    4. Use the anchor positions to extract the candidate 24-bp barcode.
    5. Match the candidate barcode against all 96 references by Levenshtein distance.
    6. Accept the assignment if the best distance is below max_distance (default 5).
    7. If no assignment below threshold, try again on the reverse complement.

    Read structure expected by the barcode library:
      SPACER(15bp) – BARCODE(24bp) – LINKER(15bp) – …amplicon DNA… – LINKER_rc – BARCODE_rc – SPACER_rc
    """

    def __init__(self, config):
        seq = config["sequence"]
        amp = config["amplicones"]

        self.min_qual = seq["minimum_average_phred_quality"]
        self.min_len = seq["min_sequence_length"]
        self.max_len = seq["max_sequence_length"]
        self.headcrop = seq["sequence_headcrop"]
        self.tailcrop = seq["sequence_tailcrop"]
        self.window = seq["barcode_window_size"]
        self.max_distance = config.get("max_barcode_distance", 5)

        self.barcodes = Barcodes(
            barcodes_path=amp["barcodes_path"],
            barcode_length=amp["barcode_length"],
            spacer_length=amp["spacer_length"],
            linker_length=amp["linker_length"],
        )

    # ------------------------------------------------------------------
    # Quality / length filter
    # ------------------------------------------------------------------

    def _filter(self, seq, qual):
        min_len = self.min_len + self.headcrop + self.tailcrop
        if ave_qual(qual) < self.min_qual or not (min_len <= len(seq) <= self.max_len):
            return None, None
        end = len(seq) - self.tailcrop if self.tailcrop > 0 else len(seq)
        return seq[self.headcrop:end], qual[self.headcrop:end]

    # ------------------------------------------------------------------
    # Anchor search (spacer / linker)
    # ------------------------------------------------------------------

    def _find_spacer(self, part, rev_complemented):
        """Slide a window over `part` to locate the nearest spacer; return (barcode_seq, anchor_dist)."""
        spacers = self.barcodes.get_spacers(rev_complemented)
        slen = self.barcodes.spacer_length
        blen = self.barcodes.barcode_length

        best_dist = np.inf
        best_idx = -1
        for spacer in spacers:
            for i in range(len(part) - slen):
                d = levenshtein(part[i:i + slen], spacer)
                if d < best_dist:
                    best_dist = d
                    best_idx = i

        # Spacer precedes barcode (fwd) or follows barcode (rev)
        barcode_start = (best_idx + slen) if not rev_complemented else (best_idx - blen)
        barcode = part[barcode_start:barcode_start + blen]
        return barcode, best_dist

    def _find_linker(self, part, rev_complemented):
        """Slide a window over `part` to locate the nearest linker; return (barcode_seq, anchor_dist)."""
        linkers = self.barcodes.get_linkers(rev_complemented)
        llen = self.barcodes.linker_length
        blen = self.barcodes.barcode_length

        best_dist = np.inf
        best_idx = -1
        for linker in linkers:
            for i in range(len(part) - llen):
                d = levenshtein(part[i:i + llen], linker)
                if d < best_dist:
                    best_dist = d
                    best_idx = i

        # Linker follows barcode (fwd) or precedes barcode (rev)
        barcode_start = (best_idx - blen) if not rev_complemented else (best_idx + llen)
        barcode = part[barcode_start:barcode_start + blen]
        return barcode, best_dist

    # ------------------------------------------------------------------
    # Barcode matching
    # ------------------------------------------------------------------

    def _classify_part(self, part, rev_complemented):
        """Return (barcode_name, distance) for the closest reference barcode in this window."""
        refs = self.barcodes.get_barcodes(rev_complemented)

        bc_spacer, _ = self._find_spacer(part, rev_complemented)
        bc_linker, _ = self._find_linker(part, rev_complemented)

        best_name, best_dist = None, np.inf
        for name, ref in refs.items():
            d = min(levenshtein(ref, bc_spacer), levenshtein(ref, bc_linker))
            if d < best_dist:
                best_name, best_dist = name, d
        return best_name, best_dist

    def _classify_sequence(self, seq):
        """Check both ends of the read; return the assignment with lower distance."""
        start_part = seq[:self.window]
        end_part = seq[-self.window:]

        bc_start, d_start = self._classify_part(start_part, rev_complemented=False)
        bc_end, d_end = self._classify_part(end_part, rev_complemented=True)

        if d_end < d_start:
            return bc_end, d_end, "end"
        return bc_start, d_start, "start"

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def classify_read(self, row: FastqRow):
        """Return (barcode, distance, location, cropped_seq) or (None, None, None, None) if filtered/unassigned."""
        seq, qual = self._filter(row.seq, row.qual)
        if seq is None:
            return None, None, None, None

        bc, dist, loc = self._classify_sequence(seq)

        # Try reverse complement if initial distance is above threshold
        if dist >= self.max_distance:
            seq_rc = str(Seq(seq).reverse_complement())
            bc_rc, dist_rc, loc_rc = self._classify_sequence(seq_rc)
            if dist_rc < dist:
                bc, dist, loc = bc_rc, dist_rc, loc_rc

        if dist < self.max_distance:
            return bc, dist, loc, seq
        return None, None, None, seq

    def demultiplex_file(self, fastq_file: FastqFile, output_dir: str) -> ClassifiedFastqFile:
        """Demultiplex all reads in `fastq_file`, appending each read to
        `output_dir/{barcode}.fastq` in standard 4-line FASTQ format.

        Append mode lets multiple input files accumulate into the same
        per-barcode output files within one run.

        Returns a ClassifiedFastqFile containing only successfully assigned reads.
        """
        classified = ClassifiedFastqFile(fastq_file)
        os.makedirs(output_dir, exist_ok=True)

        handles = {}
        try:
            for row in fastq_file.rows:
                bc, dist, loc, cropped = self.classify_read(row)
                if bc is not None:
                    classified.add_row(ClassifiedFastqRow(bc, row))
                    if bc not in handles:
                        handles[bc] = open(
                            os.path.join(output_dir, f"{bc}.fastq"), "a"
                        )
                    seq_id = row.attributes.get("seq_id", "unknown")
                    handles[bc].write(
                        f"@{seq_id} barcode={bc}\n{row.seq}\n+\n{row.qual}\n"
                    )
        finally:
            for h in handles.values():
                h.close()

        return classified
