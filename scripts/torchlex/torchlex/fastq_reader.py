import gzip
import io
import os

from Bio.SeqIO.QualityIO import FastqGeneralIterator


class FastqFile:
    """Read a FASTQ file into a list of FastqRow objects.

    strict_filenames: enforce Oxford Nanopore naming {flow_cell_id}_pass_{runid}_{idx}.fastq
    """

    def __init__(self, path, strict_filenames=False):
        self.path = path
        self.filename = os.path.basename(path)
        # Strip .gz so the stem is the same whether compressed or not
        stem = self.filename
        if stem.endswith(".gz"):
            stem = stem[:-3]
        self.filename_without_extension = os.path.splitext(stem)[0]
        self.strict_filenames = strict_filenames
        self.rows = []

        self._read_file()

        if self.rows:
            self.run_id = self.rows[0].attributes.get("runid", "unknown")
            self.flow_cell_id = self.rows[0].attributes.get("flow_cell_id", "unknown")
        else:
            self.run_id = "unknown"
            self.flow_cell_id = "unknown"

    def _read_file(self):
        if self.path.endswith(".gz"):
            opener = lambda: io.TextIOWrapper(gzip.open(self.path, "rb"), encoding="utf-8")
        else:
            opener = lambda: open(self.path, "r")
        with opener() as f:
            for title, seq, qual in FastqGeneralIterator(f):
                attrs = self._parse_title_line(title)
                if self.strict_filenames:
                    self._check_filename(self.filename, attrs)
                self.rows.append(FastqRow(attrs, seq, qual))

    @staticmethod
    def _parse_title_line(title_line):
        attrs = {}
        for part in title_line.split(" "):
            try:
                key, val = part.split("=", 1)
                attrs[key] = val
            except ValueError:
                attrs["seq_id"] = part
        return attrs

    @staticmethod
    def _check_filename(filename, attrs):
        parts = filename.split("_")
        if len(parts) < 3:
            raise ValueError(
                f"Filename '{filename}' does not match "
                "{flow_cell_id}_pass_{runid}_{idx}.fastq"
            )
        expected_runid = parts[2]
        expected_flow_cell_id = parts[0]
        if not attrs.get("runid", "").startswith(expected_runid):
            raise ValueError(
                f"runid mismatch in '{filename}': "
                f"header has '{attrs.get('runid')}', filename implies '{expected_runid}'"
            )
        if attrs.get("flow_cell_id", "") != expected_flow_cell_id:
            raise ValueError(
                f"flow_cell_id mismatch in '{filename}': "
                f"header has '{attrs.get('flow_cell_id')}', filename implies '{expected_flow_cell_id}'"
            )

    def __len__(self):
        return len(self.rows)


class FastqRow:
    def __init__(self, attributes, seq, qual):
        self.attributes = attributes
        self.seq = seq
        self.qual = qual

    def get_attribute(self, key):
        return self.attributes[key]


class ClassifiedFastqRow(FastqRow):
    def __init__(self, barcode, row):
        self.barcode = barcode
        self.attributes = row.attributes
        self.seq = row.seq
        self.qual = row.qual


class ClassifiedFastqFile:
    """Container for demultiplexed reads from a single FASTQ file."""

    def __init__(self, source_file):
        self.path = source_file.path
        self.filename = source_file.filename
        self.filename_without_extension = source_file.filename_without_extension
        self.run_id = source_file.run_id
        self.flow_cell_id = source_file.flow_cell_id
        self.rows = []
        self._fasta_rows = []

    def add_row(self, row):
        self.rows.append(row)
        seq_id = row.attributes.get("seq_id", "unknown")
        self._fasta_rows.append(f">{row.barcode}_{seq_id}\n{row.seq}\n")

    def get_fasta_string(self):
        return "".join(self._fasta_rows)

    def __len__(self):
        return len(self.rows)
