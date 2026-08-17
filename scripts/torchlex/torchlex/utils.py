import math
import sys


def _errs_tab(n):
    return [10 ** (q / -10) for q in range(n + 1)]


_ERRS_TAB = _errs_tab(128)


def ave_qual(quals, qround=False):
    """Convert ASCII Phred quals to error probabilities, average, convert back to Phred scale."""
    if quals:
        mq = -10 * math.log(sum(_ERRS_TAB[ord(q) - 33] for q in quals) / len(quals), 10)
        return round(mq) if qround else mq
    sys.stderr.write("Warning: zero-length read quality string\n")
    return 0


def fasta_reader_for_barcodes(handle):
    """Yield (title, [line0, line1, line2]) for barcodes FASTA.

    Each entry has exactly 3 sequence lines: spacer / barcode / linker.
    """
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        return

    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, lines
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, lines
