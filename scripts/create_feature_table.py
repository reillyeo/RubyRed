import sys
from collections import defaultdict

def parse_fasta(fasta_file):
    """
    Parse a FASTA file and group sequences by sample ID.
    Assumes sequence IDs are in the format: SAMPLEID_x (e.g., BRK01_1).
    """
    sample_counts = defaultdict(lambda: defaultdict(int))  # {sample_id: {seq_id: count}}
    current_seq_id = None

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # Header line
                current_seq_id = line[1:]  # Remove '>'
                sample_id = current_seq_id.split("_")[0]  # Extract sample ID (e.g., BRK01 from BRK01_1)
                sample_counts[sample_id][current_seq_id] += 1
            else:
                # Optionally validate sequence content here
                pass

    return sample_counts


def write_feature_table(sample_counts, output_file):
    """
    Write a feature table in TSV format from the parsed sample counts.
    """
    # Collect all unique feature IDs
    all_features = sorted({seq_id for sample in sample_counts.values() for seq_id in sample})

    # Write the header
    with open(output_file, "w") as f:
        f.write("#OTU ID\t" + "\t".join(sorted(sample_counts.keys())) + "\n")  # Header line

        # Write the counts for each feature
        for feature in all_features:
            counts = [sample_counts[sample].get(feature, 0) for sample in sorted(sample_counts.keys())]
            f.write(feature + "\t" + "\t".join(map(str, counts)) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_feature_table.py <input_fasta> <output_tsv>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_tsv = sys.argv[2]

    # Step 1: Parse the FASTA file
    sample_counts = parse_fasta(input_fasta)

    # Step 2: Write the feature table
    write_feature_table(sample_counts, output_tsv)

    print(f"Feature table written to {output_tsv}")
