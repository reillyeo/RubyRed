import os
import sys
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

fasta_file = sys.argv[1]
taxonomy_file = sys.argv[2]
output_dir = sys.argv[3]

# Output directory
os.makedirs(output_dir, exist_ok=True)

# Step 1: Read taxonomy.tsv
feature_to_taxon = {}
with open(taxonomy_file, "r") as tax_file:
    header = next(tax_file)
    for line in tax_file:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            feature_id, taxon = parts[0], parts[1]
            feature_to_taxon[feature_id] = taxon

# Step 2: Group sequences by taxon
taxon_to_seqs = defaultdict(list)
for record in SeqIO.parse(fasta_file, "fasta"):
    feature_id = record.id
    taxon = feature_to_taxon.get(feature_id)
    if taxon:
        taxon_to_seqs[taxon].append(record)

# Step 3: Write each taxon to its own FASTA file
def sanitize_filename(name):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)

for taxon, records in taxon_to_seqs.items():
    filename = sanitize_filename(taxon) + ".fasta"
    output_path = Path(output_dir) / filename
    SeqIO.write(records, output_path, "fasta")

print(f"Created {len(taxon_to_seqs)} FASTA files in '{output_dir}'")
