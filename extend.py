import argparse
import os
import subprocess
import sys

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Extend protein fragments by finding matching gene sequences using BLAST.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    "--genes",
    default="genes_e_coli_new.fa",
    help="Path to the input FASTA file containing gene sequences."
)
parser.add_argument(
    "--fragments",
    default="protein_fragments.fa",
    help="Path to the input FASTA file containing protein fragment sequences."
)
parser.add_argument(
    "--output",
    default="extended_proteins.fa",
    help="Path to the output FASTA file for extended protein sequences."
)
parser.add_argument(
    "--db_name",
    default="ecoli_db",
    help="Name for the BLAST database to be created."
)
parser.add_argument(
    "--blast_results",
    default="blast_results.tsv",
    help="Path to the output file for BLAST results."
)

args = parser.parse_args()

# Fasta filenames from arguments
genes_fasta_file = args.genes
fragments_fasta_file = args.fragments
output_fasta_file = args.output

# BLAST filenames from arguments
blast_db_name = args.db_name
blast_results_file = args.blast_results

# Check if provided files exist
if not os.path.exists(genes_fasta_file):
    print(f"ERROR: input file '{genes_fasta_file}' not found.", file=sys.stderr)
    sys.exit(1)

if not os.path.exists(fragments_fasta_file):
    print(f"ERROR: input file '{fragments_fasta_file}' not found.", file=sys.stderr)
    sys.exit(1)

# Create BLAST database from the input file
print(f"1. Creating BLAST database from file '{genes_fasta_file}'...")
makeblastdb_cmd = ["makeblastdb", "-in", genes_fasta_file, "-dbtype", "nucl", "-out", blast_db_name]
subprocess.run(makeblastdb_cmd, check=True, capture_output=True)

# Search the created BLAST database using the input fragments file
print(f"2. Searching created BLAST database using '{fragments_fasta_file}'...")
blast_results_file = "blast_results.tsv"
tblastn_cmd = [
    "tblastn",
    "-query", fragments_fasta_file,
    "-db", blast_db_name,
    "-out", blast_results_file,
    "-outfmt", "6 qseqid sseqid sframe pident length evalue bitscore",
    "-max_target_seqs", "1"
]
subprocess.run(tblastn_cmd, check=True, capture_output=True)

# Process the results and save them to the output file.
print(f"3. Processing results and saving to '{output_fasta_file}'...")
gene_records = SeqIO.to_dict(SeqIO.parse(genes_fasta_file, "fasta"))
output_records = []

with open(blast_results_file, 'r') as results:
    for line in results:
        parts = line.strip().split('\t')
        if len(parts) < 7: continue  # Skip malformed lines

        # Unpack row
        query_id, subject_id, frame_str, pident, length, evalue, bitscore = parts
        frame = int(frame_str)

        # Find the right gene
        if subject_id in gene_records:
            gene_record = gene_records[subject_id]
            dna_seq = gene_record.seq

            # Use the right dna strand and perform translation
            if frame < 0:
                dna_seq = dna_seq.reverse_complement()

            start_pos = abs(frame) - 1
            protein_seq = dna_seq[start_pos:].translate(to_stop=True)

            # Prepare record for saving
            output_record = SeqIO.SeqRecord(
                seq=protein_seq,
                id=query_id,
                description=f"| matched_gene:{subject_id} frame:{frame} pident:{pident} eval:{evalue}"
            )
            output_records.append(output_record)

if output_records:
    SeqIO.write(output_records, output_fasta_file, "fasta")
    print(f"\nSaved {len(output_records)} sequences to file '{output_fasta_file}'.")
else:
    print("\nNo matches found.")