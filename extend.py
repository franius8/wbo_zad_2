import os
import sys
import subprocess
from shutil import which
from Bio.Seq import Seq
from Bio import SeqIO

# Input file names
GENES_FASTA_FILE = "genes_e_coli_new.fa"
FRAGMENTS_FASTA_FILE = "protein_fragments.fa"
OUTPUT_FASTA_FILE = "extended_proteins.fa"

# Output file names
BLAST_DB_NAME = "ecoli_db"
BLAST_RESULTS_FILE = "blast_results.tsv"

# Check if provided files exist
if not os.path.exists(GENES_FASTA_FILE):
    print(f"BŁĄD: Plik wejściowy z genami '{GENES_FASTA_FILE}' nie został znaleziony.", file=sys.stderr)
    sys.exit(1)

if not os.path.exists(FRAGMENTS_FASTA_FILE):
    print(f"BŁĄD: Plik wejściowy z fragmentami '{FRAGMENTS_FASTA_FILE}' nie został znaleziony.", file=sys.stderr)
    sys.exit(1)

print(f"1. Tworzenie bazy danych BLAST z pliku '{GENES_FASTA_FILE}'...")
db_name = "ecoli_db"
makeblastdb_cmd = ["makeblastdb", "-in", GENES_FASTA_FILE, "-dbtype", "nucl", "-out", db_name]
subprocess.run(makeblastdb_cmd, check=True, capture_output=True)

print(f"2. Przeszukiwanie bazy danych przy użyciu '{FRAGMENTS_FASTA_FILE}'...")
blast_results_file = "blast_results.tsv"
tblastn_cmd = [
    "tblastn",
    "-query", FRAGMENTS_FASTA_FILE,
    "-db", db_name,
    "-out", blast_results_file,
    "-outfmt", "6 qseqid sseqid sframe pident length evalue bitscore",
    "-max_target_seqs", "1"
]
subprocess.run(tblastn_cmd, check=True, capture_output=True)

print(f"3. Przetwarzanie wyników i zapisywanie do '{OUTPUT_FASTA_FILE}'...")
gene_records = SeqIO.to_dict(SeqIO.parse(GENES_FASTA_FILE, "fasta"))
output_records = []

with open(blast_results_file, 'r') as results:
    for line in results:
        parts = line.strip().split('\t')
        if len(parts) < 7: continue  # Pomijaj niepoprawnie sformatowane linie

        # Rozpakuj wszystkie kolumny, nawet jeśli nie wszystkie są używane
        query_id, subject_id, frame_str, pident, length, evalue, bitscore = parts
        frame = int(frame_str)

        # Znajdź odpowiedni gen
        if subject_id in gene_records:
            gene_record = gene_records[subject_id]
            dna_seq = gene_record.seq

            # Użyj odpowiedniej nici DNA i dokonaj translacji
            if frame < 0:
                dna_seq = dna_seq.reverse_complement()

            start_pos = abs(frame) - 1
            protein_seq = dna_seq[start_pos:].translate(to_stop=True)

            # Przygotuj rekord do zapisu, dodając więcej informacji do nagłówka
            output_record = SeqIO.SeqRecord(
                seq=protein_seq,
                id=query_id,
                description=f"| matched_gene:{subject_id} frame:{frame} pident:{pident} eval:{evalue}"
            )
            output_records.append(output_record)

if output_records:
    SeqIO.write(output_records, OUTPUT_FASTA_FILE, "fasta")
    print(f"\nZakończono. Zapisano {len(output_records)} sekwencji do pliku '{OUTPUT_FASTA_FILE}'.")
else:
    print("\nZakończono. Nie znaleziono żadnych dopasowań.")