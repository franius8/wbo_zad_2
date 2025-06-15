import argparse
import csv
import os
import sys
import time

import requests
from Bio import SeqIO

# HMMER servers
HMMER_SCAN_URL = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan"
HMMER_RESULT_URL_BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/"

# Request status polling delay (in seconds)
POLL_INTERVAL_SECONDS = 5

def submit_hmmscan_job(sequence):
    """
   Sends a single hmmscan job to the HMMER server.

    Args:
        sequence (str): protein sequence

    Returns:
        str: job ID if succesfully sent, else None
    """
    try:
        request_payload = {'database': 'pfam', 'input': sequence}
        headers = {"Accept": "application/json"}
        response = requests.post(HMMER_SCAN_URL, json=request_payload, headers=headers)

        # Exception for HTTP errors
        response.raise_for_status()

        return response.json().get('id')
    except requests.exceptions.RequestException as e:
        print(f"NETWORK ERROR: could not send job {e}", file=sys.stderr)
        return None


def get_job_results(job_id):
    """
    Get results from HMMER server.

    Args:
        job_id (str): ID of job to check.

    Returns:
        dict: results in JSON format or None if not successful.
    """
    result_url = HMMER_RESULT_URL_BASE + job_id
    headers = {"Accept": "application/json"}

    while True:
        try:
            response = requests.get(result_url, headers=headers)

            # The job finished and the results are ready
            if response.status_code == 200:
                print("    > Results are ready.")
                return response.json()

            # The job is still ongoing
            elif response.status_code == 202:
                print(f"    > Job ongoing, retrying in {POLL_INTERVAL_SECONDS}s...")
                time.sleep(POLL_INTERVAL_SECONDS)

            # Other codes treated as errors
            else:
                response.raise_for_status()

        except requests.exceptions.RequestException as e:
            print(f"NETWORK ERROR: could not retrieve results for job {job_id}. {e}", file=sys.stderr)
            return None


def main():
    parser = argparse.ArgumentParser(
        description="Analyzes protein sequences from a FASTA file, finds Pfam domains using HMMER, and creates a presence/absence matrix."
    )
    parser.add_argument(
        "-i", "--input",
        default="extended_proteins.fa",
        help="Path to the input FASTA file. Default: 'extended_proteins.fa'"
    )
    parser.add_argument(
        "-o", "--output",
        default="pfam_domain_matrix.csv",
        help="Path to the output CSV file. Default: 'pfam_domain_matrix.csv'"
    )
    args = parser.parse_args()

    # Use filenames from the parsed arguments
    input_fasta_file = args.input
    output_csv_file = args.output

    # 1. Check if the input file exists
    if not os.path.exists(input_fasta_file):
        print(f"ERROR: '{input_fasta_file}' not found", file=sys.stderr)
        sys.exit(1)

    # 2. Read sequences from FASTA file
    print(f"1. Reading sequences from file '{input_fasta_file}'...")
    try:
        protein_records = list(SeqIO.parse(input_fasta_file, "fasta"))
        if not protein_records:
            print(f"WARNING: '{input_fasta_file}' file is empty or does not contain any valid FASTA records.",
                  file=sys.stderr)
            sys.exit(0)
    except Exception as e:
        print(f"ERROR could not process FASTA file {e}", file=sys.stderr)
        sys.exit(1)

    # 3. Process each protein
    print("\n2. Sending requests to HMMER server and getting results...")
    protein_domains = {}
    all_domains = set()

    for record in protein_records:
        protein_id = record.id
        sequence = str(record.seq)
        print(f"  - Processing protein: {protein_id} ({len(sequence)} aa)")

        job_id = submit_hmmscan_job(sequence)
        if not job_id:
            continue  # Process next protein if any errors occur

        print(f"    > Job sent, ID: {job_id}")

        results_json = get_job_results(job_id)
        print(results_json)
        if not results_json:
            continue  # Process next protein if any errors occur

        # Processing results
        found_domains = set()
        if 'result' in results_json and results_json['result'] and 'hits' in results_json['result']:
            for hit in results_json['result']['hits']:
                domain_acc = hit['acc'].split('.')[0]
                found_domains.add(domain_acc)

        if found_domains:
            print(f"    > {len(found_domains)} unique domains found.")
            protein_domains[protein_id] = found_domains
            all_domains.update(found_domains)  # Add found domains to a global set
        else:
            print("    > No domains found.")
            protein_domains[protein_id] = set()

    # 4. Save results to a CSV file
    if not protein_domains:
        print("\nNo results to save.", file=sys.stderr)
        return

    print(f"\n3. Saving results to file '{output_csv_file}'...")

    # Sorting IDs in the CSV
    sorted_protein_ids = sorted(protein_domains.keys())
    sorted_domain_ids = sorted(list(all_domains))

    try:
        with open(output_csv_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)

            # Save header
            header = ['protein_id'] + sorted_domain_ids
            writer.writerow(header)

            # Save data for each protein
            for protein_id in sorted_protein_ids:
                row = [protein_id]
                domains_for_this_protein = protein_domains.get(protein_id, set())
                for domain_id in sorted_domain_ids:
                    row.append(1 if domain_id in domains_for_this_protein else 0)
                writer.writerow(row)

        print(
            f"\nSaved results for {len(sorted_protein_ids)} proteins and {len(sorted_domain_ids)} unique domains.")

    except IOError as e:
        print(f"ERROR: could not save'{output_csv_file}': {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
