import requests
import sys
import os
import time
import json
import csv
from Bio import SeqIO

# --- Konfiguracja ---
# Zmień nazwę pliku wejściowego, jeśli używasz innego niż w poprzednim kroku.
INPUT_FASTA_FILE = "extended_proteins.fa"
OUTPUT_CSV_FILE = "pfam_domain_matrix.csv"

# Adresy API serwera HMMER
HMMER_SCAN_URL = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan"
HMMER_RESULT_URL_BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/"

# Odstęp czasowy (w sekundach) między kolejnymi zapytaniami o status zadania
POLL_INTERVAL_SECONDS = 5


def submit_hmmscan_job(sequence):
    """
    Wysyła pojedyncze zadanie hmmscan do serwera HMMER.

    Args:
        sequence (str): Sekwencja białkowa.

    Returns:
        str: ID zadania, jeśli zostało pomyślnie wysłane, w przeciwnym razie None.
    """
    try:
        request_payload = {'database': 'pfam', 'input': sequence}
        headers = {"Accept": "application/json"}
        response = requests.post(HMMER_SCAN_URL, json=request_payload, headers=headers)

        # Zgłoś wyjątek dla kodów błędów HTTP (4xx lub 5xx)
        response.raise_for_status()

        return response.json().get('id')
    except requests.exceptions.RequestException as e:
        print(f"BŁĄD SIECIOWY: Nie można wysłać zadania. {e}", file=sys.stderr)
        return None


def get_job_results(job_id):
    """
    Pobiera wyniki zakończonego zadania HMMER, cyklicznie sprawdzając jego status.

    Args:
        job_id (str): ID zadania do sprawdzenia.

    Returns:
        dict: Wyniki w formacie JSON, jeśli zadanie zakończyło się pomyślnie, w przeciwnym razie None.
    """
    result_url = HMMER_RESULT_URL_BASE + job_id
    headers = {"Accept": "application/json"}

    while True:
        try:
            response = requests.get(result_url, headers=headers)

            # 200 OK: Zadanie zakończone, wyniki są gotowe
            if response.status_code == 200:
                print("    > Wyniki gotowe.")
                return response.json()

            # 202 Accepted: Zadanie jest wciąż w toku
            elif response.status_code == 202:
                print(f"    > Zadanie w toku, ponawianie za {POLL_INTERVAL_SECONDS}s...")
                time.sleep(POLL_INTERVAL_SECONDS)

            # Inne kody są traktowane jako błędy
            else:
                response.raise_for_status()

        except requests.exceptions.RequestException as e:
            print(f"BŁĄD SIECIOWY: Nie można pobrać wyników dla zadania {job_id}. {e}", file=sys.stderr)
            return None


def main():
    """Główna funkcja programu."""
    # 1. Sprawdź, czy plik wejściowy istnieje
    if not os.path.exists(INPUT_FASTA_FILE):
        print(f"BŁĄD: Plik wejściowy '{INPUT_FASTA_FILE}' nie został znaleziony.", file=sys.stderr)
        sys.exit(1)

    # 2. Wczytaj sekwencje białkowe z pliku FASTA
    print(f"1. Wczytywanie sekwencji z pliku '{INPUT_FASTA_FILE}'...")
    try:
        protein_records = list(SeqIO.parse(INPUT_FASTA_FILE, "fasta"))
        if not protein_records:
            print(f"OSTRZEŻENIE: Plik '{INPUT_FASTA_FILE}' jest pusty lub nie zawiera poprawnych rekordów FASTA.",
                  file=sys.stderr)
            sys.exit(0)
    except Exception as e:
        print(f"BŁĄD: Nie można przetworzyć pliku FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    # 3. Przetwórz każde białko: wyślij zadanie, pobierz wyniki i zapisz domeny
    print("\n2. Wysyłanie zapytań do serwera HMMER i pobieranie wyników...")
    protein_domains = {}  # Słownik do przechowywania domen: {'id_białka': {'domena1', 'domena2'}}
    all_domains = set()  # Zbiór wszystkich unikalnych domen znalezionych we wszystkich białkach

    for record in protein_records:
        protein_id = record.id
        sequence = str(record.seq)
        print(f"  - Przetwarzanie białka: {protein_id} ({len(sequence)} aa)")

        job_id = submit_hmmscan_job(sequence)
        if not job_id:
            continue  # Przejdź do następnego białka w przypadku błędu

        print(f"    > Zadanie wysłane, ID: {job_id}")

        results_json = get_job_results(job_id)
        print(results_json)
        if not results_json:
            continue  # Przejdź do następnego białka w przypadku błędu

        # Przetwarzanie wyników i wyodrębnianie identyfikatorów domen Pfam
        found_domains = set()
        if 'result' in results_json and results_json['result'] and 'hits' in results_json['result']:
            for hit in results_json['result']['hits']:
                # 'acc' to akcesja Pfam, np. "PF00595.22". Używamy części przed kropką.
                domain_acc = hit['acc'].split('.')[0]
                found_domains.add(domain_acc)

        if found_domains:
            print(f"    > Znaleziono {len(found_domains)} unikalnych domen.")
            protein_domains[protein_id] = found_domains
            all_domains.update(found_domains)  # Dodaj znalezione domeny do globalnego zbioru
        else:
            print("    > Nie znaleziono żadnych domen.")
            protein_domains[protein_id] = set()

    # 4. Zapisz wyniki do pliku CSV
    if not protein_domains:
        print("\nNie zebrano żadnych wyników do zapisania.", file=sys.stderr)
        return

    print(f"\n3. Zapisywanie macierzy obecności/nieobecności domen do pliku '{OUTPUT_CSV_FILE}'...")

    # Sortowanie identyfikatorów dla spójnej kolejności w pliku CSV
    sorted_protein_ids = sorted(protein_domains.keys())
    sorted_domain_ids = sorted(list(all_domains))

    try:
        with open(OUTPUT_CSV_FILE, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)

            # Zapisz nagłówek (pierwsza komórka pusta, potem ID domen)
            header = ['protein_id'] + sorted_domain_ids
            writer.writerow(header)

            # Zapisz dane dla każdego białka
            for protein_id in sorted_protein_ids:
                row = [protein_id]
                domains_for_this_protein = protein_domains.get(protein_id, set())
                for domain_id in sorted_domain_ids:
                    row.append(1 if domain_id in domains_for_this_protein else 0)
                writer.writerow(row)

        print(
            f"\nZakończono. Zapisano wyniki dla {len(sorted_protein_ids)} białek i {len(sorted_domain_ids)} unikalnych domen.")

    except IOError as e:
        print(f"BŁĄD: Nie można zapisać pliku CSV '{OUTPUT_CSV_FILE}': {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
