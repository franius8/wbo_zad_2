import sys
import os
import csv

# --- Konfiguracja ---
# Plik wejściowy z wynikami BLAST, aby zidentyfikować, które geny nas interesują.
# Powinien to być ten sam plik, który został wygenerowany przez pierwszy skrypt.
BLAST_RESULTS_FILE = "blast_results.tsv"

# Plik wejściowy z adnotacjami Gene Ontology (w formacie GAF).
# Musisz pobrać ten plik dla E. coli, np. z http://current.geneontology.org/annotations/
GAF_FILE = "e_coli.gaf"

# Nazwa pliku wyjściowego CSV.
OUTPUT_CSV_FILE = "go_term_matrix.csv"


def get_relevant_gene_ids(blast_file):
    """
    Odczytuje plik wyników BLAST i zwraca zbiór unikalnych identyfikatorów genów (sseqid),
    które miały dopasowanie.
    """
    if not os.path.exists(blast_file):
        print(f"BŁĄD: Plik z wynikami BLAST '{blast_file}' nie został znaleziony.", file=sys.stderr)
        return None

    gene_ids = set()
    print(f"1. Odczytywanie identyfikatorów genów z '{blast_file}'...")
    with open(blast_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # Kolumna sseqid (subject id) jest druga (indeks 1)
            if len(parts) > 1:
                gene_ids.add(parts[1])

    if not gene_ids:
        print("OSTRZEŻENIE: Nie znaleziono żadnych identyfikatorów genów w pliku BLAST.", file=sys.stderr)
    else:
        print(f"   > Znaleziono {len(gene_ids)} unikalnych genów do analizy.")

    return gene_ids


def parse_gaf_file(gaf_file, relevant_genes):
    """
    Przetwarza plik GAF i tworzy mapowanie genów na ich terminy GO.

    Args:
        gaf_file (str): Ścieżka do pliku GAF.
        relevant_genes (set): Zbiór ID genów, które nas interesują.

    Returns:
        tuple: Słownik mapujący ID genu na zbiór jego terminów GO oraz
               zbiór wszystkich unikalnych terminów GO znalezionych dla tych genów.
    """
    if not os.path.exists(gaf_file):
        print(f"BŁĄD: Plik GAF '{gaf_file}' nie został znaleziony.", file=sys.stderr)
        return None, None

    gene_to_go = {}  # Słownik: {'id_genu': {'GO:00001', 'GO:00002'}}
    all_go_terms = set()

    print(f"2. Przetwarzanie pliku adnotacji '{gaf_file}'...")
    with open(gaf_file, 'r') as f:
        for line in f:
            # Pomiń linie komentarza, które zaczynają się od '!'
            if line.startswith('!'):
                continue

            parts = line.strip().split('\t')
            # Format GAF ma co najmniej 17 kolumn. Potrzebujemy:
            # - Kolumna 2 (indeks 1): DB Object Symbol (ID genu)
            # - Kolumna 5 (indeks 4): GO ID
            if len(parts) >= 5:
                gene_id = parts[2]
                go_id = parts[4]

                # Sprawdź, czy gen jest na naszej liście interesujących genów
                if gene_id in relevant_genes:
                    # Dodaj mapowanie, jeśli jeszcze nie istnieje
                    if gene_id not in gene_to_go:
                        gene_to_go[gene_id] = set()

                    gene_to_go[gene_id].add(go_id)
                    all_go_terms.add(go_id)

    print(f"   > Znaleziono adnotacje GO dla {len(gene_to_go)} genów.")
    print(f"   > Łącznie zidentyfikowano {len(all_go_terms)} unikalnych terminów GO.")
    return gene_to_go, all_go_terms


def main():
    """Główna funkcja programu."""
    # Krok 1: Pobierz listę genów, które miały dopasowanie w BLAST
    relevant_gene_ids = get_relevant_gene_ids(BLAST_RESULTS_FILE)
    if relevant_gene_ids is None:
        sys.exit(1)
    if not relevant_gene_ids:
        print("\nZakończono. Brak genów do dalszego przetwarzania.")
        sys.exit(0)

    # Krok 2: Przetwórz plik GAF, aby znaleźć funkcje GO dla tych genów
    gene_go_mappings, all_go_terms = parse_gaf_file(GAF_FILE, relevant_gene_ids)
    if gene_go_mappings is None:
        sys.exit(1)
    if not gene_go_mappings:
        print("\nZakończono. Nie znaleziono żadnych adnotacji GO dla podanych genów.")
        sys.exit(0)

    # Krok 3: Zapisz wyniki do pliku CSV
    print(f"\n3. Zapisywanie macierzy obecności/nieobecności terminów GO do pliku '{OUTPUT_CSV_FILE}'...")

    # Sortowanie identyfikatorów dla spójnej kolejności w pliku CSV
    sorted_gene_ids = sorted(list(gene_go_mappings.keys()))
    sorted_go_ids = sorted(list(all_go_terms))

    try:
        with open(OUTPUT_CSV_FILE, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)

            # Zapisz nagłówek (pierwsza komórka pusta, potem ID terminów GO)
            header = ['gene_id'] + sorted_go_ids
            writer.writerow(header)

            # Zapisz dane dla każdego genu
            for gene_id in sorted_gene_ids:
                row = [gene_id]
                go_terms_for_this_gene = gene_go_mappings.get(gene_id, set())
                for go_id in sorted_go_ids:
                    row.append(1 if go_id in go_terms_for_this_gene else 0)
                writer.writerow(row)

        print(
            f"\nZakończono. Zapisano wyniki dla {len(sorted_gene_ids)} genów i {len(sorted_go_ids)} unikalnych terminów GO.")

    except IOError as e:
        print(f"BŁĄD: Nie można zapisać pliku CSV '{OUTPUT_CSV_FILE}': {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
