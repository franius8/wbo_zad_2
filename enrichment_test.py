import sys
import argparse
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def run_enrichment_analysis(matrix_file, id_list_file, alpha):
    """
    Główna funkcja przeprowadzająca analizę wzbogacenia.

    Args:
        matrix_file (str): Ścieżka do pliku CSV z macierzą cech (domeny/GO).
        id_list_file (str): Ścieżka do pliku z listą ID do testowania.
        alpha (float): Próg istotności statystycznej.
    """
    # --- 1. Wczytywanie danych ---
    print(f"1. Wczytywanie macierzy z pliku: '{matrix_file}'")
    try:
        # Wczytaj macierz, używając pierwszej kolumny jako indeksu wierszy
        df = pd.read_csv(matrix_file, index_col=0)
    except FileNotFoundError:
        print(f"BŁĄD: Nie znaleziono pliku macierzy '{matrix_file}'.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"BŁĄD: Nie można wczytać pliku CSV: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"2. Wczytywanie listy identyfikatorów z pliku: '{id_list_file}'")
    try:
        with open(id_list_file, 'r') as f:
            # Wczytaj ID i usuń puste linie oraz białe znaki
            interesting_ids = {line.strip() for line in f if line.strip()}
    except FileNotFoundError:
        print(f"BŁĄD: Nie znaleziono pliku z listą ID '{id_list_file}'.", file=sys.stderr)
        sys.exit(1)

    print(f"   > Znaleziono {len(interesting_ids)} unikalnych ID w liście.")

    # Sprawdzenie, które z podanych ID faktycznie istnieją w macierzy
    existing_ids = set(df.index).intersection(interesting_ids)
    print(f"   > {len(existing_ids)} z podanych ID istnieje w macierzy.")
    if not existing_ids:
        print("BŁĄD: Żaden z podanych ID nie został znaleziony w pliku macierzy.", file=sys.stderr)
        sys.exit(1)

    # --- 2. Przeprowadzanie testu Fishera dla każdej kolumny ---
    print("\n3. Przeprowadzanie testu Fishera dla każdej cechy...")
    results = []

    # Podział na grupę "interesującą" i tło
    interesting_df = df.loc[list(existing_ids)]
    background_df = df.drop(index=list(existing_ids))

    # Pętla po wszystkich cechach (kolumnach)
    for feature in df.columns:
        # Budowanie tabeli kontyngencji 2x2
        #         | Ma cechę | Nie ma cechy
        # --------|----------|-------------
        # Grupa   |    a     |      b
        # Tło     |    c     |      d

        a = interesting_df[feature].sum()
        b = len(interesting_df) - a
        c = background_df[feature].sum()
        d = len(background_df) - c

        # Pomiń cechy, które nie występują wcale w całym zbiorze
        if a + c == 0:
            continue

        contingency_table = [[a, b], [c, d]]

        # Oblicz test Fishera
        odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')

        results.append({
            'feature': feature,
            'p_value': p_value,
            'group_hits': int(a),
            'group_total': int(a + b),
            'background_hits': int(c),
            'background_total': int(c + d)
        })

    if not results:
        print("Nie przeprowadzono żadnych testów (może to oznaczać brak cech w danych).")
        return

    # --- 3. Poprawka na wielokrotne testowanie (Bonferroni) ---
    print("\n4. Stosowanie poprawki Bonferroniego...")
    raw_p_values = [res['p_value'] for res in results]

    reject, pvals_corrected, _, _ = multipletests(raw_p_values, alpha=alpha, method='bonferroni')

    # Dodaj poprawione p-wartości do wyników
    for i, res in enumerate(results):
        res['p_value_corrected'] = pvals_corrected[i]
        res['is_significant'] = reject[i]

    # --- 4. Wyświetlanie wyników ---
    significant_results = [res for res in results if res['is_significant']]

    print("\n" + "=" * 80)
    print(f"WYNIKI ANALIZY WZBOGACENIA (próg alpha={alpha} z poprawką Bonferroniego)")
    print("=" * 80)

    if not significant_results:
        print("Nie znaleziono statystycznie istotnego wzbogacenia dla żadnej cechy.")
    else:
        # Sortuj wyniki po poprawionej p-wartości
        significant_results.sort(key=lambda x: x['p_value_corrected'])

        print(f"Znaleziono {len(significant_results)} istotnie wzbogaconych cech:\n")

        # Wyświetl nagłówek tabeli
        print(
            f"{'Cecha':<20} | {'P-wartość (popr.)':<20} | {'P-wartość (surowa)':<20} | {'Trafienia w grupie':<20} | {'Trafienia w tle'}")
        print("-" * 110)

        for res in significant_results:
            group_ratio = f"{res['group_hits']}/{res['group_total']}"
            background_ratio = f"{res['background_hits']}/{res['background_total']}"
            print(
                f"{res['feature']:<20} | {res['p_value_corrected']:<20.4e} | {res['p_value']:<20.4e} | {group_ratio:<20} | {background_ratio}")

    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description="Przeprowadza test wzbogacenia Fishera z poprawką Bonferroniego na macierzy cech.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--matrix',
        required=True,
        help="Ścieżka do pliku CSV z macierzą (np. 'pfam_domain_matrix.csv' lub 'go_term_matrix.csv')."
    )
    parser.add_argument(
        '--id_list',
        required=True,
        help="Ścieżka do pliku tekstowego z listą identyfikatorów do analizy (jeden ID w linii)."
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.05,
        help="Próg istotności statystycznej (alfa), domyślnie: 0.05."
    )

    args = parser.parse_args()

    run_enrichment_analysis(args.matrix, args.id_list, args.alpha)


if __name__ == "__main__":
    main()
