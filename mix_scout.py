import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from itertools import combinations
from dataclasses import dataclass
import eutectic as eut
import logging

logging.basicConfig(level=logging.INFO)

@dataclass
class MixtureResult:
    combination: tuple
    melting_points: np.ndarray
    concentrations: np.ndarray
    clearing_points: np.ndarray


def parse_args():
    parser = argparse.ArgumentParser(description='Find optimal mixtures from (many) single components')
    parser.add_argument('-i', '--input', required=True, type=str, help='CSV file with compound data')
    parser.add_argument('-mp', '--max_melt', default=100.0, type=float, help='Maximum acceptable melt temperature (deg C)')
    parser.add_argument('-n', '--max_comps', default=2, type=int, help='Maximum number of components in mixture')
    parser.add_argument('-s', '--search_size', default=50000, type=int, help='Base number of search points per mixture')
    parser.add_argument('--export', type=str, help='Optional CSV filename to export good mixtures')
    return parser.parse_args()


def load_compound_data(csv_file):
    logging.info(f"Reading data from {csv_file}")
    column_names = ['Compound Name', 'Melting Point', 'T_Niso', 'Enthalpy']
    df = pd.read_csv(csv_file, header=None, names=column_names)
    df.set_index('Compound Name', inplace=True)
    return df


def generate_combinations(compounds, max_components):
    names = compounds.index.tolist()
    return [combo for r in range(2, max_components + 1) for combo in combinations(names, r)]


def dynamic_search_size(base_size, num_compounds, min_size=100, max_size=1e6):
    return int(min(max_size, max(min_size, base_size * (num_compounds ** 2))))


def evaluate_combinations(compounds, args):
    good_mixtures = []
    combos = generate_combinations(compounds, args.max_comps)

    for combo in tqdm(combos, desc="Evaluating combinations"):
        mp = compounds.loc[list(combo), 'Melting Point'].values

        enthalpies = compounds.loc[list(combo), 'Enthalpy'].values
        t_clear = compounds.loc[list(combo), 'T_Niso'].values
        search_size = dynamic_search_size(args.search_size, len(combo))

        concs, melts, clears = eut.make_phase_diagram(mp, enthalpies, t_clear, search_size=search_size, plotting=False)

        if np.min(melts) <= args.max_melt:
            good_mixtures.append(MixtureResult(combo, melts, concs, clears))

    return good_mixtures


def summarise_results(mixture_list):
    best_low, best_high, best_wide = None, None, None
    melt_min, clear_max, range_max = float('inf'), float('-inf'), float('-inf')

    for mix in mixture_list:
        melt, clear, span = eut.print_composition(mix.concentrations, mix.melting_points, mix.clearing_points, quiet=True)
        if melt < melt_min:
            best_low, melt_min = mix, melt
        if clear > clear_max:
            best_high, clear_max = mix, clear
        if span > range_max:
            best_wide, range_max = mix, span

    return best_low, best_high, best_wide

def print_table_format(title, mixture):
    melt, clear, span = eut.print_composition(mixture.concentrations, mixture.melting_points, mixture.clearing_points, quiet=True)
    print("\n" + "*" * 25)
    print(f"{title}:")
    print(" | ".join(mixture.combination))
    idx_min = np.argmin(mixture.melting_points)
    comps = mixture.concentrations[idx_min]
    print(" | ".join(f"{c:.3f}%" for c in comps))
    print(f"T(melt) = {melt:.2f} °C")
    print(f"T(clear) = {clear:.2f} °C")
    print(f"dT(lc) = {span:.2f} °C")

def main():
    args = parse_args()
    df = load_compound_data(args.input)
    good_mixtures = evaluate_combinations(df, args)

    logging.info(f"Found {len(good_mixtures)} good mixtures")
    best_low, best_high, best_wide = summarise_results(good_mixtures)


    if best_low:
        print_table_format("Lowest Melting Point", best_low)
    if best_high:
        print_table_format("Highest Clearing Point", best_high)
    if best_wide:
        print_table_format("Widest Enantiotropic Range", best_wide)
        
    if args.export:
        export_data = [{
            'Combination': ", ".join(m.combination),
            'Lowest Melt': np.min(m.melting_points),
            'Highest Clear': np.max(m.clearing_points),
        } for m in good_mixtures]
        pd.DataFrame(export_data).to_csv(args.export, index=False)
        logging.info(f"Exported good mixtures to {args.export}")


if __name__ == "__main__":
    main()