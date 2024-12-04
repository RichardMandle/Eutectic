import argparse

import pandas as pd
import numpy as np
from tqdm import tqdm
from itertools import combinations

import eutectic as eut # get our own module

'''
Little script for running through a big pile of DSC data for PURE COMPOUNDS and generating estimated eutectic mix data

Current status:
well, it works. 

Wants your data as a .csv file like this:

NAME | MELT (deg C) | CLEAR (deg C) | dH_fusion (kJ mol-1)

Problems:
None?

Forseen problems:
    * The number of search-points per mixture (-s, args.search_size) is very poorly implemented; 
    it should probably be dynamic, becasue a 6- component mixture would need many more data points
    than a 2- component one.
'''

def initialize():
    '''
    accept commandline arguments:
    '''
    parser = argparse.ArgumentParser(description='Find optimal mixtures from (many) single components')

    parser.add_argument('-i', '--input', default='', type=str, help='csv file with data arranged as name | T_melt (dC)| T_iso (dC)| dH (kJ mol^-1)')
    parser.add_argument('-mp', '--max_melt', default=100, type=float, help='melting point cutoff (max)')
    parser.add_argument('-n', '--max_comps', default=2, type=int, help='maximum number of components')
    parser.add_argument('-s', '--search_size', default=10000, type=int, help='number of points to evaluate per mixture')

    return parser
    
def load_compound_data(csv_file):
    '''
    load DSC data from a .csv file; has to be in the order name | T_melt | T_iso | enthalpy
    
    Temps in degrees C; enthalpy in kJ mol-1.
    '''
    
    column_names = ['Compound Name', 'Melting Point', 'T_Niso', 'Enthalpy']
    data = pd.read_csv(csv_file, header=None, names = column_names)
    data.set_index('Compound Name', inplace=True)
    return data

def generate_combinations(dataframe, max_components):
    '''
    use itertools to generate all possible combinations from our components list, with maximum set by our input argument
    '''
    compound_names = dataframe.index.tolist()
    all_combinations = []
    for r in range(2,  args.max_comps + 1):
        all_combinations.extend(combinations(compound_names, r))
    return all_combinations

def dynamic_search_size(args, num_compounds, n = 4, min_size = 100, max_size = 1e7):
    '''
    Dynamically calculate a search size based on the number of components
    '''
    
    search_size = max(min_size, min(max_size, args.search_size * (n ** num_compounds)))
    
    return search_size

def evaluate_combinations(dataframe):
    '''
    retrieve the DSC data from the pd.dataframe for the materials in each combination and evaluate
    
    if they are below our max permitted melt, store them in "good mixtures"
    '''
    good_mixtures = []

    combinations = generate_combinations(dataframe, args.max_comps)
    num_combinations = len(combinations)

    data = np.empty((num_combinations, 4), dtype=object) # np array to store results

    for x, combination in enumerate(tqdm(combinations, desc="Evaluating combinations")):
        num_compounds = len(combination)
        
        melting_points = np.zeros(num_compounds)
        clearing_points = np.zeros(num_compounds)
        enthalpies = np.zeros(num_compounds)

        for idx, compound in enumerate(combination):
            #compound_data = dataframe.loc[compound]
            melting_points[idx] = dataframe.loc[compound]['Melting Point']
            clearing_points[idx] = dataframe.loc[compound]['T_Niso']
            enthalpies[idx] = dataframe.loc[compound]['Enthalpy']
        
        search_size = dynamic_search_size(args, num_compounds)
          
        concs, melt, clear = eut.make_phase_diagram(melting_points, clearing_points, enthalpies, search_size=search_size, plotting=False)

        data[x, 0] = combination
        data[x, 1] = melt
        data[x, 2] = concs
        data[x, 3] = clear

    # outsid ethe main loop, iterate over the array and return 'good' mixtres to the list.
    for x in range(num_combinations):
        if np.min(data[x, 1]) <= args.max_melt:
            good_mixtures.append({
                'combination': data[x, 0],
                'melting_points': data[x, 1],
                'concentrations': data[x, 2],
                'clearing_points': data[x, 3]
            })

    return good_mixtures

if __name__ == "__main__":
    '''
    main thing
    '''
    args = initialize().parse_args()

    if args.input:
        compounds = load_compound_data(args.input)
        good_mixtures = evaluate_combinations(compounds)
        melt = 9999 # use some large/small values for melt/high/wide; we use these to track the lowest melting point, highest clearing point and widest phase range in loops later.
        wide = high = -9999
        
        lowest_melt, widest_range =[], []
        print('Found ' + str(len(good_mixtures)) + ' suitable mixtures!')
        
        for mixture in good_mixtures:
            #print(mixture['combination'])
            melting_point, clearing_point, enantiotropic_range = eut.print_composition(mixture['concentrations'],mixture['melting_points'],mixture['clearing_points'], quiet = True)
            
            if melting_point <= melt:
                lowest_melt = mixture
                melt = melting_point
                
            if clearing_point >= high:
                highest_point = mixture
                high = clearing_point
                
            if enantiotropic_range >= wide:
                widest_range = mixture
                wide = enantiotropic_range
                        
        if lowest_melt != []:
            print('*' * 25)    
            print('Lowest Melting Point:')
            print(lowest_melt['combination'])
            _, _, _ = eut.print_composition(lowest_melt['concentrations'],lowest_melt['melting_points'],lowest_melt['clearing_points'])
            print('\n')
            
        if highest_point != []:
            print('*' * 25)    
            print('Highest Clearing Point:')
            print(highest_point['combination'])
            _, _, _ = eut.print_composition(highest_point['concentrations'],highest_point['melting_points'],highest_point['clearing_points'])
            print('\n')   
            
        if widest_range != []:
            print('*' * 25)    
            print('Widest Enantiotropic Range:')
            print(widest_range['combination'])
            _, _, _ = eut.print_composition(widest_range['concentrations'],widest_range['melting_points'],widest_range['clearing_points'])
            print('*' * 25)    
    else:
        print("No input file specified. Please provide a CSV file with compound data.")