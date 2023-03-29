from npktools import run_all, convert_copy_number_to_concentration, is_nucleic_acid

default_inputs = {
    'strands': None,
    'nucleic_acid_type': 'DNA',
    'temperature_C': [37],
    'complexes': 2,
    'concentration_dict': None,
    'Na_M': [0.150],
    'Mg_M': [0.010],
    'species': None,
    'sampleVolume_uL': 25,
}

user_strands = {'A': 'CCTAGCTCTGACCACTTCACACCTACGACCACAGATGGTGAGGATGGCTGACCTGACGTACGGCTCTCA',
                'B': 'GTTCACCTTCAAGAGTTTCTTCCTATGAGAGCCGTACGTCAGGTC'}
# make sure that the keys are strings and all the values are nucleic acid sequences
assert all([isinstance(key, str) and is_nucleic_acid(value) for key, value in user_strands.items()])

# assert all([isinstance(key, str) for key in user_strands.keys()])
# assert all([is_nucleic_acid(value, default_inputs['nucleic_acid_type']) for value in user_strands.values()])


user_concentration_dict = {'A': [1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13],
                           'B': [convert_copy_number_to_concentration(25, 1e7)]}

user_species = 'B'
assert user_species in user_strands

default_inputs['strands'] = user_strands

default_inputs['concentration_dict'] = user_concentration_dict

default_inputs['species'] = user_species




if __name__ == "__main__":
    # Run the analysis
    results, df_Nupack, temp_concBound, t_result, c_result = run_all(**default_inputs, savefig=True)
