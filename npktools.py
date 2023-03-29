import itertools
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
from itertools import cycle
from datetime import datetime
from nupack import *
config.parallelism = True
config.cache = 8.0 # GB

'''
These are Josef's notes:

nupack-4.0.0.26

Nupack download: http://www.nupack.org/downloads

Nupack documentation: https://piercelab-caltech.github.io/nupack-docs/analysis/

Note: As of 14-Feb-2022, the NUPACK web app and the NUPACK python library run on different backend models. 
The python model (herein) is the newest version.

If you want the results of this code to match the NUPACK web app, use the specific model-material designators <br> 'dna04-nupack3' or 'rna95-nupack3'.  
See below.

(NOTE: I suggest using ctr+f to find the correct lines to change. Ex: ctrl + F -> search -> Model(

For DNA <br>
my_model = Model(material='dna04-nupack3', ensemble='some-nupack3', celsius=37)<br>

For RNA:<br>
model1 = Model(material='rna95-nupack3', ensemble='some-nupack3', celsius=37)


# specify strands
A = Strand('CTGATCGAT', name='Strand A')
B = Strand('GATCGTAGTC', name='Strand B')

# specify tube
t1 = Tube(strands={A: 1e-8, B: 1e-9},
    complexes=SetSpec(max_size=2), name='Tube 1') # all complexes of up to 2 strands

# run tube analysis job
my_results = tube_analysis(tubes=[t1], model=model1)
my_results
'''

#%% Supporting functions
def make_list(concentration_A):
    """This puts concentration_A in a list.
       Input: Numerical input. Can be list of/or int/float
       Output: list of concentration_A
    """
    return concentration_A if isinstance(concentration_A, list) else [concentration_A]


def is_string(nucleic_acid_type):
    """This puts nucleic_acid_type and checks it as a string. CODE 1
       Input: String input. Can be 'rna' or 'dna'
       Output: output is input string
    """
    if not isinstance(nucleic_acid_type, str):
        nucleic_acid_type = str(nucleic_acid_type).lower()
    return nucleic_acid_type


def print_strands(A, B):
    """This prints two strands of nucleotides.
       Input: Two string inputs representing nucleotide sequences.
       Output: Two printed strings showing the 5' and 3' ends of each strand.
    """
    print(f"5'-{A}-3'")
    print(f"3'-{B}-5'")


# unused function
def is_nucleic_acid(nucleic_acid_string):
    """Validate a string of nucleobases.

    Args:
        nucleic_acid_string (str): A string of nucleobases (A, T, G, C, U).

    Raises:
        ValueError: If the input is not a string or contains non-nucleobase characters.
    """
    if not isinstance(nucleic_acid_string, str):
        raise ValueError("Nucleic acid input must be a string.")

    valid_nucleobases = {"A", "T", "G", "C", "U"}
    if not all(nucleobase.upper() in valid_nucleobases for nucleobase in nucleic_acid_string):
        raise ValueError(f"Invalid nucleobases in input string: {nucleic_acid_string}")

    logging.info("Input nucleic acid: %s", nucleic_acid_string)


def report_parameters(results):
    """This reports the number of conditions and returns the keys of a dictionary.
       Input: A dictionary `results` with condition names as keys.
       Output: Prints the number of conditions and returns a list of the condition names.
    """
    print(f"{len(results)} number of conditions")
    return list(results.keys())
    # print(str(len(results.keys())) + ' number of conditions')
    # index = list(results.keys())
    # return index


def report_summary(results, index_number):
    """This reports a summary of a specific condition based on its index and returns its corresponding value from a dictionary.
       Input: A dictionary `results` with condition names as keys, and an integer `index_number` representing the index of the condition to report.
       Output: Prints a summary of the specified condition and returns its corresponding value from the dictionary.
    """
    index = report_parameters(results)
    condition = index[index_number]
    print(f"Summary of condition '{condition}': {results[condition]}")
    return results[condition]


def convertK_to_C(tempK):
    """Converts temperature in Kelvin to Celsius and returns the result rounded to the nearest integer.
       Input: A temperature `tempK` in Kelvin.
       Output: The converted temperature `tempC` in Celsius rounded to the nearest integer.
    """
    return round(tempK - 273.15)


def print_tempC(temperature_C):
    """Reports the number of temperatures being evaluated.
       Input: A list of temperature values in Celsius.
       Output: Prints the number of temperature values being evaluated.
    """
    print(f"Number of temperatures being evaluated: {len(temperature_C)}")


def print_tubes(tubes):
    """Reports the number of tube conditions being evaluated.
       Input: A list of tube conditions.
       Output: Prints the number of tube conditions being evaluated.
    """
    num_tubes = len(make_list(tubes)[0])
    print(f"Number of tubes being evaluated: {num_tubes}")


def convert_copy_number_to_concentration(sample_volume_uL=25.0, target_copy_number=1000):
    """Calculates the molar concentration of a target molecule in solution
    given its estimated copy number and the volume of the solution.

    Args:
        sample_volume: the volume of the sample in liters (default: 25 μL)
        target_copy_number: the estimated number of target molecules in the sample (default: 1000)

    Returns:
        The molar concentration of the target molecule in the sample.
    """
    avogadro_number = 6.022e23  # Avogadros number (mol^-1)
    sample_volume_L = sample_volume_uL * 1e-6  # Convert sample size to L
    target_concentration = (target_copy_number / sample_volume_L) / avogadro_number  # convert concentration to molarity
    return float('{:.3g}'.format(target_concentration))


def convert_concentration_to_copy_number(concentration_M, sample_size_uL=25.0):
    """
    Finds copy number from molarity of target in solution
        Input: sampleSize in μL, estimated target molecule amount
        Output: copy number
    """
    avogadro_number = 6.022e23  # Avogadros number (mol^-1)
    sample_size_L = sample_size_uL * 1e-06  # Convert sample size to L
    copy_number = concentration_M * sample_size_L * avogadro_number
    return copy_number
    # return round(float('{:.1g}'.format(copy_number)), 2)


def convert_rna_to_dna(RNA_string='CAUACuAUCaU', material='dna'):
    """
    This function converts RNA strings to DNA strings
    """
    if material.lower() == 'dna':
        DNA_string = RNA_string.upper().replace('U', 'T')
        # DNA_string = ''
        # if "U" in RNA_string or "u" in RNA_string:
        #     DNA_string = (RNA_string.replace('U', 'T').replace('u', 'T')).upper()
    else:
        DNA_string = RNA_string

    return DNA_string


def strand_names_catalog(num_strands=2):
    """Returns a list of strand names based on the input number of strands.
       Input: An integer representing the number of strands.
       Output: A list of strand names.
    """
    try:
        alphabet = [chr(i) for i in range(ord('A'), ord('Z')+1)]  # generate list of letters A-Z
        return [alphabet[i % 26] if i < 26 else alphabet[(i // 26) - 1] + alphabet[i % 26] for i in range(num_strands)]
    except IndexError as e:
        print(f"Error: {e}")
        return None
    # # strand_names = strand_names_catalog()
    # # 104 strand input names. Hopefully no one will input more 104 than strands
    # strand_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
    #                 'U', 'V', 'W', 'X', 'Y', 'Z',
    #                 'AA', 'AB', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AJ', 'AK', 'AL', 'AM', 'AN', 'AO', 'AP',
    #                 'AQ', 'AR', 'AS', 'AT', 'AU', 'AV', 'AW', 'AX', 'AY', 'AZ',
    #                 'BA', 'BB', 'BC', 'BD', 'BE', 'BF', 'BG', 'BH', 'BI', 'BJ', 'BK', 'BL', 'BM', 'BN', 'BO', 'BP',
    #                 'BQ', 'BR', 'BS', 'BT', 'BU', 'BV', 'BW', 'BX', 'BY', 'BZ',
    #                 'CA', 'CB', 'CC', 'CD', 'CE', 'CF', 'CG', 'CH', 'CI', 'CJ', 'CK', 'CL', 'CM', 'CN', 'CO', 'CP',
    #                 'CQ', 'CR', 'CS', 'CT', 'CU', 'CV', 'CW', 'CX', 'CY', 'CZ']
    # return strand_names


def save_fig(fig, name='figures/NUPACK_results_'):
    """ Saviing the figure
        Input: fig object and figure name
        Output: Saved png figure
    """
    date = datetime.now()
    date_str = str(date).replace('-', '').replace(' ', '_').replace(':', '')[:-7]
    fig.savefig(name + date_str + '.png', dpi=100, bbox_inches='tight')
    return print("Figure saved")

#%% Main Functions

def models_maker(nucleic_acid_type='dna', temperature_C=37, Na_M=0.05, Mg_M=0):
    """Input: material: (DNA as 'dna' or RNA as 'rna'),
              temperature (C): temperatures to evaluate as list
              sodium [Na+]: in Molarity (M), (Ex: 1e-9) aka 1nM. Minimum = 0.05
              magnesium [Mg++]: in Molarity (M), (Ex: .000001) aka μM
              Note: This function is set up for only RNA salt concentrations.
              DNA would also need potassium and ammonium. Not sure how to get that functionality

       Output: Models for each temperature.
    """
    if not isinstance(Na_M, list):
        if Na_M < 0.05:
            print('Minimum [Na] >= 0.05')
            print('Resetting Na concentration to 0.05 M')
            Na_M = 0.05
        else:
            pass
    elif isinstance(Na_M, list):
        for C_Na in Na_M:
            if C_Na < 0.05:
                print('Minimum [Na] >= 0.05')
                print('Reseting [Na]: ' + str(C_Na) + ' to 0.05 M')
                index = Na_M.index(C_Na)
                Na_M[index] = 0.05
        else:
            pass

    else:
        pass
        if not isinstance(Mg_M, list):
            if Mg_M < 0.05:
                print('Minimum [Mg] >= 0.05')
                print('Resetting Mg concentration to 0.05 M')
                Mg_M = 0.05
            else:
                pass
        elif isinstance(Mg_M, list):
            for C_Mg in Mg_M:
                if C_Mg < 0.05:
                    print('Minimum [Mg] >= 0.05')
                    print('Reseting [Mg]: ' + C_Mg + ' to 0.05 M')
                    index = Mg_M.index('C_Mg')
                    Mg_M[index] = 0.05
            else:
                pass

        else:
            pass

    oligonucleotide = is_string(nucleic_acid_type).lower()
    models = {}  # collection of models based on temperature
    i = 0
    for temp in make_list(temperature_C):
        for conc_Na in make_list(Na_M):
            for conc_Mg in make_list(Mg_M):
                label = str(temp) + 'C, ' + 'Na+: ' + str(conc_Na) + 'M, ' + ' Mg++: ' + str(conc_Mg) + 'M, '

                #####  II        IIIII    IIIII   II III     II  II IIIIII IIIIII IIIIII    II
                #####  II       II   II  II   II  IIII       IIIIII II---  II__I  II---     II
                #####  IIIIII    IIIII    IIIII   II III     II  II IIIIII II  II IIIIII    ..

                #####  III    III   IIIII   III   IIIII  II       IIIIIIII  III  II  IIIIII  II  II  IIIIIIII
                #####  II I  I II  II   II  II I  II---  II          II     II I II  II__II  II  II     II
                #####  II  I   II   IIIII   III   IIIII  IIIIII   IIIIIIII  II  III  II      IIIIII     II

                # The model below is the must up-to-date
                # models[label] = Model(material=oligonucleotide, celsius=temp, sodium=conc_Na, magnesium=conc_Mg) # salt concentrations in M

                # If you want your output to match the NUPACK web app, comment out the model line about and uncomment
                # the line below.
                # FOR DNA
                models[label] = Model(material='dna04-nupack3', ensemble='some-nupack3', celsius=temp, sodium=conc_Na,
                                      magnesium=conc_Mg)
                # FOR RNA
                # models[label] = Model(material='rna95-nupack3', celsius=temp, sodium=conc_Na, magnesium=conc_Mg)
                i += 1

    # check if list
    if not isinstance(temperature_C, list):
        temperature_C = [temperature_C]
    elif isinstance(temperature_C, list):
        pass
    else:
        pass
    return models


def specify_strands(fiveprime_strands):
    """fiveprime_strand: Needs to be string of A's, T's, G's, C's:   (CODE 3)
    """
    if not isinstance(fiveprime_strands, list):  ######################### NEED TO UPDATE FOR MULTIPLE STRANDS ###
        fiveprime_strands = [fiveprime_strands]

    if len(fiveprime_strands) == 1:
        fiveprime_strand = fiveprime_strands[0]

        threeprime_strand = ''
        if "U" not in fiveprime_strand and "T" in fiveprime_strand:
            for nb in fiveprime_strand:

                if nb.upper() == 'A':
                    threeprime_strand = threeprime_strand + 'T'
                elif nb.upper() == 'T':
                    threeprime_strand = threeprime_strand + 'A'
                elif nb.upper() == 'G':
                    threeprime_strand = threeprime_strand + 'C'
                elif nb.upper() == 'C':
                    threeprime_strand = threeprime_strand + 'G'
                else:
                    pass

        elif "U" in fiveprime_strand and "T" not in fiveprime_strand:
            for nb in fiveprime_strand:
                if nb.upper() == 'A':
                    threeprime_strand = threeprime_strand + 'U'
                elif nb.upper() == 'U':
                    threeprime_strand = threeprime_strand + 'A'
                elif nb.upper() == 'G':
                    threeprime_strand = threeprime_strand + 'C'
                elif nb.upper() == 'C':
                    threeprime_strand = threeprime_strand + 'G'
                else:
                    pass
        else:
            print('It looks like there is a mix up with RNA and DNA. Check inputs. (CODE 3)')

    elif len(fiveprime_strands) == 2:
        fiveprime_strand = fiveprime_strands[0]
        threeprime_strand = fiveprime_strands[1]

    # print("[" + fiveprime_strand + ">")
    # print("<" + threeprime_strand + "]")
    # print(" ")

    #     print("Strand A: ["+fiveprime_strand+">")
    #     print("Strand B: ["+threeprime_strand[::-1]+">")

    A = Strand(fiveprime_strand, name='Strand A')
    B = Strand(threeprime_strand[::-1], name='Strand B')

    # print_strands(A,B)
    return A, B


def specify_tube(strand_dict, setSpec=2):
    """Inputs: strand_dict: dict of strand names and concentrations (in M)
                 setSpec: int, all complexes of up to setSpec strands
       Output: tube object and complexes object
    """
    strand_names = list(strand_dict.keys())
    strands = {name: Strand(seq, material='dna') for name, seq in strand_dict.items()}

    dict_of_concentrations = {}
    for name, conc in strand_dict.items():
        dict_of_concentrations[name] = make_list(conc)

    tube_dict = {}
    i = 0
    for strand_combination in itertools.combinations_with_replacement(strand_names, setSpec):
        concs = [dict_of_concentrations[strand_name] for strand_name in strand_combination]
        tube_name = ', '.join([f'{name}: {conc} M' for name, conc in zip(strand_combination, concs)])
        tube_strands = {name: conc for name, conc in zip(strand_combination, concs)}
        tube_dict[tube_name] = Tube(strands=tube_strands, complexes=SetSpec(max_size=setSpec), name=tube_name)
        i += 1

    return tube_dict


# def specify_tube(strand_dict, setSpec=2):
#     """Inputs: A & B are <Strand Strand A>, can be list, must be same size,
#                concentrations as molar: e.g. input 1e-9 as 1nm
#                setSpec: int(), all complexes of up to 2 strands
#                name: string of tube name/number.
#        Output: tube object and complexes object
#     """
#     # get the names of the strands in the input dictionary
#     strand_names = list(strand_dict.keys())
#
#     # generate all possible combinations of strand concentrations
#     strand_concentrations = [make_list(concentration) for concentration in strand_dict.values()]
#     strand_combinations = list(itertools.product(*strand_concentrations))
#
#     # create a dictionary of strand names and concentrations
#     strand_concentration_dict = {name: conc for name, conc in zip(strand_names, strand_combinations)}
#
#     # create a list of Tube objects
#     tubes = {}
#     for strand_combination in strand_combinations:
#         tube_strands = {name: conc for name, conc in zip(strand_names, strand_combination)}
#         tube_name = ', '.join([f"{name}: {conc:.2g}M" for name, conc in tube_strands.items()])
#         tubes[tube_name] = Tube(strands=tube_strands, complexes=SetSpec(max_size=setSpec), name=tube_name)
#
#     print_tubes(tubes)
#     return tubes
#
#
# # #### NEED TO IMPROVE THIS FUNCTION TO ANALYSE MORE STRANDS
# # def specify_tube(strand_dict, setSpec=2):  # name=''
# #     # def specify_tube(strand_conc_dict, setSpec=2):
# #     """Inputs: A & B are <Strand Strand A>, can be list, must be same size,
# #                concentrations as molar: e.g. input 1e-9 as 1nm
# #                setSpec: int(), all complexes of up to 2 strands
# #                name: string of tube name/number.
# #        Output: tube object and complexes object
# #     """
# #     # strand_dict= {A: [1e-09, 1e-8], B: [3e-08]}
# #
# #     strand_names = strand_names_catalog()
# #
# #     dict_of_concentrations = {}
# #
# #     strands = list(strand_dict.keys())
# #     A = strands[0]  ## Hard code needs to be updated for generlized cases
# #     B = strands[1]
# #
# #     concentrations = list(strand_dict.values())
# #
# #     ###############################################
# #
# #     for i in range(0, len(concentrations)):
# #         dict_of_concentrations[strand_names[i]] = make_list(concentrations[i])
# #
# #     concentration_A = dict_of_concentrations['A']  ## Hard code needs to be updated for generlized cases
# #     concentration_B = dict_of_concentrations['B']  ## Hard code needs to be updated for generlized cases
# #
# #     ########################### THIS NEEDS WORK
# #     #     function Recurse (y, number)
# #     #        if (number > 1)
# #     #           Recurse ( y, number - 1 )
# #     #        else
# #     #           for x in range (y)
# #     #           whatever()
# #
# #     ##################### This is hard coded for two strands. Need to update to generalize
# #
# #     tubes = {}
# #     i = 0  # counter
# #     for concA in concentration_A:
# #         for concB in concentration_B:
# #             tube_name = 'A: ' + str(concA) + 'M, ' + ' B: ' + str(concB) + 'M'
# #
# #             tubes[str(tube_name)] = Tube(strands={A: concA, B: concB},
# #                                          ######### Hard code needs to be updated for generlized cases
# #
# #                                          complexes=SetSpec(max_size=setSpec),
# #                                          name=tube_name)  # all complexes of up to set spec # of strands
# #             i += 1
# #
# #             # check if list
# #     #     tubes_ = make_list(tubes)
# #     #     print("Number of tubes: "+str(len(tubes_[0])))   ##### debugging
# #
# #     print_tubes(tubes)
# #     return tubes


def analysis_job(tubes, models):
    """Input: tubes: <Tube Tube 1>, can be list of tubes
              models: <nupack.model.Model at XXXX>, can be dictionary
       Output: results from analysis: partition func, ΔG(kcal/mol),concnetration results
    """
    # print("Tube type: "+str(type(tubes)))
    # convert dict or object to list
    if isinstance(models, dict):
        models_list = list(models.values())
    elif not isinstance(models, list):
        models_list = [models]
    else:
        pass

    # convert dict or object to list
    if isinstance(tubes, dict):
        tubes_list = list(tubes.values())
    elif not isinstance(tubes, list):
        tubes_list = [tubes]
    else:
        pass

    results = {}

    i = 0  # counter

    kwargs = {'max_iters': 20000, 'tolerance': 1e-06}

    for model in models_list:
        label = str(list(models.keys())[i])

        ii = 0
        for tube in tubes_list:
            label_full = label + str(list(tubes.keys())[ii])

            # run tube analysis job # Debugging
            # tube = tube[list(tube.keys())[0]] # Debugging
            # model = model[list(model.keys())[0]] # Debugging
            results[label_full] = tube_analysis(tubes=[tube], model=model)  #, **kwargs)

            ii += 1
        i += 1
    return results


def run_analysis(strands_input=['ATGC'], nucleic_acid_type='dna', temperature_C=37,
                 Na_M=1, Mg_M=0,
                 concentrations=[[1e-09], [3e-09]],
                 complex_numb=2):
    """Wrapper function for analysis functions
    """
    if not isinstance(strands_input, list):
        strands_input = [strands_input]

    strand_names = strand_names_catalog()

    # Build models
    models = models_maker(nucleic_acid_type, temperature_C, Na_M, Mg_M)  # rna or dna, temp, salt

    # print(len(strands_input))
    strand_dict = {}
    # print(len(strands_input))
    if len(strands_input) == 1:
        # Build complement strand
        A, B = specify_strands(strands_input[0])
        strand_dict = {A: concentrations[0], B: concentrations[1]}

    #     elif len(strands_input) == 2:
    #         A = Strand(strands_input[0], name='Strand A')
    #         B = Strand(strands_input[1], name='Strand B')

    else:

        # print("Strand A initial concentrations: " + str(concentrations[0]))
        print("len(list(strands_input)): ", len(list(strands_input)))
        for i in range(0, len(list(strands_input))):
            # for ii in range(0, len(concentrations[0])):
            # print(ii)
            # Build a dictionary of all the input strands
            print("strands_input[i]: ", strands_input[i], "\nstrand_names[i]: ", strand_names[i], "\nconcentrations[i]: ", concentrations[i])
            print('Strand: ',Strand(strands_input[i], name='Strand ' + strand_names[i]))
            strand_dict[Strand(strands_input[i], name='Strand_' + strand_names[i])] = concentrations[i]
            # print(strand_dict)

    # Build test tubes
    #     tubes = Tube(strands=strand_dict,
    #                complexes=SetSpec(max_size=2),
    #                     name='tube')

    #     t2 = Tube(strands={A: 1e-10, B: 1e-9}, complexes=SetSpec(max_size=2), name='t2')
    #     # Build test tubes
    #     Tube(strands={A: 1e-6, B: 1e-8, C: 1e-12},
    #          complexes=SetSpec(max_size=complexes))
    ############# NEED to UPDATE HERE FOR MULTIPLE STRAND FUNCTIONALITY
    strandA = str(list(strand_dict.keys())[0])
    strandB = str(list(strand_dict.keys())[1])

    if nucleic_acid_type.lower() == 'dna':
        strandA = convert_rna_to_dna(str(strandA), nucleic_acid_type)
        strandB = convert_rna_to_dna(str(strandB), nucleic_acid_type)

    print("Strand A: [" + strandA + ">")
    print("Strand B: [" + strandB + ">")
    # print("Strand B: ["+strandB[::-1]+">")

    print("Strand A initial concentrations: " + str(concentrations[0]))
    print("Strand B initial concentrations: " + str(concentrations[1]))

    # Build test tubes
    tubes = specify_tube(strand_dict, setSpec=complex_numb)
    # tube_analysis
    # tubes = Tube(strands=strand_dict, complexes=SetSpec(max_size=complexes), name='t1')
    # print((tubes))

    return analysis_job(tubes, models)  # , tubes


#     return tube_analysis(tubes, models)


def make_dataFrame_from_results(Results=pd.DataFrame()):
    """This function creates a simple dataframe from the Nupack query data.
        This function is currently static, and should be rewritten to generalize.
        Inputs: results from make_dataFrame_from_results()
    """
    strand_names = strand_names_catalog()
    column_names = ['Material',
                    'TemperatureC',
                    'Na_molarity (M)',
                    'Mg_molarity (M)',
                    'Complex A Conc (M)',
                    'Complex B Conc (M)',
                    'Complex AA Conc (M)',
                    'Complex AB Conc (M)',
                    'Complex BB Conc (M)',
                    'ΔG A kcal/mol',
                    'ΔG B kcal/mol',
                    'ΔG AA kcal/mol',
                    'ΔG AB kcal/mol',
                    'ΔG BB kcal/mol',
                    'Partition func A kcal/mol',
                    'Partition func B kcal/mol',
                    'Partition func AA kcal/mol',
                    'Partition func AB kcal/mol',
                    'Partition func BB kcal/mol']

    df = pd.DataFrame(columns=column_names)

    for i in range(0, len(list(Results.keys()))):

        t_result = Results[list(Results.keys())[i]].tubes
        c_result = Results[list(Results.keys())[i]].complexes

        # dir(c_result[list(c_result.keys())[0]].model) # Debugging

        ############ These indices need to be updated!

        free_energy_B = '{:.2g}'.format(c_result[list(c_result.keys())[0]].free_energy)
        free_energy_A = '{:.2g}'.format(c_result[list(c_result.keys())[1]].free_energy)
        free_energy_AA = '{:.2g}'.format(c_result[list(c_result.keys())[2]].free_energy)
        free_energy_AB = '{:.2g}'.format(c_result[list(c_result.keys())[3]].free_energy)
        free_energy_BB = '{:.2g}'.format(c_result[list(c_result.keys())[4]].free_energy)

        pfunc_B = '{:.2g}'.format(c_result[list(c_result.keys())[0]].pfunc)
        pfunc_A = '{:.2g}'.format(c_result[list(c_result.keys())[1]].pfunc)
        pfunc_AA = '{:.2g}'.format(c_result[list(c_result.keys())[2]].pfunc)
        pfunc_AB = '{:.2g}'.format(c_result[list(c_result.keys())[3]].pfunc)
        pfunc_BB = '{:.2g}'.format(c_result[list(c_result.keys())[4]].pfunc)

        # Other options
        #     mfe_B  = c_result[list(c_result.keys())[0]].mfe
        #     mfe_A  = c_result[list(c_result.keys())[1]].mfe
        #     mfe_AA = c_result[list(c_result.keys())[2]].mfe
        #     mfe_AB = c_result[list(c_result.keys())[3]].mfe
        #     mfe_BB = c_result[list(c_result.keys())[4]].mfe

        material = c_result[list(c_result.keys())[0]].model.material

        temperatureK = c_result[list(c_result.keys())[0]].model.temperature  # In kelvin
        temperatureC = convertK_to_C(temperatureK)

        Na_molarity = c_result[list(c_result.keys())[0]].model.conditions.na_molarity
        Mg_molarity = c_result[list(c_result.keys())[0]].model.conditions.mg_molarity

        concentrations = t_result[list(t_result.keys())[0]].complex_concentrations  # iteratable

        ############################################## Consider making this a function
        key_value_pair = list(concentrations.items())

        for oldkey, conc in key_value_pair:
            # print(str(oldkey.strands)) # verification
            string = str((oldkey.strands)).replace('Strand', '')

            for element in string:
                # print(string)
                if element not in strand_names:
                    # print(element)

                    string = string.replace(str(element), '')

                string = ''.join(sorted(string))
            # print(string)   # verification

            concentrations[string] = concentrations.pop(oldkey)
        ##############################################

        complex_A = '{:.3g}'.format(concentrations['A'])
        complex_B = '{:.3g}'.format(concentrations['B'])
        complex_AA = '{:.3g}'.format(concentrations['AA'])
        complex_BB = '{:.3g}'.format(concentrations['BB'])
        complex_AB = '{:.3g}'.format(concentrations['AB'])

        df.loc[i] = [material, temperatureC, Na_molarity, Mg_molarity,
                     complex_A, complex_B, complex_AA, complex_AB, complex_BB,
                     free_energy_A, free_energy_B, free_energy_AA, free_energy_AB, free_energy_BB,
                     pfunc_A, pfunc_B, pfunc_AA, pfunc_AB, pfunc_BB]
    return df, t_result, c_result


# Building a dictionary of each experimental database
def make_sorted_df(df_Nupack, concentrations):
    """This function is used to get and sort the same experimental conditions in the df_Nupack
        dataframe, and make a new simple dataframe to be used with plotting. The plotting is
        specifically set up for "TempC","percent_bound_A","percent_bound_B","AB_Conc_M", "A_Conc_M", "B_Conc_M"
        but this can be changed at a later date.
        Input: df_Nupack- dataframe of nupack results
        Output: Dataframe of sorted chosen data
    """

    concentration_A = make_list(concentrations[0])
    concentration_B = make_list(concentrations[1])

    # making the column names easier to handle.
    column_names = ['Material',
                    'TempC',
                    'Na_M',
                    'Mg_M',
                    'A_Conc_M',
                    'B_Conc_M',
                    'AA_Conc_M',
                    'AB_Conc_M',
                    'BB_Conc_M',
                    'ΔG_A',
                    'ΔG_B',
                    'ΔG_AA',
                    'ΔG_AB',
                    'ΔG_BB',
                    'Pfunc_A',
                    'Pfunc_B',
                    'Pfunc_AA',
                    'Pfunc_AB',
                    'Pfunc_BB']

    df_Nupack.columns = column_names  # setting column names

    conditions = len(df_Nupack.TempC.unique())  # +len(df_Nupack.Na_M.unique()+len(df_Nupack.Mg_M.unique()
    dfsize = int(df_Nupack.shape[0])
    unique_rows = int(dfsize / conditions)

    sorted_df_Nupack = {}
    for i in range(0, unique_rows):
        sorted_df_Nupack[i] = df_Nupack.iloc[i::unique_rows, :]  # Splicing out rows for each experiment

    temp_concBound = {}

    # Building a new sorted & simple dictionary for plotting
    # Also adding percent bound of species
    for i in range(0, len(sorted_df_Nupack)):
        sorted_df_Nupack[i].reset_index(inplace=True)

        # Making new dataframe to store data I want to graph
        df_temp_concBound = pd.DataFrame(columns=["TempC",
                                                  "Na_M",
                                                  "Mg_M",
                                                  "percent_bound_A",
                                                  "percent_bound_B",
                                                  "AB_Conc_M",
                                                  "A_Conc_M",
                                                  "B_Conc_M"])

        for iii in range(0, len(concentration_A)):  # this for loop cycles through the concentrations of A
            ## Iterating over each set of ordered experiments
            for ii in range(0, len(sorted_df_Nupack[i])):
                TempC = sorted_df_Nupack[i]["TempC"][ii]
                Na_M = sorted_df_Nupack[i]["Na_M"][ii]
                Mg_M = sorted_df_Nupack[i]["Mg_M"][ii]
                AB_Conc_M = float(sorted_df_Nupack[i]["AB_Conc_M"][ii])

                A_Conc_M = float(sorted_df_Nupack[i]["A_Conc_M"][ii])
                B_Conc_M = float(sorted_df_Nupack[i]["B_Conc_M"][ii])

                # Note: This is an important operation for the outcome of this code  ###################################
                # These are the percent bound [AB/A] and [AB/B]
                percent_bound_A = ((AB_Conc_M) / float(concentration_A[iii])) * 100  # [AB/A]
                percent_bound_B = ((AB_Conc_M) / float(concentration_B[0])) * 100  # [AB/B]

                ## Writing new row in the concentration bound dataframe
                df_temp_concBound.loc[ii] = [TempC, Na_M, Mg_M, percent_bound_A, percent_bound_B, AB_Conc_M, A_Conc_M,
                                             B_Conc_M]

            # Dictonary for df_temp_concBound values
            temp_concBound[i] = df_temp_concBound

    return temp_concBound

from itertools import cycle

def plot_experiment(temp_concBound, concentration_A, sampleVolume_uL, targetAmount, species="B"):
    """This function plots the analyzed results found from Nupack
        Input: dictionary of temperatures and % bound species, sample volume in μL,
               target molc amount
        Output: Plot of results
    """

    # handling species input
    if species.upper() == "B":
        dependent_var = 'Percent bound [AB]/[B$]_{0}$'
        bound_species = "percent_bound_B"

    elif species.upper() == "A":
        dependent_var = "Percent bound [AB]/[A]"
        bound_species = "percent_bound_A"

    # Custom title

    if species == "B":
        nout = "o"
    else:
        nout = ""
    title_seg1 = ("Percent bound ([AB]/[" + str(species) + "$]_{0}$" + "). " + str(
        "{:.2e}".format(targetAmount)) + " copies of target in " +
                  str(sampleVolume_uL) + " μL")

    fig, axes = plt.subplots(figsize=(12, 8), nrows=1, ncols=1)

    colors = cycle(['tab:blue', 'tab:orange', 'tab:red', 'tab:green', 'tab:purple',
                    'tab:pink', 'tab:cyan', 'tab:olive', 'tab:brown', 'k'])

    for i in range(0, len(temp_concBound)):

        Na_M = temp_concBound[i]["Na_M"]
        Mg_M = temp_concBound[i]["Mg_M"]

        flag_1 = True  # A generic boolean for formatting new lines in the title bar
        ions = ''

        if len(np.unique(Na_M.squeeze())) == 1:
            title_seg2 = '   $Na^+ (M): $' + str(Na_M[i])
            title = title_seg1 + '.' + '\n' + title_seg2
        else:
            flag_1 = False
            ions = ions + title_seg2  # part of legend label for Na salts
            title_seg2 = '   $Na^+ (M): $' + ', '.join(map(str, list(Na_M.squeeze().unique())))
            title = title_seg1 + '\n' + title_seg2

        if len(np.unique(Mg_M.squeeze())) == 1 and flag_1 == True:
            title_seg3 = ' $Mg^{++} (M): $' + str(Mg_M[i])
            title = title + '\n' + title_seg3
        else:
            ions = ions + title_seg3  # part of legend label for Na salts
            title = title + '\n' + '$Mg^{++} (M): $' + ', '.join(map(str, list(Mg_M.squeeze().unique())))

        legend_label = 'Probe Conc: ' + str(concentration_A[i]) + ' M. ' + ions

        X = temp_concBound[i]["TempC"]
        Y = (temp_concBound[i][bound_species])
        axes.plot(X, Y, color=next(colors), markersize=10, marker='o', fillstyle='none', linewidth=2, label=legend_label)

    axes.set_xlabel('Temperature [C]', fontsize=18)
    axes.set_ylabel(dependent_var, fontsize=18)
    axes.set_title(title, fontsize=20, pad=20)
    axes.grid(True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='large')

    fig.tight_layout()
    plt.show()
    # plt.savefig('figures/test.png')
    return fig

# def plot_experiment(temp_concBound, concentration_A, sampleVolume_uL, targetAmount, species="B"):
#     """This function plots the analyzed results found from Nupack
#         Input: dictionary of temperatures and % bound species, sample volume in μL,
#                target molc amount
#         Output: Plot of results
#     """
#
#     # handling species input
#     if species.upper() == "B":
#         dependent_var = 'Percent bound [AB]/[B$]_{0}$'
#         bound_species = "percent_bound_B"
#
#     elif species.upper() == "A":
#         dependent_var = "Percent bound [AB]/[A]"
#         bound_species = "percent_bound_A"
#
#     # Custom title
#     if species == "B":
#         nout = "o"
#     else:
#         nout = ""
#     title_seg1 = ("Percent bound ([AB]/[" + str(species) + "$]_{0}$" + "). " + str(
#         "{:.2e}".format(targetAmount)) + " copies of target in " +
#                   str(sampleVolume_uL) + " μL")
#
#     fig, axes = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
#
#     colors = cycle(['tab:blue', 'tab:orange', 'tab:red', 'tab:green', 'tab:purple',
#                     'tab:pink', 'tab:cyan', 'tab:olive', 'tab:brown', 'k'])
#
#     ####
#
#     ##
#     for i in range(0, len(temp_concBound)):
#
#         Na_M = temp_concBound[i]["Na_M"]
#         Mg_M = temp_concBound[i]["Mg_M"]
#
#         flag_1 = True  # A generic boolian for formatting new lines in the title bar
#         ions = ''
#
#         if len(np.unique(Na_M.squeeze())) == 1:
#             title_seg2 = '   $Na^+ (M): $' + str(Na_M[i])
#             #             ions = ions + title_seg2 # part of legend label for Na salts
#
#             title = title_seg1 + '.' + '\n' + title_seg2
#         else:
#             flag_1 = False
#             ions = ions + title_seg2  # part of legend label for Na salts
#             title_seg2 = '   $Na^+ (M): $' + ', '.join(map(str, list(Na_M.squeeze().unique())))
#             title = title_seg1 + '\n' + title_seg2
#
#         if len(np.unique(Mg_M.squeeze())) == 1 and flag_1 == True:
#             title_seg3 = ' $Mg^{++} (M): $' + str(Mg_M[i])
#
#             #             ions = ions + title_seg3 # part of legend label for Na salts
#             title = title + '\n' + title_seg3
#         else:
#             ions = ions + title_seg3  # part of legend label for Na salts
#             title = title + '\n' + '$Mg^{++} (M): $' + ', '.join(map(str, list(Mg_M.squeeze().unique())))
#
#         legend_label = 'Probe Conc: ' + str(concentration_A[i]) + ' M. ' + ions
#
#         X = temp_concBound[i]["TempC"]
#         Y = (temp_concBound[i][bound_species])
#         axes.plot(X, Y, color=next(colors), markersize=8, marker='o', fillstyle='none', linewidth=2, label=legend_label)
#     axes.set_xlabel('Temperature [C]', fontsize=16)
#     axes.set_ylabel(dependent_var, fontsize=16)
#     axes.set_title(title, fontsize=16)
#     plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='large')
#
#     fig.tight_layout()
#     plt.show()
#     plt.savefig('figures/test.png')
#     return fig


def run_all(strands='ATGC', nucleic_acid_type='dna', temperature_C=37, sampleVolume_uL=25,
            concentrations=[[1e-09], [3e-09]], Na_M=[0.15], Mg_M=[0.01], complexes=2, species="A",
            strand_dict=None, savefig=False):
    """Runs a NUPACK analysis and returns the results.

    Args:
        strands (list of str): The nucleic acid sequences of the strands.
        nucleic_acid_type (str): The type of nucleic acid ('RNA' or 'DNA').
        temperature_C (list of float): The temperatures at which to run the analysis.
        sample_volume_uL (float): The sample volume in microliters.
        concentrations (list of lists of float): The concentrations of each strand in molar.
        Na_M (list of float): The sodium ion concentrations in molar.
        Mg_M (list of float): The magnesium ion concentrations in molar.
        complexes (int): The maximum size of complexes to consider.
        species (str): The species for which to calculate concentrations.
        strand_dict (dict): A dictionary of strand names and concentrations.
        savefig (bool): Whether to save a figure of the analysis.

    Returns:
        A tuple containing the results of the analysis.
    """
    concentration_A = make_list(concentrations[0])
    # print(concentration_A)
    concentration_B = make_list(concentrations[1])

    #     print('concentration_B' + str(concentration_B))
    #     print(temperature_C)
    results = run_analysis(strands, nucleic_acid_type, temperature_C,
                           Na_M, Mg_M,
                           concentrations, complexes)

    print(report_parameters(results)[0][5:-26])  # Displays results and information

    df_Nupack, t_result, c_result = make_dataFrame_from_results(results)  ###################################

    temp_concBound = make_sorted_df(df_Nupack, concentrations)

    targetAmount = convert_concentration_to_copy_number(concentration_B[0], sampleVolume_uL)

    ##### Comment out plotting for debugging
    fig = plot_experiment(temp_concBound, concentration_A, sampleVolume_uL, targetAmount, species)

    if savefig:  #.lower():
        save_fig(fig)

    return results, df_Nupack, temp_concBound, t_result, c_result

    #
    # if strands is None and strand_dict is None:
    #     raise ValueError("Either 'strands' or 'strand_dict' must be provided")
    #
    # if strand_dict is not None:
    #     strands = list(strand_dict.keys())
    #
    # concentration_lists = []
    # for i in range(len(strands)):
    #     if concentrations is not None and i < len(concentrations):
    #         concentration_lists.append(make_list(concentrations[i]))
    #     else:
    #         concentration_lists.append([1e-9])
    #
    # tubes = specify_tube(strand_dict, setSpec=complexes)
    #
    # Results = run_analysis(strands, nucleic_acid_type, temperature_C,
    #                        Na_M, Mg_M,
    #                        concentration_lists, complexes, tubes)
    #
    # print(report_parameters(Results)[0][5:-26])  # Displays results and information
    #
    # df_Nupack, t_result, c_result = make_dataFrame_from_results(Results)
    # temp_concBound = make_sorted_df(df_Nupack, concentration_lists)
    #
    # targetAmount = convert_concentration_to_copy_number(concentration_lists[-1][0], sampleVolume_uL)
    # fig = plot_experiment(temp_concBound, concentration_lists[:-1], sampleVolume_uL, targetAmount, species)
    #
    #
    # if savefig:  #.lower():
    #     save_fig(fig)
    #
    # return Results, df_Nupack, temp_concBound, t_result, c_result

