#%% Set up parameters and run analysis
'''
Type in your parameters here and run the cell to get your results.

Strand A: (Note: right now the reverse complement stand is automatically generated with only one input)
Strand B:
Nucleic acid type: RNA or DNA
Complex number
Temperature start (˚C)
Temperature stop (˚C)
Temperature interval (˚C)
Concentration A (M):
Sample volume (μL)
Target/template copy number for Strand B
Sodium concentration (M). Does not yet have list compatibility.
Magnesium concentration (M). Does not yet have list compatibility.
Species: Use "A" for percent bound A
Species: Use "B" for percent bound target/template
Saving figures: If you want to save the figure, include "save", otherwise it doesn't matter what you put.
'''
import math

from npktools import run_all, convert_copy_number_to_concentration

# Insert your sequence here! NOTE: You do not need to type in a complementary strand_2 ...
# ... if the second strand is fully complementary to the first.
# Strand A  5'- 3'
# Strand B  5'- 3'      Note: Nupack web app is 5'- 3' for both inputs
strand_A = 'CCTAGCTCTGACCACTTCACACCTACGACCACAGATGGTGAGGATGGCTGACCTGACGTACGGCTCTCA'  # primer or probe  ####  v23 MS2 G probe 11/4/22
strand_B = 'GTTCACCTTCAAGAGTTTCTTCCTATGAGAGCCGTACGTCAGGTC'  # target/template  ### v23 MS2 Target template (complement to probe sequence) 11/4/22

strands = [strand_A, strand_B]
#strands = [strand_A]
nucleic_acid_type ='DNA'    ####  Input here

complexes = 2               ####  Input here

temp_startC = 50            ####  Input here
temp_stopC = 100            ####  Input here
temp_intervalC = 3          ####  Input here

temperatureC = list(range(temp_startC,temp_stopC+1,temp_intervalC))

concentration_A_M = [1e-5, 1e-6, 1e-7, 1e-8] #, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]  #, 1e-15, 1e-16, 1e-17] # Capture
# probe

# Strand B
sampleVolume_uL = 25        ####  Input here
targetCopyNumber = 10000    ####  Input here # 2.5e6

# Don't change! This code assumes you are inputting the target copy number.
# If you know the exact concentration(s), make a list and add them into the concentrations list variable below.
#concentration_B_M = [1e-6]#, 1e-10, 1e-11, 1e-12] # Target
concentration_B_M = [convert_copy_number_to_concentration(sampleVolume_uL, targetCopyNumber)] # Default: Target 1000 DNA/RNA's in 25μL

concentrations = [concentration_A_M, concentration_B_M]

# Na:  MIN: 0.05  MAX value: 1.1
Na_M = [0.150]               ####  Input here

# Mg:  MIN: 0.0   MAX value: 0.2
Mg_M = [0.00]                ####  Input here

# [AB]/[species]
species = "B"                ####  Input here

# run_all() # Works with no inputs for default values.

if __name__ == "__main__":
    # Run the analysis
    Results, df_Nupack, temp_concBound, t_result, c_result = run_all(strands,
                                                                     nucleic_acid_type,
                                                                     temperatureC,
                                                                     sampleVolume_uL,
                                                                     concentrations,
                                                                     Na_M, Mg_M,
                                                                     complexes,
                                                                     species,savefig=True)