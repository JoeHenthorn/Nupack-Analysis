from npktools_dict import run_all, convert_copy_number_to_concentration, is_nucleic_acid

default_inputs = {
    'num_strands': 2,
    'strands': None,
    'nucleic_acid_type': 'DNA',
    'temperature_C': [37],
    'complexes': 2,
    'concentration_dict': None,
    'Na_M': [0.150],
    'Mg_M': [0.010],
    'species': None,
    'sampleVolume_uL': 25
}

user_inputs = {}


def validate_user_inputs(user_inputs):
    # Validate that the number of strands is an integer
    assert isinstance(user_inputs['num_strands'], int)

    # Validate that each strand name is a string
    for name in user_inputs:
        if name != 'num_strands':
            assert isinstance(name, str)

    # Validate that each strand sequence is a DNA or RNA sequence
    nucleic_acid_type = user_inputs['nucleic_acid_type']
    assert nucleic_acid_type in ['DNA', 'RNA']
    for sequence in user_inputs.values():
        if isinstance(sequence, str):
            assert is_nucleic_acid(sequence, nucleic_acid_type)

    # Additional validation...

    # If nucleic acid type is RNA, validate U's in sequence
    if nucleic_acid_type == 'RNA':
        for sequence in user_inputs.values():
            assert 'U' in sequence

    # Validate that concentration_dict keys match the strand names
    if 'concentration_dict' in user_inputs:
        assert set(user_inputs['concentration_dict'].keys()) == set(user_inputs.keys()) - {'num_strands'}

    # Validate the temperature range
    assert isinstance(user_inputs['temp_startC'], float)
    assert isinstance(user_inputs['temp_stopC'], float)
    assert isinstance(user_inputs['temp_intervalC'], float)
    assert user_inputs['temp_startC'] < user_inputs['temp_stopC']

    # Validate the species is one of the strand names
    assert user_inputs['species'] in user_inputs.keys()

    # Additional validation...

    # If all valid, return the validated user_inputs
    return user_inputs

def run_analysis(default_inputs, user_inputs):
    # Update default_inputs with validated user_inputs
    default_inputs.update(user_inputs)

    # Run the analysis
    results, df_Nupack, temp_concBound, t_result, c_result = run_all(**default_inputs, savefig=True)

    return results, df_Nupack, temp_concBound, t_result, c_result


if __name__ == "__main__":
    # Print default number of strands and prompt user
    print(f'Default number of strands: {default_inputs["num_strands"]}')
    num_strands = input('Enter new value or press Enter to accept: ')

    # Check if user entered new value or accepted default
    if num_strands:
        num_strands = int(num_strands)
    else:
        num_strands = default_inputs['num_strands']

    # Add number of strands to user_inputs
    user_inputs['num_strands'] = num_strands

    # Dynamically get strand info from user, accepting default values
    for i in range(num_strands):
        # Print default name and prompt user
        print(f'Default name for strand {i + 1}: {chr(i + 65)}')
        name = input('Enter new name or press Enter to accept: ')

        # Check if user entered new name or accepted default
        if name:
            user_inputs[name] = input(f'Enter the sequence of strand {name}: ')
        else:
            user_inputs[chr(i + 65)] = input(f'Enter the sequence of strand {chr(i + 65)}: ')

    # Print default nucleic acid type
    print(f'Default nucleic acid type: {default_inputs["nucleic_acid_type"]}')

    # Prompt user to enter new value or press Enter
    nucleic_acid_type = input('Enter new value or press Enter to accept: ')

    # Check if user entered new value or pressed Enter
    if nucleic_acid_type:
        user_inputs['nucleic_acid_type'] = nucleic_acid_type
    else:
        user_inputs['nucleic_acid_type'] = default_inputs['nucleic_acid_type']

    # Validate and store user inputs
    user_inputs = validate_user_inputs(user_inputs)

    # Run analysis with default and user inputs
    results, df_Nupack, temp_concBound, t_result, c_result = run_analysis(default_inputs, user_inputs)