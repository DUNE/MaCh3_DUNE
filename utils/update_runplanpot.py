import os
import re

def update_pot_values_inplace(directory, total_pot, custom_pot_file, custom_pot_percentage):
    yaml_files = [f for f in os.listdir(directory) if f.endswith(".yaml")]
    
    if custom_pot_file not in yaml_files:
        print(f"Error: {custom_pot_file} not found in {directory}")
        return
    
    other_files = [f for f in yaml_files if f != custom_pot_file]
    num_other_files = len(other_files)
    
    if num_other_files == 0:
        print("Error: No other YAML files found to distribute POT.")
        return
    
    custom_pot_value = total_pot * (custom_pot_percentage / 100)
    remaining_pot = total_pot - custom_pot_value
    pot_per_file = remaining_pot / num_other_files
    
    def modify_pot_in_file(file_path, new_pot):
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Match POT line like: POT: <value> (with optional spaces)
        new_content, count = re.subn(r'^(POT\s*:\s*).+$', r'\g<1>{:.6e}'.format(new_pot),
                                     content, flags=re.MULTILINE)
        if count == 0:
            print(f"Warning: No POT key found in {file_path}, skipping.")
            return
        
        with open(file_path, 'w') as f:
            f.write(new_content)
    
    # Update the custom POT file
    modify_pot_in_file(os.path.join(directory, custom_pot_file), custom_pot_value)
    
    # Update other files
    for file in other_files:
        modify_pot_in_file(os.path.join(directory, file), pot_per_file)
    
    print("POT values updated successfully.")

# Example usage
update_pot_values_inplace(
    "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/Samples/OA_samples_subsamples/Off_axis_studies/Off_axis_q0q3_ptpzEnu/onaxis75perecent_offaxis25perecnt/",
    6.6e21,
    "0m_all.yaml",
    75
)
