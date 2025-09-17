#!/usr/bin/env python3
import yaml

# Input/output files
input_file = "/vols/dune/jmm224/MaCh3/MaCh3_DUNE/Configs/EventRates_Beam_NDGAr_clean.yaml"
output_file = "/vols/dune/jmm224/MaCh3/MaCh3_DUNE/Configs/NewEventRates_Beam_NDGAr.yaml"

# Recursive function to convert simple lists to flow style
def convert_lists(obj):
    if isinstance(obj, list):
        # If all elements are simple scalars, mark for flow style
        if all(isinstance(i, (int, float, str, bool, type(None))) for i in obj):
            obj.fa_flow_style = True  # note: PyYAML doesn't support this directly, handled below
            return obj
        else:
            return [convert_lists(i) for i in obj]
    elif isinstance(obj, dict):
        return {k: convert_lists(v) for k, v in obj.items()}
    else:
        return obj

# Load YAML
with open(input_file, "r") as f:
    data = yaml.safe_load(f)

# Dump YAML using flow style for simple lists
def represent_list(dumper, data):
    if all(isinstance(i, (int, float, str, bool, type(None))) for i in data):
        return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=False)

yaml.add_representer(list, represent_list)

# Dump cleaned YAML
with open(output_file, "w") as f:
    yaml.safe_dump(data, f, default_flow_style=False, sort_keys=False)

print(f"Clean YAML written to {output_file}")
