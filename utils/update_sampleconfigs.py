import os
import ruamel.yaml

# --- CONFIG ---
yaml_file = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/Samples/OA_samples_subsamples/AlltheMC/pTpzEnu/0m_all.yaml"
data_dir = "/project/rpp-nilic/abipeake/PRISM_CAFs/FHC/FHC_skimmed/"
prefix = "CAFv7_0m_"
suffix = ".root"
splinefile = ""
nutype = 14
oscnutype = 14
signal = False
# ---------------

# Use ruamel.yaml (preserves formatting, indentation, quotes)
yaml = ruamel.yaml.YAML()
yaml.preserve_quotes = True
yaml.indent(mapping=2, sequence=4, offset=2)

# Load YAML
with open(yaml_file, "r") as f:
    data = yaml.load(f)

# Make sure SubSamples exists
if "SubSamples" not in data or not isinstance(data["SubSamples"], list):
    data["SubSamples"] = []

# Find all matching ROOT files
files = sorted([
    f for f in os.listdir(data_dir)
    if f.startswith(prefix) and f.endswith(suffix)
])
if not files:
    raise RuntimeError(f"❌ No files found in {data_dir} matching {prefix}*{suffix}")

# Determine existing entries
existing = {s["mtuplefile"] for s in data["SubSamples"] if "mtuplefile" in s}
next_index = max(
    (s.get("samplevecno", -1) for s in data["SubSamples"] if isinstance(s.get("samplevecno"), int)),
    default=-1
) + 1

new_entries = []
skipped_entries = []

# Add new minimal entries
for fname in files:
    base = os.path.splitext(fname)[0]
    if base in existing:
        skipped_entries.append(base)
        continue
    entry = {
        "Name": base,
        "LatexName": base,
        "mtuplefile": base,
        "splinefile": splinefile,
        "samplevecno": next_index,
        "nutype": nutype,
        "oscnutype": oscnutype,
        "signal": signal
    }
    data["SubSamples"].append(entry)
    next_index += 1
    new_entries.append(base)

# Optional: sort SubSamples by samplevecno
data["SubSamples"] = sorted(data["SubSamples"], key=lambda x: x["samplevecno"])

# Update NSubSamples
data["NSubSamples"] = len(data["SubSamples"])

# Write back the updated YAML
with open(yaml_file, "w") as f:
    yaml.dump(data, f)

# Print summary
if new_entries:
    print(f"✅ Added {len(new_entries)} new SubSamples:")
    for e in new_entries:
        print(f"  - {e}")
else:
    print("ℹ️ No new CAFv7_0m_*.root files to add.")

if skipped_entries:
    print(f"⚠️ Skipped {len(skipped_entries)} files already present in YAML:")
    for e in skipped_entries:
        print(f"  - {e}")

print(f"✅ YAML file '{yaml_file}' now has {data['NSubSamples']} total SubSamples.")
