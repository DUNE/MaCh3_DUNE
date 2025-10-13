#!/usr/bin/env python3
import os
import re
import shutil
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq

# --- CONFIG ---
yaml_dir = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/Samples/OA_samples_subsamples/AlltheMC/pTpzEnu/"
data_dir = "/project/rpp-nilic/abipeake/PRISM_CAFs/FHC/"   # where the CAFv7_*.root files live
output_dir = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/updated_yamls/"

suffix = ".root"
splinefile = ""
nutype = 14
oscnutype = 14
signal = False

recursive = False   # set True to walk yaml_dir recursively and preserve subdirs in output_dir
# ---------------

yaml = YAML()
yaml.preserve_quotes = True
yaml.indent(mapping=2, sequence=4, offset=2)

os.makedirs(output_dir, exist_ok=True)

# find YAML files
yaml_files = []
for root, _, files in os.walk(yaml_dir):
    for fn in files:
        if fn.endswith(".yaml"):
            yaml_files.append(os.path.join(root, fn))
    if not recursive:
        break

if not yaml_files:
    raise RuntimeError(f"No .yaml files found in {yaml_dir!r}")

# build a SAFE textual Binning block exactly as you want it to appear
DESIRED_BINNING_TEXT = """Binning:
  Axes:
    - VarStr: pT
      Uniform: [200,0,10]
      Title: "p_{T} [GeV]"
    - VarStr: pz
      Uniform: [25,-0.5,2.5]
      Title: "p_{z} [GeV]"
    - VarStr: RecoNeutrinoEnergy
      VarBins: [0., 0.5, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 5., 6., 10.]
      Title: "E^{Reco.}_{#nu} [GeV]"
"""

def replace_binning_in_text_file(path, desired_block):
    """Replace existing top-level 'Binning:' block (if present) with desired_block.
       If not present, append the block at the end of file.
    """
    with open(path, "r") as fh:
        txt = fh.read()

    # regex: from a line starting with 'Binning:' up to (but not including) the next top-level key '^\S.*?:' or end-of-file
    pat = re.compile(r'(?ms)^Binning:.*?^(?=\S.*?:|\Z)')
    if pat.search(txt):
        new_txt = pat.sub(desired_block.rstrip() + "\n", txt, count=1)
        action = "replaced"
    else:
        # append at EOF (keep a blank line before)
        new_txt = txt.rstrip() + "\n\n" + desired_block.rstrip() + "\n"
        action = "appended"

    with open(path, "w") as fh:
        fh.write(new_txt)
    return action

# process each YAML
total_new = 0
total_skipped = 0

for yf in yaml_files:
    print(f"\n--- Processing {yf} ---")

    # load yaml
    with open(yf, "r") as fh:
        data = yaml.load(fh)

    # make a backup next to original
    backup = yf + ".bak"
    shutil.copy2(yf, backup)
    print(f"Backup created: {backup}")

    # derive token to limit CAFv7 matches (tries filename, parent dir, then regex search)
    base = os.path.splitext(os.path.basename(yf))[0]
    token = None
    token_cand = base.split("_")[0]
    if re.match(r'^\d+(\.\d+)?m$', token_cand):
        token = token_cand
    else:
        parent = os.path.basename(os.path.dirname(yf))
        if re.match(r'^\d+(\.\d+)?m$', parent):
            token = parent
        else:
            m = re.search(r'\d+(\.\d+)?m', base)
            if m:
                token = m.group(0)

    if token:
        pattern = re.compile(rf"^CAFv7_{re.escape(token)}_.*{re.escape(suffix)}$")
        print(f"Using token '{token}' -> will only add files matching CAFv7_{token}_*.root")
    else:
        pattern = re.compile(r"^CAFv7_\d.*\.root$")
        print("WARNING: could not derive a '<number>m' token from YAML name; falling back to broad CAFv7_<digit> matcher")

    # list matching files for THIS yaml only
    try:
        file_list = sorted([f for f in os.listdir(data_dir) if pattern.match(f)])
    except FileNotFoundError:
        raise RuntimeError(f"data_dir {data_dir!r} not found or not accessible on this node")

    print(f"Found {len(file_list)} matching CAF files in {data_dir}")

    # ensure SubSamples exists
    if "SubSamples" not in data or not isinstance(data["SubSamples"], list):
        data["SubSamples"] = []

    existing = {s["mtuplefile"] for s in data["SubSamples"] if "mtuplefile" in s}
    next_index = max(
        (s.get("samplevecno", -1) for s in data["SubSamples"] if isinstance(s.get("samplevecno"), int)),
        default=-1
    ) + 1

    new_entries = []
    skipped_entries = []

    for fname in file_list:
        base_noext = os.path.splitext(fname)[0]
        if base_noext in existing:
            skipped_entries.append(base_noext)
            continue
        entry = {
            "Name": base_noext,
            "LatexName": base_noext,
            "mtuplefile": base_noext,
            "splinefile": splinefile,
            "samplevecno": next_index,
            "nutype": nutype,
            "oscnutype": oscnutype,
            "signal": signal
        }
        data["SubSamples"].append(entry)
        new_entries.append(base_noext)
        next_index += 1

    # sort & count
    data["SubSamples"] = sorted(data["SubSamples"], key=lambda x: x["samplevecno"])
    data["NSubSamples"] = len(data["SubSamples"])

    # set a Binning structure so ruamel includes the key (actual formatting will be replaced afterwards)
    data["Binning"] = CommentedMap()
    data["Binning"]["Axes"] = CommentedSeq()  # placeholder; exact textual block will be injected

    # determine output_path (preserve subdir structure if recursive)
    rel = os.path.relpath(yf, yaml_dir)
    out_path = os.path.join(output_dir, rel)
    out_dir = os.path.dirname(out_path)
    os.makedirs(out_dir, exist_ok=True)

    # dump with ruamel (so other formatting/preserving details are handled)
    with open(out_path, "w") as fh:
        yaml.dump(data, fh)

    # Now forcibly replace the Binning text with the exact desired inline-list block:
    action = replace_binning_in_text_file(out_path, DESIRED_BINNING_TEXT)
    print(f"Binning block {action} in {out_path}")

    print(f"Added {len(new_entries)} new entries, skipped {len(skipped_entries)} already-present")
    print(f"Saved updated YAML -> {out_path}")

    total_new += len(new_entries)
    total_skipped += len(skipped_entries)

# final summary
print("\n=== Batch Summary ===")
print(f"Processed {len(yaml_files)} YAML files")
print(f"Total new SubSamples added: {total_new}")
print(f"Total skipped (already present): {total_skipped}")
print(f"Updated YAMLs written to: {output_dir}")
print("======================\n")
