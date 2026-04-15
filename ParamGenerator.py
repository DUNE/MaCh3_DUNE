import yaml

# -----------------------------
# User configuration
# -----------------------------
NumBins = 40
Emin = 0.0
Emax = 10.0

# Map PDG → human-readable neutrino names
flavour_names = {
    12:  "nue",
    14:  "numu",
    16:  "nutau",
    -12: "nuebar",
    -14: "numubar",
    -16: "nutaubar",
}

# List of (Unosc, Osc) flavour pairs
channels = [
    (14, 14),
    (12, 12),
    (-14, -14),
    (-12, -12),
    (14, 12),
    (-14, -12),
    (12, 14),
    (-12, -14),
    (14, 16),
    (12, 16),
    (-14, -16),
    (-12, -16),
]

# Samples — each gets a full 480-param set
# samples = ["FD_FHC_numu", "FD_FHC_nue", "FD_RHC_numu", "FD_RHC_nue"]

# -----------------------------
# Build systematics
# -----------------------------
systematics = []
bin_width = (Emax - Emin) / NumBins

# for sample in samples:

for (nu_unosc, nu_osc) in channels:

    unosc_name = flavour_names[nu_unosc]
    osc_name   = flavour_names[nu_osc]

    for b in range(NumBins):
        E_low = Emin + b * bin_width
        E_high = E_low + bin_width

        # Format energy nicely for names
        E_low_str = f"{E_low:.2f}"
        E_high_str = f"{E_high:.2f}"

        # Construct descriptive parameter name:
        #   sample_unosc_osc_0.00_0.25
        # pname = f"{sample}_{unosc_name}_{osc_name}_{E_low_str}_{E_high_str}"
        pname = f"{unosc_name}_{osc_name}_{E_low_str}_{E_high_str}"

        block = {
            "Systematic": {
                "SampleNames": "FD_*",
                "Error": 1.0,
                "FlatPrior": True,
                "FixParam": False,
                "KinematicCuts": [
                    {"TrueNeutrinoEnergy": [round(E_low, 6), round(E_high, 6)]}
                ],
                "Names": {
                    "FancyName": pname,
                    "ParameterName": pname,
                },
                "NeutrinoFlavourUnosc": [nu_unosc],
                "NeutrinoFlavour": [nu_osc],
                "ParameterBounds": [0.0, 1.0],
                "ParameterGroup": "EParam",
                "ParameterValues": {
                    "Generated": 1.0,
                    "PreFitValue": 1.0,
                },
                "StepScale": {"MCMC": 1.0},
                "Type": "Norm",
            }
        }
    

        systematics.append(block)

# -----------------------------
# Write YAML
# -----------------------------
with open("EParams.yaml", "w") as f:
    yaml.dump(systematics, f, sort_keys=False)

print(f"Saved {len(systematics)} parameters to EParams.yaml")
