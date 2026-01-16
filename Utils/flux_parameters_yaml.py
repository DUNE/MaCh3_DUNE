import yaml

# ============================
#  Insert your parameter list
# ============================
flux_params = [
    "TargetUpstreamDegredation",
    "TargetTiltTransverseY",
    "TargetTiltTransverseX",
    "TargetLength",
    "TargetDisplaceTransverseY",
    "TargetDisplaceTransverseX",
    "TargetDensity",
    "ProtonBeamTransverseY",
    "ProtonBeamTransverseX",
    "ProtonBeamRadius",
    "ProtonBeamAngleY",
    "ProtonBeamAngleX",
    "HornWaterLayerThickness",
    "HornCurrent",
    "HornCTiltTransverseY",
    "HornCTiltTransverseX",
    "HornCEllipticityXInducedBField",
    "HornCEccentricityXInducedBField",
    "HornCDisplaceLongitudinalZ",
    "HornBTiltTransverseY",
    "HornBTiltTransverseX",
    "HornBEllipticityXInducedBField",
    "HornBDisplaceLongitudinalZ",
    "HornATiltTransverseY",
    "HornATiltTransverseX",
    "HornAEllipticityXInducedBField",
    "HornAEccentricityXInducedBField",
    "HornADisplaceLongitudinalZ",
    "DecayPipeTiltY",
    "DecayPipeTiltX",
    "DecayPipeRadius",
    "DecayPipeLength",
    "DecayPipeGeoBField",
    "DecayPipeEllipticalCrossSectionYB",
    "DecayPipeEllipticalCrossSectionXA",
    "DecayPipeDisplaceTransverseY",
    "DecayPipeDisplaceTransverseX",
    "DecayPipe3SegmentBowingY",
    "DecayPipe3SegmentBowingX",
    "HornCDisplaceTransverseY",
    "HornBDisplaceTransverseY",
    "HornADisplaceTransverseY",
    "HornCDisplaceTransverseX",
    "HornBDisplaceTransverseX",
    "HornADisplaceTransverseX",
] + [
    f"Flux_HadProd_Param_{i}" for i in range(21)
]

# ============================
# YAML entry template
# ============================

def make_block(name):
    return {
        "Systematic": {
            "Sample_Name": ["ND*"],
            "Error": 0.1,
            "FlatPrior": False,
            "Names": {
                "FancyName": name,
                "ParameterName": name
            },
            "ParameterBounds": [-1.0, 5.0],
            "ParameterGroup": "FluxSys",
            "ParameterValues": {
                "Generated": 0.0,
                "PreFitValue": 0.0
            },
            "StepScale": {
                "MCMC": 0.05
            },
            "Type": "Functional"
        }
    }

# ============================
# Build YAML structure
# ============================

yaml_dict = {"Systematics": [make_block(p) for p in flux_params]}

# ============================
# Write YAML file
# ============================

with open("flux_systematics.yaml", "w") as f:
    yaml.dump(yaml_dict, f, sort_keys=False)

print("Generated flux_systematics.yaml with", len(flux_params), "parameters.")
