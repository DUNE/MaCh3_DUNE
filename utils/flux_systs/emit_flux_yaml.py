# ---
# Systematics:
# - Systematic:
#     SampleNames: ["ND_*", "FD_*"]
#     Error: 1.0
#     FlatPrior: false
#     Names:
#       FancyName: MAQE
#       ParameterName: MAQE
#     ParameterBounds:
#     - -4.0
#     - 4.0
#     ParameterGroup: Xsec
#     ParameterValues:
#       Generated: 0.0
#       PreFitValue: 0.0
#     SplineInformation:
#       Mode:
#       - 0
#       SplineName: maqe
#     StepScale:
#       MCMC: 0.002
#     Type: Spline

import yaml

import ROOT

Systematics = []

fin = ROOT.TFile.Open("flux_variations_FD_and_PRISM_2023.root")

dfp = fin.GetDirectory("FluxParameters")
for pgroups in dfp.GetListOfKeys():
  for p in dfp.GetDirectory(pgroups.GetName()).GetListOfKeys():
    Systematic = {}
    Systematic["SampleNames"] = ["*"]
    Systematic["Error"] = 1.0
    Systematic["FlatPrior"] = False
    Systematic["Names"] = { "FancyName": p.GetName(), "ParameterName": p.GetName() }
    Systematic["ParameterBounds"] = [-4.0, 4.0]
    Systematic["ParameterGroup"] = "Flux"
    Systematic["ParameterValues"] = { "Generated": 0, "PreFitValue": 0 }
    Systematic["StepScale"] = { "MCMC": 0.01 }
    Systematic["Type"] = "Functional"
    Systematics.append({"Systematic": Systematic})

print(yaml.dump({"Systematics": Systematics}))