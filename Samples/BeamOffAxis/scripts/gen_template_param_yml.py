#!/usr/bin/env python

enu_bins = [ 0.0, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5,
             3.75, 4.0, 5.0, 6.0, 10.0, 100.0 ]
enubias_bins = [ -2.0, -0.6, -0.581, -0.5595, -0.538, -0.5165, -0.495, -0.4735,
                 -0.452, -0.4305, -0.409, -0.3875, -0.366, -0.3445, -0.323,
                 -0.3015, -0.28, -0.2585, -0.237, -0.2155, -0.194, -0.1725,
                 -0.151, -0.1295, -0.108, -0.0865, -0.065, -0.0435, -0.022,
                 0.0, 0.1 ]

print("""---
Systematics:""")

for ie in range(len(enu_bins)-1):
  for ieb in range(len(enubias_bins)-1):

    pnm = f"tmplw_e{ie}_eb{ieb}"

    print(f"""- Systematic:
    SampleNames: ["*"]
    Error:  0.5
    FlatPrior: false
    KinematicCuts:
    - TrueNeutrinoEnergy:
      - {enu_bins[ie]}
      - {enu_bins[ie+1]}
    - Enubias:
      - {enubias_bins[ieb]}
      - {enubias_bins[ieb+1]}
    Names:
      FancyName: {pnm}
    ParameterBounds:
    - 0.0
    - 4.0
    ParameterValues:
      Generated: 1.0
      PreFitValue: 1.0
    StepScale:
      MCMC: 1.0
    Type: Norm
    ParameterGroup: Xsec""")
