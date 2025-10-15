import uproot
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

Systematics = []

with uproot.open("flux_variations_FD_and_PRISM_2023.root") as file:
  for pgk in file["FluxParameters"].keys(recursive=False):
    pg = pgk.split(";")[0]
    for pk in file["FluxParameters"][pg].keys(recursive=False):
      p = pk.split(";")[0]
      Systematics.append(p)

with uproot.open(sys.argv[1]) as file:

  pp = PdfPages('bla.pdf')
  for syst in Systematics:

    fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(6, 12))
    fig.text(0.5,0.95,syst)
    for i in range(5):
      bin_heights, bin_edges = file[syst]["ND_FHC_CCnumu"][f"Variation_{i}"].to_numpy()
      axs[0].stairs(bin_heights, bin_edges)

      if i != 2:
        bin_heights_ref, bin_edges_ref = file[syst]["ND_FHC_CCnumu"]["Variation_2"].to_numpy()
        axs[1].stairs(bin_heights/bin_heights_ref, bin_edges)

    axs[1].set_ylim([0.8,1.2])
    pp.savefig(fig)
    fig.clf()
  pp.close()