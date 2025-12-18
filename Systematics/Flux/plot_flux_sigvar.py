import uproot
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np

Systematics = []

with uproot.open("flux_variations_FD_and_PRISM_2023.root") as inpfile:
  with uproot.open("flux_1sigma_weights.root") as valifile:
    with uproot.open(sys.argv[1]) as varfile:

      for pgk in inpfile["FluxParameters"].keys(recursive=False):
        pg = pgk.split(";")[0]
        for pk in inpfile["FluxParameters"][pg].keys(recursive=False):
          p = pk.split(";")[0]
          pdir = inpfile["FluxParameters"][pg][p]
          offaxis_axis = pdir["OffAxisTAxis"] \
                            if (pg == "Focussing") else \
                            pdir["ND_nu_numu"]["OffAxisTAxis"]
          
          num_hists = np.count_nonzero([ x for x in pdir["ND_nu_numu"].keys(recursive=False) if (x.find("ND_nu_numu") == 0) ])
          hists = [ pdir["ND_nu_numu"][f"ND_nu_numu_{x}"] for x in range(num_hists) ]
          Systematics.append({ "name": p, 
                               "offaxis_axis": offaxis_axis, 
                               "hists": hists})


      Systematics.reverse()
      pp = PdfPages('FluxVariationsValidation.pdf')
      for syst in Systematics:
        print(syst["name"])

        vardir_0m = varfile[syst["name"]]["ND_FHC_CCnumu_0m"]
        vardir_12m = varfile[syst["name"]]["ND_FHC_CCnumu_12m"]

        fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(18, 12))

        fig.text(0.5,0.95,"Varied Parameter: " + syst["name"], horizontalalignment="center", size="xx-large")

        labels = [r"$-3\sigma$",r"$-1\sigma$",r"Nominal",r"$+1\sigma$",r"$+3\sigma$"]
        chweel = ["firebrick", "tomato", "black", "deepskyblue", "mediumblue"]
        lswheel = ["solid", "dashed", "dotted", "dashdot"]
        for i in range(5):
          bin_heights, bin_edges = vardir_0m[f"Variation_{i}"].to_numpy()
          axs[0,0].stairs(bin_heights, bin_edges, label=labels[i], color=chweel[i])

          if i != 2:
            bin_heights_ref, bin_edges_ref = vardir_0m["Variation_2"].to_numpy()
            axs[1,0].stairs(bin_heights/bin_heights_ref, bin_edges, color=chweel[i])

          bin_heights, bin_edges = vardir_12m[f"Variation_{i}"].to_numpy()
          axs[0,1].stairs(bin_heights, bin_edges, label=labels[i], color=chweel[i])

          if i != 2:
            bin_heights_ref, bin_edges_ref = vardir_12m["Variation_2"].to_numpy()
            axs[1,1].stairs(bin_heights/bin_heights_ref, bin_edges, color=chweel[i])

        for i,oad in enumerate([0, 12]):
          off_axis_edges = syst["offaxis_axis"].edges()
          oa_bin = np.argmax(off_axis_edges >= oad)
          bin_heights, bin_edges = syst["hists"][oa_bin].to_numpy()
          enu_lim = np.argmax(bin_edges > 5)
          axs[0,2].stairs(bin_heights[:enu_lim-1], bin_edges[:enu_lim], color=chweel[3], \
            linestyle=lswheel[i], label=f"${off_axis_edges[oa_bin]} -- {off_axis_edges[oa_bin+1]}$ m")

          bin_heights, bin_edges = valifile[syst["name"]][f"ND_Numu_{oad}m"].to_numpy()
          axs[1,2].stairs(bin_heights, bin_edges, color=chweel[3], \
            linestyle=lswheel[i], label=f"{oad} m")

        axs[1,0].set_ylim([0.8,1.2])

        axs[0,2].set_ylim([-0.2, 0.2])
        axs[1,1].set_ylim([0.8,1.2])
        axs[1,2].set_ylim([0.8,1.2])

        axs[0,0].set_xlim([0,5])
        axs[1,0].set_xlim([0,5])

        axs[0,1].set_xlim([0,5])
        axs[1,1].set_xlim([0,5])

        axs[0,2].set_xlim([0,5])
        axs[1,2].set_xlim([0,5])

        fig.text(0.2,0.89,"MaCh3 SigVar On Axis",horizontalalignment="left", size="large")
        fig.text(0.2,0.47,"MaCh3 SigVar On Axis Ratio",horizontalalignment="left", size="large")

        fig.text(0.4,0.89,"MaCh3 SigVar 12 m Off Axis",horizontalalignment="left", size="large")
        fig.text(0.4,0.47,"MaCh3 SigVar 12 m Off Axis Ratio",horizontalalignment="left", size="large")

        fig.text(0.65,0.89,"Histograms from inputs",horizontalalignment="left", size="large")
        fig.text(0.65,0.47,"Weights directly from Weight Interface",horizontalalignment="left", size="large")

        axs[0,0].legend()
        axs[0,2].legend()
        axs[1,2].legend()

        axs[0,0].set_ylabel("Event Rate", size="x-large")
        axs[1,0].set_ylabel("Varied/Nominal Event Rate", size="x-large")
        axs[0,1].set_ylabel("Event Rate", size="x-large")
        axs[1,1].set_ylabel("Varied/Nominal Event Rate", size="x-large")
        
        axs[0,2].set_ylabel("Input Ratio", size="x-large")
        axs[1,2].set_ylabel("Systematic Weight", size="x-large")

        axs[0,0].set_xlabel(r"$E_\nu$", size="x-large")
        axs[1,0].set_xlabel(r"$E_\nu$", size="x-large")
        axs[0,1].set_xlabel(r"$E_\nu$", size="x-large")
        axs[1,1].set_xlabel(r"$E_\nu$", size="x-large")
        axs[0,2].set_xlabel(r"$E_\nu$", size="x-large")
        axs[1,2].set_xlabel(r"$E_\nu$", size="x-large")
        pp.savefig(fig)
        plt.close(fig)
      pp.close()