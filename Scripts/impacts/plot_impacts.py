import matplotlib.pyplot as plt
import mplhep
mplhep.set_style("CMS")
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import json

def get_impacts_fig(impacts, pulls, num_offset=0):
  systs = list(pulls.keys())
  variations = impacts["metric_values"]["variations"]
  
  fig, axs = plt.subplots(1, 3, sharey=True, layout="constrained")
  fig.set_constrained_layout_pads(wspace=0, hspace=0, w_pad=0, h_pad=0)
  fig.get_layout_engine().set(
    rect=(0.01, 0.01, 0.98, 0.98)
  )
  axs[0].set_xticks([])
  axs[0].set_yticks([])
  axs[0].set_ylim(-0.5, len(systs)-0.5)
  axs[0].set_xlim(0, 1)
  axs[1].set_xlim(-3, 3)
  axs[1].set_xticks([-2, -1, 0, 1, 2])
  axs[0].yaxis.set_inverted(True)
  axs[0].set_axis_off()
  axs[1].set_xlabel(r"$(\nu_{\text{fit}} - \nu_{\text{prior}}) / \sigma_{\text{prior}}$")

  

  nominal_metric = impacts["metric_values"]["nominal"]

  min_metric = min(*[v["down"] for v in variations.values()], *[v["up"] for v in variations.values()])
  max_metric = max(*[v["down"] for v in variations.values()], *[v["up"] for v in variations.values()])
  
  metric_range = max(max_metric - nominal_metric,  nominal_metric - min_metric)
  axs[2].set_xlim(-1.1*metric_range, 1.1*metric_range)
  
  #axs[2].set_xlabel(r"$\Delta \delta_{\text{CP}}$")
  metric_label = impacts["metric_definition"]["type"]
  if impacts["metric_definition"]["param"] is not None:
    metric_label += f"({impacts['metric_definition']['param']})"
  axs[2].set_xlabel(rf"$\Delta${metric_label}", fontsize=18)
  
  #axs[2].text(0, -0.5, r"$\hat{\delta}_{\text{CP}} = %0.2f$" % nominal_metric, ha="center", va="bottom", fontsize=18)
  axs[2].text(0, -0.5, rf"{metric_label}$ = %0.2f$" % nominal_metric, ha="center", va="bottom", fontsize=18)
  
  for x in [-2, -1, 0, 1, 2]:
    axs[1].axvline(x, color="lightgray", linestyle="--")

  # for x in [-metric_range*0.75, 0, metric_range]:
  #   axs[2].axvline(x, color="lightgray", linestyle="--")
  axs[2].axvline(0, color="lightgray", linestyle="--")

  for i, syst in enumerate(systs):
    print(syst)
    if i % 2 == 1:
      axs[0].fill_between([0, 1], i-0.5, i+0.5, color="lightgray", alpha=0.5)
      axs[1].fill_between([-3, 3], i-0.5, i+0.5, color="lightgray", alpha=0.5)
      axs[2].fill_between([-1.5*metric_range, 1.5*metric_range], i-0.5, i+0.5, color="lightgray", alpha=0.5)
      
    #name_str = f"{impact['FancyName']}\n({impact['name']})"
    name_str = syst
    axs[0].text(0.1, i, str(i+num_offset+1), ha="left", va="center", fontweight="bold", fontsize=16)
    axs[0].text(0.9, i, name_str, ha="right", va="center", fontsize=12)  
    
    pull, posterior_std = pulls[syst]
    axs[1].errorbar(pull, i, xerr=posterior_std, fmt="o", capsize=5, color="black")

    up_label = "Up" if i == 0 else ""
    down_label = "Down" if i == 0 else ""
    axs[2].barh(i, variations[syst]["up"]-nominal_metric, height=0.4, color="tab:blue", alpha=0.5, label=up_label)
    axs[2].barh(i, variations[syst]["down"]-nominal_metric, height=0.4, color="tab:red", alpha=0.5, label=down_label)

  fig.legend(loc="outside lower left", ncol=2, fontsize=18)
  
  return fig
  
def make_impact_plot(impacts, pulls, save_path=None):
  # find ordering according to impact size
  nom = impacts["metric_values"]["nominal"]
  variations = impacts["metric_values"]["variations"]
  impact_size = [max(abs(v["up"]-nom), abs(v["down"]-nom)) for v in variations.values()]
  ordering = np.argsort(impact_size)[::-1]
  systs = list(pulls.keys())
  
  with PdfPages(save_path) as pdf:
    for i in range(0, len(systs), 15):
      pulls_subset = {systs[j]: pulls[systs[j]] for j in ordering[i:i+15]}
      impacts_subset = {
        "metric_definition": impacts["metric_definition"],
        "metric_values": {
          "nominal": impacts["metric_values"]["nominal"],
          "variations": {systs[j]: variations[systs[j]] for j in ordering[i:i+15]}
        }
      }
      
      fig = get_impacts_fig(impacts_subset, pulls_subset, num_offset=i)
      pdf.savefig(fig)
      plt.close(fig)    
      
def main(impacts_path, pull_paths, save_path="impacts.pdf"):
  with open(impacts_path, "r") as f:
    impacts = json.load(f)
  with open(pull_paths) as f:
    pulls = json.load(f)
    
  assert set(impacts["metric_values"]["variations"].keys()) == set(pulls.keys()), "Impacts and pulls must have the same systematics"
  
  make_impact_plot(impacts, pulls, save_path)
  
if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("impacts_path")
  parser.add_argument("pull_paths")
  parser.add_argument("--save-path", default="impacts.pdf")
  args = parser.parse_args()

  main(args.impacts_path, args.pull_paths, args.save_path)