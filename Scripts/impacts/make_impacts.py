import json
import ROOT
from read_config import ChainConfig
import metrics

def main(chain_path, metric_type, metric_param=None, save_path="impacts.json", burn_in=0):
  if metric_type in ["Mean"] and not metric_param:
    raise ValueError(f"metric_param must be provided for metric_type {metric_type}")
  
  config = ChainConfig(chain_path)

  rdf_all = ROOT.RDataFrame("posteriors", chain_path)
  ROOT.RDF.Experimental.AddProgressBar(rdf_all)
  rdf = rdf_all.Filter(f"step > {burn_in}")

  rdf_nominal = rdf.Define("weight", "1.0")
  metrics_pointers = {
    "nominal": getattr(metrics, metric_type)(rdf_nominal, metric_param)
  }

  for syst in config.systematics[:20]:
    mu, sigma = syst["PreFitValue"], syst["Error"]
    mu_up = syst["PreFitValue"] + sigma
    mu_down = syst["PreFitValue"] - sigma
    sigma_shrink = 0.5 * sigma
    
    rdf_up = rdf.Define("weight", f"TMath::Gaus({syst["PosteriorsName"]}, {mu_up}, {sigma_shrink}) / TMath::Gaus({syst["PosteriorsName"]}, {mu}, {sigma})")
    rdf_down = rdf.Define("weight", f"TMath::Gaus({syst["PosteriorsName"]}, {mu_down}, {sigma_shrink}) / TMath::Gaus({syst["PosteriorsName"]}, {mu}, {sigma})")
    
    metrics_pointers[syst["FancyName"]+"_up"] = getattr(metrics, metric_type)(rdf_up, metric_param)
    metrics_pointers[syst["FancyName"]+"_down"] = getattr(metrics, metric_type)(rdf_down, metric_param)

  metrics_values = {k: v.GetValue() for k, v in metrics_pointers.items()}

  metric_definition = {
    "type": metric_type,
    "param": metric_param,
  }

  output = {
    "metric_definition": metric_definition,
    "metric_values": {
      "nominal": metrics_values["nominal"],
      "variations": {}
    },
  }
  for syst in config.systematics[:20]:
    output["metric_values"]["variations"][syst["FancyName"]] = {
      "up": metrics_values[syst["FancyName"]+"_up"],
      "down": metrics_values[syst["FancyName"]+"_down"],
    }

  with open(save_path, "w") as f:
    json.dump(output, f, indent=2)

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("chain_path")
  parser.add_argument("metric_type")
  parser.add_argument("--metric-param", "-p", default=None)
  parser.add_argument("--save-path", default="impacts.json")
  parser.add_argument("--burn-in", type=int, default=0)
  args = parser.parse_args()

  ROOT.EnableImplicitMT(16)
  #chain_path = "/home/vol04/scarf1534/Liban/AnaChains/FD_Fit_1000_ktmwyr_dcp_mpi2/HaddedChains/FD_Fit_1000_ktmwyr_dcp_mpi2.root"
  main(args.chain_path, args.metric_type, args.metric_param, args.save_path, args.burn_in)

# model = ROOT.RDF.TH1DModel("", "", 100, -3.14, 0)

# delta_cp_hists = {
#   "nominal": rdf.Histo1D(model, "delta_cp")
# }

# for syst in config.systematics[:5]:
#   mu, sigma = syst["PreFitValue"], syst["Error"]
#   mu_up = syst["PreFitValue"] + sigma
#   mu_down = syst["PreFitValue"] - sigma
#   sigma_shrink = 0.5 * sigma
  
#   rdf = rdf.Define(f"weight_{syst['PosteriorsName']}_up", f"TMath::Gaus({syst['PosteriorsName']}, {mu_up}, {sigma_shrink}) / TMath::Gaus({syst['PosteriorsName']}, {mu}, {sigma})")
#   rdf = rdf.Define(f"weight_{syst['PosteriorsName']}_down", f"TMath::Gaus({syst['PosteriorsName']}, {mu_down}, {sigma_shrink}) / TMath::Gaus({syst['PosteriorsName']}, {mu}, {sigma})")
  
#   delta_cp_hists[syst["FancyName"]+"_up"] = rdf.Histo1D(model, "delta_cp", f"weight_{syst['PosteriorsName']}_up")
#   delta_cp_hists[syst["FancyName"]+"_down"] = rdf.Histo1D(model, "delta_cp", f"weight_{syst['PosteriorsName']}_down")

# delta_cp_means = {k: v.GetMean() for k, v in delta_cp_hists.items()}

# with open("impacts.json", "w") as f:
#   json.dump(delta_cp_means, f, indent=2)

