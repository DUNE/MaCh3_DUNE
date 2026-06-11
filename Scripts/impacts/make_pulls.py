import ROOT
import json
from read_config import ChainConfig

def main(chain_path, save_path="pulls.json", burn_in=0):
  config = ChainConfig(chain_path)

  rdf_all = ROOT.RDataFrame("posteriors", chain_path)
  ROOT.RDF.Experimental.AddProgressBar(rdf_all)
  rdf = rdf_all.Filter(f"step > {burn_in}")

  results_pointers = {}

  for syst in config.systematics[:]:
    results_pointers[syst["FancyName"]] = [rdf.Mean(syst["PosteriorsName"]), rdf.StdDev(syst["PosteriorsName"])]
    
  pulls = {}
  for syst in config.systematics[:]:
    mean = results_pointers[syst["FancyName"]][0].GetValue()
    stddev = results_pointers[syst["FancyName"]][1].GetValue()
    pull = (mean - syst["PreFitValue"]) / syst["Error"]
    pulls[syst["FancyName"]] = [pull, stddev / syst["Error"]]
    
  with open(save_path, "w") as f:
    json.dump(pulls, f, indent=2)
  
if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("chain_path")
  parser.add_argument("--save-path", default="pulls.json")
  parser.add_argument("--burn-in", type=int, default=0)
  args = parser.parse_args()

  ROOT.EnableImplicitMT()
  #chain_path = "/home/vol04/scarf1534/Liban/AnaChains/FD_Fit_1000_ktmwyr_dcp_mpi2/HaddedChains/FD_Fit_1000_ktmwyr_dcp_mpi2.root"
  main(args.chain_path, args.save_path, args.burn_in)