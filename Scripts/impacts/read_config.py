import ROOT
import yaml

class Version:
  def __init__(self, version_str):
    self.version_str = version_str
    self.major, self.minor, self.patch = [int(x) for x in version_str.split(".")]
  
  def __gt__(self, other):
    if self.major != other.major:
      return self.major > other.major
    if self.minor != other.minor:
      return self.minor > other.minor
    return self.patch > other.patch
  
  def __ge__(self, other):
    return self > other or self == other
  
  def __str__(self):
    return self.version_str

class ChainConfig:
  def __init__(self, file_path):
    f = ROOT.TFile(file_path, "READ")
    self.config = self.__TMacroToYAML(f.Get("MaCh3_Config"))
    self.mach3_version = Version(self.__getMaCh3Version(f))
    self.xsec_cov = self.__TMacroToYAML(f.Get("CovarianceFolder/Config_xsec_cov"))
    self.systematics = self.__getSystematics()
    f.Close()    

  def __getSystematics(self):
    systematics = []
    
    posteriors_param_prefix = "param_" if self.mach3_version >= Version("2.4.0") else "xsec_"
    
    for i, param in enumerate(self.xsec_cov["Systematics"]):
      if param["Systematic"]["ParameterGroup"] == "Osc":
        continue
      systematics.append({
        "PosteriorsName": f"{posteriors_param_prefix}{i}",
        "ParameterName": param["Systematic"]["Names"]["ParameterName"],
        "FancyName": param["Systematic"]["Names"]["FancyName"],
        "Error": param["Systematic"]["Error"],
        "PreFitValue": param["Systematic"]["ParameterValues"]["PreFitValue"],
        "GeneratedValue": param["Systematic"]["ParameterValues"]["Generated"],
      })
    return systematics

  def __getMaCh3Version(self, f):
    line = str(f.Get("MaCh3Engine/version_header").GetLineWith("const char* MaCh3_VERSION="))
    version = line.split("=")[-1].split(";")[0].strip('"')
    return version

  def __TMacroToYAML(self, macro):
    linesList = macro.GetListOfLines()
    s = "\n".join([str(line) for line in linesList])
    return yaml.safe_load(s)

if __name__ == "__main__":
  chain_path = "/work4/ppd/scarf1488/mach3_stuff/MaCh3_DUNE/outputs/runs/180_fit/180_fit_hadded_300.root"
  #chain_path = "~/MaCh3Things/MaCh3_DUNE_FD_Fit/FD_Fit_downsampled_use_adaptive_test_2.root"
  config = ChainConfig(chain_path)