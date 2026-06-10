import ROOT

class HistBasedMetric:
  def __init__(self, rdf, param_name, model=None):
    model = ROOT.RDF.TH1DModel("", "", 100, -3.14, 0) # hard-coded for deltacp for the moment
    if model:
      self.hist = rdf.Histo1D(model, param_name, "weight")
    else:
      self.hist = rdf.Histo1D(param_name, "weight")

class Mean(HistBasedMetric):
  def GetValue(self):
    return self.hist.GetMean()

# have to be careful with size of bins otherwise it looks like most systematics have no impact
# class HPD(HistBasedMetric):    
#   def GetValue(self):
#     maxBin = self.hist.GetMaximumBin()
#     return self.hist.GetBinCenter(maxBin)
    

