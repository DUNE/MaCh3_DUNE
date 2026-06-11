import ROOT
import numpy as np
import warnings
import scipy.stats as stats

class HistBasedMetric:
  model = None
  
  def __init__(self, rdf, param_name):
    if HistBasedMetric.model is None:
      print("Initializing model for HistBasedMetric")
      stat = rdf.Stats(param_name)
      HistBasedMetric.model = ROOT.RDF.TH1DModel("", "", 10000, stat.GetMin(), stat.GetMax())
    self.hist = rdf.Histo1D(HistBasedMetric.model, param_name, "weight")
    
  def getQuantile(self, q):
    xp = np.array([0], dtype=np.double)
    p = np.array([q], dtype=np.double)
    self.hist.GetQuantiles(1, xp, p)
    return xp[0]

class Mean(HistBasedMetric):
  def GetValue(self):
    return self.hist.GetMean()

# have to be careful with size of bins otherwise it looks like most systematics have no or very similar impact
class HPD(HistBasedMetric):
  def __init__(self, rdf, param_name):
    warnings.warn("HPD metric is sensitive to the number of bins in the histogram. It can lead to unreliable ordering of impacts. Use with caution.")
    super().__init__(rdf, param_name)
      
  def GetValue(self):
    #self.hist.Smooth(1) # might want to smooth
    maxBin = self.hist.GetMaximumBin()
    return self.hist.GetBinCenter(maxBin)

class Median(HistBasedMetric):
  def GetValue(self):
    return self.getQuantile(0.5)
  
class StdDev(HistBasedMetric):
  def GetValue(self):
    return self.hist.GetStdDev()
  
class IntervalWidth68(HistBasedMetric):
  def GetValue(self):
    return self.getQuantile(0.84) - self.getQuantile(0.16)

class IntervalMid68(HistBasedMetric):
  def GetValue(self):
    return 0.5 * (self.getQuantile(0.84) + self.getQuantile(0.16))
  
class IntervalWidth95(HistBasedMetric):
  def GetValue(self):
    return self.getQuantile(0.975) - self.getQuantile(0.025)
  
class IntervalMid95(HistBasedMetric):
  def GetValue(self):
    return 0.5 * (self.getQuantile(0.975) + self.getQuantile(0.025))

class MassOrderingBF:
  def __init__(self, rdf, metric_param=None, rw_factor=(1/np.exp(10))):
    self.n_io = rdf.Filter("delm2_23 < 0").Sum("weight")
    self.n_no = rdf.Filter("delm2_23 > 0").Sum("weight")
    self.rw_factor = rw_factor # hard-coded rw factor that accounts for Hank's upweighting of the IO region in his fits
    
  def GetValue(self):
    return self.rw_factor * self.n_io.GetValue() / self.n_no.GetValue()
  
class MassOrderingBFSigma(MassOrderingBF):
  def GetValue(self):
    bayes_factor = super().GetValue()
    return stats.norm.ppf(1 - 0.5 * min(bayes_factor, 1/bayes_factor))
    