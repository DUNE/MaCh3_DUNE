import ROOT
import numpy as np
import matplotlib.pyplot as plt

# --- File path ---
filename = "/scratch/abipeake/October_chains/EventRates_1D_70bins_afteradapted_wed/poteriorpredictivenew.root"

# --- Open ROOT file ---
f = ROOT.TFile.Open(filename)
if not f or f.IsZombie():
    raise RuntimeError(f"Cannot open file: {filename}")

# --- ROOT 6.26-safe way to get all keys ---
all_keys = []
keys_list = f.GetListOfKeys()
for i in range(keys_list.GetEntries()):  # Use GetEntries(), not GetSize()
    key = keys_list.At(i)                 # TKey object
    all_keys.append(key.GetName())

print(f"Found {len(all_keys)} objects in the ROOT file")
print(all_keys[:10])

# --- Filter histograms ending with "_pdf0" ---
hist_names = [name for name in all_keys if name.endswith("_pdf0")]
if not hist_names:
    raise RuntimeError("No histograms ending with '_pdf0' found in the file.")

print(f"Found {len(hist_names)} histograms. First few:")
print(hist_names[:5])

# --- Load histograms ---
samples = [f.Get(name) for name in hist_names]

# --- Prepare binning ---
nBins = samples[0].GetNbinsX()
bin_centers = np.array([samples[0].GetBinCenter(b) for b in range(1, nBins + 1)])

# --- Build array of bin contents ---
values = np.array([[h.GetBinContent(b) for h in samples] for b in range(1, nBins + 1)])

# --- Compute mean and standard deviation ---
mean = np.mean(values, axis=1)
std = np.std(values, axis=1)

# --- Plot ---
plt.figure(figsize=(8,5))
plt.fill_between(bin_centers, mean - std, mean + std, color='skyblue', alpha=0.5, label='±1σ')
plt.plot(bin_centers, mean, color='blue', lw=2, label='Mean')

plt.xlabel('Reco neutrino energy [GeV]')
plt.ylabel('Events')
plt.title('Posterior predictive mean ± 1σ')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

# --- Save figure ---
output_file = "/scratch/abipeake/October_chains/EventRates_1D_70bins_afteradapted_wed/posterior_mean_std.png"
plt.savefig(output_file, dpi=300)
print(f"Plot saved to {output_file}")

plt.show()
