from matplotlib.backends.backend_pdf import PdfPages
import uproot
import numpy as np
import matplotlib.pyplot as plt
import math

# ---------------------------
# USER SETTINGS
# ---------------------------
root_file = "/scratch/abipeake/October_chains/Tuesday_041125_enubuas_postadaptivetestnew/Tuesday_041125_enubuas_postadaptivetestnew_chain_6_job_0.root" #"/scratch/abipeake/October_chains/Tuesday_041125_enubias_newattempt_withinitialproposal/Tuesday_041125_enubias_newattempt_withinitialproposal_chain_0_job_0copy2.root"
tree_name = "posteriors"
n_params = 200
bin_size = 1000
prior_mean_val = 1.0
prior_cv_val = 0.5
params_per_page = 4
output_pdf = "MCMC_convergence_new.pdf"  # Single PDF for all pages
# ---------------------------

# Generate parameter names
param_names = [f"xsec_{i}" for i in range(n_params)]
prior_mean = [prior_mean_val]*n_params
prior_cv = [prior_cv_val]*n_params

# Load ROOT tree
file = uproot.open(root_file)
tree = file[tree_name]

# Extract samples
samples = [tree[name].array(library="np") for name in param_names]
samples = np.column_stack(samples)
n_steps, n_params = samples.shape

# Compute cumulative mean and CV
n_bins = n_steps // bin_size
x_vals = np.arange(n_bins) * bin_size

binned_mean = np.zeros((n_bins, n_params))
binned_cv = np.zeros((n_bins, n_params))

for j in range(n_params):
    for k in range(n_bins):
        # Take all steps up to the end of this bin
        chunk = samples[:(k+1)*bin_size, j]
        mean = np.mean(chunk)
        std = np.std(chunk)
        binned_mean[k, j] = mean
        binned_cv[k, j] = std / abs(mean) if abs(mean) > 1e-12 else np.nan

# Number of pages
n_pages = math.ceil(n_params / params_per_page)

# Use PdfPages to save all pages into one PDF
with PdfPages(output_pdf) as pdf:
    for page in range(n_pages):
        fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
        axes = axes.flatten()
        
        for i in range(params_per_page):
            j = page*params_per_page + i
            if j >= n_params:
                axes[i].axis('off')
                continue
            
            ax = axes[i]
            ax.plot(x_vals, binned_mean[:, j], label='Binned mean', color='blue')
            ax.plot(x_vals, binned_cv[:, j], label='Binned CV', color='orange')
            ax.axhline(prior_mean[j], color='red', linestyle='--', label='Prior mean')
            ax.fill_between(x_vals,
                            prior_mean[j]*(1-prior_cv[j]),
                            prior_mean[j]*(1+prior_cv[j]),
                            color='green', alpha=0.2, label='Prior CV band')
            ax.set_title(param_names[j])
            ax.grid(True)
            ax.set_xticks(x_vals[::max(1, n_bins//10)])
            ax.set_xticklabels([f"{int(val)}" for val in x_vals[::max(1, n_bins//10)]], rotation=45)
        
        # Single legend per page
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper right', fontsize=10)
        
        plt.xlabel('Iteration')
        plt.tight_layout(rect=[0, 0, 0.95, 0.97])
        plt.suptitle(f'MCMC Diagnostics: Page {page+1}', fontsize=14)
        
        # Save this figure as a page in the PDF
        pdf.savefig(fig)
        plt.close(fig)

print(f"Saved all pages in a single PDF: {output_pdf}")
