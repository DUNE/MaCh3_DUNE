import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ---------------------------------------------------------------
# Classic academic style: serif fonts, LaTeX-like math, muted look
# ---------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["cmr10", "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True,
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "legend.fontsize": 10.5,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "figure.dpi": 150,
})

# ---------------------------------------------------------------
# Physics constants
# ---------------------------------------------------------------
m_pi = 0.13957      # GeV, charged pion mass
m_mu = 0.10566      # GeV, muon mass

# Neutrino energy in the pion rest frame (two-body decay kinematics)
E_nu_star = (m_pi**2 - m_mu**2) / (2 * m_pi)   # ~0.02979 GeV

L_ND = 574.0  # m, distance from target to DUNE Near Detector hall

# ---------------------------------------------------------------
# Actual DUNE-PRISM ND off-axis positions (metres)
# ---------------------------------------------------------------
off_axis_positions = [
    0.00, 1.75, 2.00, 4.00, 5.75, 8.00, 9.75, 12.00, 13.75,
    16.00, 17.75, 20.00, 21.75, 24.00, 25.75, 26.25, 28.00, 28.25,
]

# Pion energy range (lab frame)
E_pi = np.linspace(m_pi * 1.01, 6.0, 500)  # GeV
gamma_pi = E_pi / m_pi

# ---------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.8, 4.8))

cmap = plt.get_cmap("viridis")
norm = mpl.colors.Normalize(vmin=min(off_axis_positions), vmax=max(off_axis_positions))
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)

for x_off in off_axis_positions:
    theta = np.arctan(x_off / L_ND)  # radians
    E_nu_prime = (2 * gamma_pi * E_nu_star) / (1 + (gamma_pi * theta) ** 2)
    ax.plot(E_pi, E_nu_prime, color=cmap(norm(x_off)), lw=1.6)

ax.set_xlabel(r"Pion energy, $E_{\pi}$ [GeV]")
ax.set_ylabel(r"Neutrino energy, $E'_{\nu}$ [GeV]")
ax.set_xlim(0, 6.0)
ax.set_ylim(0, None)

ax.set_title("DUNE-PRISM: off-axis neutrino energy vs. parent pion energy", pad=10)

cbar = fig.colorbar(sm, ax=ax, pad=0.02)
cbar.set_label(r"Off-axis position [m]")

ax.grid(alpha=0.25, linewidth=0.6)
fig.tight_layout()

fig.savefig("/mnt/user-data/outputs/dune_prism_enu_vs_epi.png", dpi=300, bbox_inches="tight")
fig.savefig("/mnt/user-data/outputs/dune_prism_enu_vs_epi.pdf", bbox_inches="tight")

print("Saved.")