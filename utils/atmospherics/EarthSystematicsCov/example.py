from EarthSystematics import *

EarthSys = EarthSystematics("../../build/_deps/cudaprob3-src/modelsPREM_4layer_quad_v2.dat")

# Set 10% uncertainty for layer widths, 5% for density weights
EarthSys.set_sigmas_uniform(0.1, 0.05)

# We don't want a large uncertainty for the last layer, so let's change that
EarthSys.set_sigma('r_4', 10)

# Set the number of pulls to be accepted
EarthSys.set_N_accepted_pulls(100)

# Set the tolerance for the planet's mass and moment of inertia
EarthSys.set_M_tolerance(1e-3)
EarthSys.set_I_tolerance(1e-3)

# Print what we've set
EarthSys.print_systematics()

# Find new models
EarthSys.make_pulls()

# Get correlations, plot and add cov objects to yaml
EarthSys.find_correlations()
EarthSys.display_correlation_matrix()
EarthSys.plot_correlation_matrix()
EarthSys.add_osccov_to_yaml('../../configs/CovObjs/OscCov_PDG2021_v2_Atmospherics.yaml')