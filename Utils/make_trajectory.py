import numpy as np

np.random.seed(42)
mu = np.array([1.0, 1.0])
sigma1, sigma2, rho = 0.15, 0.15, 0.6
cov = np.array([[sigma1**2, rho*sigma1*sigma2],
                 [rho*sigma1*sigma2, sigma2**2]])
cov_inv = np.linalg.inv(cov)

def neg_log_L(theta):
    d = theta - mu
    return 0.5 * d @ cov_inv @ d

# start deliberately off-mode to show burn-in
theta = np.array([0.4, 1.6])
step_scale = 0.08
samples = [theta.copy()]

for _ in range(80):
    proposal = theta + np.random.normal(0, step_scale, size=2)
    alpha = min(1, np.exp(neg_log_L(theta) - neg_log_L(proposal)))
    if np.random.uniform() < alpha:
        theta = proposal
    samples.append(theta.copy())

samples = np.array(samples)
np.savetxt("trajectory.dat", samples, fmt="%.4f")