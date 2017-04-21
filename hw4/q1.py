from pylab import *
from scipy.stats import norm as normal

def get_tau2_posterior_samples_MCMC(deltaq):
	'''
		Get the posterior samples for tau2.
	'''
	N = 20000	
	# Posteriors before selection
	tau2 = zeros(N)
	# Initial guesses
	tau2[0] = rand() * 20 + 10
	# Generating loop
	i = 1
	while i < N:
		# Generate new sample
		tau2_i = rand() * 20 + 10
		# Determine whether to keep
		numer = exp(-deltaq**2 / 2 / (8 + 2 * tau2_i)) / sqrt(8 + 2 * tau2_i)
		denom = exp(-deltaq**2 / 2 / (8 + 2 * tau2[i-1])) / sqrt(8 + 2 * tau2[i-1])
		r = min(numer / denom, 1)
		# Keep with probability r
		if rand() > r: continue
		tau2[i] = tau2_i
		i += 1
	# Discard the first 500 to minimize the influence of the initial guess
	tau2 = tau2[500:]
	return tau2

def main():
	set_printoptions(precision=3)
	P = zeros(6)
	for i, deltaq in enumerate([30, 32, 34, 36, 38, 40]):
		tau2 = get_tau2_posterior_samples_MCMC(deltaq)
		rv = normal()
		numer = 65.2 - deltaq - deltaq * tau2 / (4 + tau2)
		denom = sqrt(8 + 8 * tau2 / (4 + tau2))
		Ps = 1 - rv.cdf(numer / denom)
		P[i] = mean(Ps)
	print P

if __name__ == '__main__':
	main()