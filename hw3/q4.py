from pylab import *
from scipy.stats import norm as normal

def get_gamma_pdf(x, alpha, beta):
	'''
		Return the pdf of custom gamma function on given points.
	'''
	return 1.0 / (beta**alpha * gamma(alpha)) * x**(alpha-1) * exp(-x/beta)

def gamma_rvs(alpha, beta, N=1):
	'''
		Return samples of the custom inverse Gamma RV.
	'''
	r = arange(0,300.1,0.1)
	gamma_pdf = get_gamma_pdf(r, alpha, beta)
	gamma_pdf = gamma_pdf / sum(gamma_pdf)
	if N == 1:
		return choice(r, p=gamma_pdf)
	else:
		return choice(r, size=N, p=gamma_pdf)

def get_beta_tau2_posterior_samples_MCMC(X, sigma2, n):
	'''
		Get the posterior samples for beta and tau2
	'''
	N = 20000
	# Posteriors before selection
	beta_j = zeros(N)
	tau2_j = zeros(N)
	# Initial guesses
	beta_j[0] = randn()*2 + 12
	tau2_j[0] = gamma_rvs(2.5, 2)
	# Generating loop
	j = 1
	while j < N:
		# Generate new sample
		beta = gamma_rvs(2.5, 2)
		tau2 = randn()*2 + 12
		# Determine whether to keep
		numerator = (randn()*sqrt(sigma2/n+tau2) + beta) * (randn()*sqrt(sigma2/n+tau2) + beta)
		denominator = (randn()*sqrt(sigma2/n+tau2_j[j-1]) + beta_j[j-1]) * (randn()*sqrt(sigma2/n+tau2_j[j-1]) + beta_j[j-1])
		r = min(numerator/denominator, 1)
		# Keep with probability r
		if rand() > r: continue
		beta_j[j] = beta
		tau2_j[j] = tau2
		j += 1
	# Discard the first 500 to minimize the influence of the initial guess
	beta_j = beta_j[500:]
	tau2_j = tau2_j[500:]
	return beta_j, tau2_j

def main():
	set_printoptions(precision=3)
	X = [array([12.7,6.6,14.7,12.2,4.4,7.8,13.8,13.7,11.1,9.1,14.0]),
		array([17.1,11.9,12.7,16.8,15.0,14.6,13.7,16.4]),
		array([5.2,4.5,10.5,15.0,5.0,14.9,7.6,8.3,10.8,14.6,15.1,7.0,9.3]),
		array([14.3,16.2,10.0,13.1,16.9,11.2,10.1,18.3,13.5,15.0,15.1,14.8,15.7,13.2,12.2,13.2]),
		array([10.5,7.5,4.7,12.5,13.1,13.5,12.2,16.1,9.0,17.9])]
	sigma2 = 4.
	n = 8.
	# Get the posteriors
	beta_post, tau2_post = get_beta_tau2_posterior_samples_MCMC(X, sigma2, n)
	# Calculate CDF
	x2_mean = mean(X[1])
	x3_mean = mean(X[2])
	m2 = (sigma2*beta_post/n + tau2_post*x2_mean) / (sigma2/n + tau2_post)
	m3 = (sigma2*beta_post/n + tau2_post*x3_mean) / (sigma2/n + tau2_post)
	v2 = (sigma2*tau2_post/n) / (sigma2/n + tau2_post)
	rv = normal()
	for b in [0, 1, 3, 5]:
		Ps = rv.cdf( (m2 - m3 - b) / (sigma2*2 + v2*2) )
		print "For b = %d, P = %f" % (b, mean(Ps))

if __name__ == '__main__':
	main()