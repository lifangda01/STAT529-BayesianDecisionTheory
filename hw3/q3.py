from pylab import *

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

def get_ac_posterior_samples_MCMC(X):
	'''
		Using a HB model, Metropolis Hasting algorithm with MCMC to generate posterior samples.
	'''
	N = 2000
	# Some preprocessing, we are only interested in AC 1 and 2
	X = X[:,:2]
	X_mean = X[1]*1.0 / X[0]
	X_num = X[0]
	# Posteriors before selection, each column represents one AC
	theta_ij = zeros((N, X.shape[1]))
	# Initial guesses
	alpha = rand()*35 + 30
	beta = rand()*1.5 - 0.057*alpha + 5
	theta_ij[0] = gamma_rvs(alpha, beta, N=X.shape[1])
	# Generating loop
	j = 1
	while j < N:
		alpha = rand()*35 + 30
		beta = rand()*1.5 - 0.057*alpha + 5
		theta_i = gamma_rvs(alpha, beta, N=X.shape[1])
		# Determine whether to keep
		r = prod( (theta_ij[j-1]/theta_i)**X_num * exp(X_num*X_mean*(1/theta_ij[j-1] - 1/theta_i)) )
		r = min(r, 1)
		# Keep with probability r
		if rand() > r: continue
		theta_ij[j] = theta_i
		j += 1
	# Discard the first 500 to minimize the influence of the intial guess
	theta_ij = theta_ij[500:]
	hist(theta_ij, bins=10, label=["AC %d"%(i) for i in range(1,3)])
	legend()
	show()
	return theta_ij

def main():
	X = array([[14,12,23,16,27,30],
		[1832,1297,2201,1312,2074,1788]])
	# We only compute posteriors for AC 1 and 2
	posteriors = get_ac_posterior_samples_MCMC(X)
	P_b = zeros(2)
	for i, b in enumerate([10,20]):
		P_b[i] = sum(posteriors[:,0] >= posteriors[:,1] + b) * 1.0 / posteriors.shape[0]
	print "P_b:", P_b

if __name__ == '__main__':
	main()