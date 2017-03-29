from pylab import *

def get_posterior_gibbs_sampler(X):
	'''
		Use Gibbs sample to generate posterior samples for A3Q1.
		X is a list of arrays for each diet.
	'''
	# Some preprocessing
	X_num = array([1.0*len(X[i]) for i in range(5)])
	X_mean = array([mean(X[i]) for i in range(5)])
	# Initial guesses
	beta = 2*randn() + 12
	tau2 = invgamma_rvs(6, 0.05)
	sigma2_i = invgamma_rvs(13.1, 0.0083)
	# Posteriors, each column represents one diet
	theta_ij = zeros((2000,5))
	theta_ij[0] = sqrt(tau2)*randn(5) + beta
	# Generating loop
	for j in range(1,2000):
		theta_i = theta_ij[j-1]
		beta = sqrt(tau2*4/(tau2+20))*randn() + (tau2*12+5*4*mean(theta_i)) / (tau2+5*4)
		tau2 = invgamma_rvs(2.5 + 6, 2*0.05 / (2 + 0.05*sum((theta_i-beta)**2)))
		sigma2_i = array([invgamma_rvs(13.1+X_num[i]/2., \
				2 / ( 2/0.0083 + sum( (X[i] - X_mean[i])**2 ) + X_num[i]*( X_mean[i] - theta_i[i] ) )) \
				for i in range(5)])
		theta_i = array([sqrt( sigma2_i[i]*tau2 / (sigma2_i[i] + X_num[i]*tau2) ) * randn() \
				+ (sigma2_i[i]*beta + tau2*X_num[i]*X_mean[i]) / (sigma2_i[i] + X_num[i]*tau2) \
				for i in range(5)])
		theta_ij[j] = theta_i
	# Discard the first 500 to minimize the influence of the intial guess
	theta_ij = theta_ij[500:]
	# hist(theta_ij, bins=50, label=["diet %d"%(i) for i in range(1,6)])
	# legend()
	# show()
	return theta_ij

def get_posterior_samples_MCMC(X):
	'''
		Using a HB model, A/R algorithm with MCMC to generate posterior samples.
	'''
	# Some initial guesses based on data
	beta = 127.
	tau2 = 94.29


def main():
	set_printoptions(precision=3)
	X = array(([138,122,121], [138,127,119], [129,124,118], [131,123,122], 
			[135,130,133], [140,140,130], [136,140,132], [130,128,128], 
			[125,120,115], [140,115,112], [140,110,110], [125,120,115],
			[130,118,141], [130,115,140], [140,112,139], [125,116,142]))


if __name__ == '__main__':
	main()