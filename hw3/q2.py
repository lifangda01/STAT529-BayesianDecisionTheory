from pylab import *

def get_corn_posterior_samples_MCMC(X):
	'''
		Using a HB model, A/R algorithm with MCMC to generate posterior samples.
	'''
	beta0 = 125.
	A = 15.
	sigma2 = 80.
	N = 2000
	# Some preprocessing
	Y_mean = mean(X.reshape(4,-1), axis=1)
	# Posteriors before selection, each column represents one species
	lambda_ij = zeros((N, 4))
	# Initial guesses
	beta = randn()*A + beta0
	tau2 = rand()*14 + 2
	lambda_ij[0] = randn(4)*sqrt(tau2/3) + beta
	# Generating loop
	j = 1
	while j < N:
		beta = randn()*A + beta0
		tau2 = rand()*14 + 2
		lambda_i = randn(4)*sqrt(tau2/3) + beta
		# Determine whether to keep
		r = min(exp( -0.5 * sum((lambda_i - Y_mean)**2 - (lambda_ij[j-1] - Y_mean)**2) / (sigma2 / 12) ),1)
		# Keep with probability r
		if rand() > r: continue
		lambda_ij[j] = lambda_i
		j += 1
	# Discard the first 500 to minimize the influence of the intial guess
	lambda_ij = lambda_ij[500:]
	# hist(lambda_ij, bins=10, label=["Corn %d"%(i) for i in range(1,5)])
	# legend()
	# show()
	return lambda_ij

def get_fertilizer_posterior_samples_MCMC(X):
	'''
		Using a HB model, A/R algorithm with MCMC to generate posterior samples.
	'''
	beta0 = 125.
	A = 15.
	sigma2 = 80.
	N = 2000
	# Some preprocessing
	T_mean = mean(X.T, axis=1)
	# Posteriors before selection, each column represents one fertilizer
	w_ij = zeros((N, 3))
	# Initial guesses
	beta = randn()*A + beta0
	tau2 = rand()*14 + 2
	w_ij[0] = randn()*sqrt(tau2/4) + beta
	# Generating loop
	j = 1
	while j < N:
		beta = randn()*A + beta0
		tau2 = rand()*14 + 2
		w_i = randn(3)*sqrt(tau2/4) + beta
		# Determine whether to keep
		r = min(exp( -0.5 * sum((w_i - T_mean)**2 - (w_ij[j-1] - T_mean)**2) / (sigma2 / 16) ),1)
		# Keep with probability r
		if rand() > r: continue
		w_ij[j] = w_i
		j += 1
	# Discard the first 500 to minimize the influence of the intial guess
	w_ij = w_ij[500:]
	# hist(w_ij, bins=10, label=["Fertilizer %d"%(i) for i in range(1,4)])
	# legend()
	# show()
	return w_ij

def get_combined_posterior_samples_MCMC(X):
	'''
		Using a HB model, A/R algorithm with MCMC to generate posterior samples.
	'''
	beta0 = 125.
	A = 15.
	sigma2 = 80.
	N = 2000
	# Some preprocessing
	X_mean = mean(X.T.flatten().reshape(-1,4),axis=1)
	X_mean = X_mean[array([1,4,7,10,2,5,8,11,3,6,9,12])-1]
	X_mean = X_mean[[3,11]]
	# Posteriors before selection, each column represents one fertilizer
	theta_ij = zeros((N, 2))
	# Initial guesses
	beta = randn()*A + beta0
	tau2 = rand()*14 + 2
	theta_ij[0] = randn()*sqrt(tau2) + beta
	# Generating loop
	j = 1
	while j < N:
		beta = randn()*A + beta0
		tau2 = rand()*14 + 2
		theta_i = randn(2)*sqrt(tau2) + beta
		# Determine whether to keep
		r = min(exp( -0.5 * sum((theta_i - X_mean)**2 - (theta_ij[j-1] - X_mean)**2) / (sigma2 / 16) ),1)
		# Keep with probability r
		if rand() > r: continue
		theta_ij[j] = theta_i
		j += 1
	# Discard the first 500 to minimize the influence of the intial guess
	theta_ij = theta_ij[500:]
	hist(theta_ij, bins=10, label=["theta %d"%(i) for i in [4,12] ])
	legend()
	show()
	return theta_ij

def main():
	set_printoptions(precision=3)
	X = array(([138,122,121], [138,127,119], [129,124,118], [131,123,122], 
				[135,130,133], [140,140,130], [136,140,132], [130,128,128], 
				[125,120,115], [140,115,112], [140,110,110], [125,120,115],
				[130,118,141], [130,115,140], [140,112,139], [125,116,142]))
	# Q2
	# Corn-wise
	lambda_ij = get_corn_posterior_samples_MCMC(X)
	P2c = 1.0*sum(argmax(lambda_ij, axis=1) == 1) / lambda_ij.shape[0]
	print"P(lambda2 | X):", P2c
	# Fertilizer-wise
	w_ij = get_fertilizer_posterior_samples_MCMC(X)
	P1f = 1.0*sum(argmax(w_ij, axis=1) == 0) / w_ij.shape[0]
	print"P(W1 | X):", P1f
	# theta4 vs theta12
	theta_ij = get_combined_posterior_samples_MCMC(X)
	P124 = 1.0*sum(argmax(theta_ij, axis=1) == 1) / theta_ij.shape[0]
	print"P(theta12 >= theta4 | X):", P124

if __name__ == '__main__':
	main()