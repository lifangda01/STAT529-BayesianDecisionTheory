from pylab import *

def get_gamma_pdf(x, alpha, beta):
	'''
		Return the pdf of custom gamma function on given points.
	'''
	return 1.0 / (beta**alpha * gamma(alpha)) * x**(alpha+1) * exp(-x/beta)

def invgamma_rvs(alpha, beta, N=1):
	'''
		Return samples of the custom inverse Gamma RV.
	'''
	r = arange(0,20.01,0.01)
	gamma_pdf = get_gamma_pdf(r, 6, 0.05)
	gamma_pdf = gamma_pdf / sum(gamma_pdf)
	if N == 1:
		return 1. / choice(r, p=gamma_pdf)
	else:
		return 1. / choice(r, size=N, p=gamma_pdf)

def get_posterior_samples_gibbs(X):
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

def main():
	set_printoptions(precision=3)
	X = [array([12.7,6.6,14.7,12.2,4.4,7.8,13.8,13.7,11.1,9.1,14.0]),
		array([17.1,11.9,12.7,16.8,15.0,14.6,13.7,16.4]),
		array([5.2,4.5,10.5,15.0,5.0,14.9,7.6,8.3,10.8,14.6,15.1,7.0,9.3]),
		array([14.3,16.2,10.0,13.1,16.9,11.2,10.1,18.3,13.5,15.0,15.1,14.8,15.7,13.2,12.2,13.2]),
		array([10.5,7.5,4.7,12.5,13.1,13.5,12.2,16.1,9.0,17.9])]
	posteriors = get_posterior_samples_gibbs(X)
	# Q1a
	# Find index of max theta for every row
	max_indices = argmax(posteriors, axis=1)
	indices, counts = unique(max_indices, return_counts=True)
	epl, P = ones(5), zeros(5)
	epl[indices] = 1 - 1.0 * counts / sum(counts)
	P[indices] = 1.0 * counts / sum(counts)
	print "EPL:", epl
	print "Posterior P:", P
	print "Bayes selected diet:", argmin(epl)+1
	# Q1b
	P_b = zeros(4)
	for i, b in enumerate([0,2,4,6]):
		P_b[i] = sum(posteriors[:,1] >= posteriors[:,2] + b) * 1.0 / posteriors.shape[0]
	print "P_b:", P_b
	# Q1c
	P_2 = sum(posteriors[:,1] >= 14) / posteriors.shape[0]
	temp = posteriors[:,1][posteriors[:,2] <= 10] # P(theta2 | theta3 <= 10)
	P_23 = sum(temp >= 14) * 1.0 / temp.shape[0]
	print "P(theta2 >= 14 | theta3 <= 10):", P_23
	print "P(theta2 >= 14):", P_2

if __name__ == '__main__':
	main()