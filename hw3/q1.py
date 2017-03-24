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
	return theta_ij[500:]

def main():
	X = [array([12.7,6.6,14.7,12.2,4.4,7.8,13.8,13.7,11.1,9.1,14.0]),
		array([17.1,11.9,12.7,16.8,15.0,14.6,13.7,16.4]),
		array([5.2,4.5,10.5,15.0,5.0,14.9,7.6,8.3,10.8,14.6,15.1,7.0,9.3]),
		array([14.3,16.2,10.0,13.1,16.9,11.2,10.1,18.3,13.5,15.0,15.1,14.8,15.7,13.2,12.2,13.2]),
		array([10.5,7.5,4.7,12.5,13.1,13.5,12.2,16.1,9.0,17.9])]
	get_posterior_gibbs_sampler(X)

if __name__ == '__main__':
	main()