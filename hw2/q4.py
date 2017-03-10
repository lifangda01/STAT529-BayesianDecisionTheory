from pylab import *
from scipy.stats import beta as std_beta

def find_probability_posterior(alpha_int, beta_int, alpha_dom, beta_dom):
	'''
		Compute the probability that theta_dom >= theta_int using posterior.
	'''
	N = 10000
	thetah_int = 2.0*16 / 33 -0.5
	thetah_dom = 2.0*16 / 27 -0.5
	thetas_int = std_beta.rvs(alpha_int, beta_int, size=N)
	thetas_dom = std_beta.rvs(alpha_dom, beta_dom, size=N)
	denom = (thetah_dom+0.5)**16 * (1.5-thetah_dom)**11 * (thetah_int+0.5)**16 * (1.5-thetah_int)**17
	P_keep = (thetas_dom+0.5)**16 * (1.5-thetas_dom)**11 * (thetas_int+0.5)**16 * (1.5-thetas_int)**17 / denom
	indices = rand(N) <= P_keep
	keepers_int = thetas_int[indices]
	keepers_dom = thetas_dom[indices]
	bools = keepers_dom >= keepers_int
	P_result = sum(bools)*1.0 / bools.size
	print "Probability of theta_dom >= theta_int is %f" % P_result
	return P_result

def find_probability_hb():
	'''
		Compute the probability that theta_dom >= theta_int using hyper priors.
	'''
	N = 1000	
	thetah_int = 2.0*16 / 33 -0.5
	thetah_dom = 2.0*16 / 27 -0.5
	num_keepers = 0
	num_total = 0
	for i in range(1000):
		# Generate prior distribution
		alpha = rand()*4.2 + 1.8 # uniform on [1.8,6]
		beta = rand()*(0.75*alpha+3 - 0.3*alpha+1.75) + 0.3*alpha+1.75 # uniform on [0.3*alpha+1.75, 0.75*alpha+3]
		thetas_int = std_beta.rvs(alpha, beta, size=N)
		thetas_dom = std_beta.rvs(alpha, beta, size=N)
		# Same as using posterior
		denom = (thetah_dom+0.5)**16 * (1.5-thetah_dom)**11 * (thetah_int+0.5)**16 * (1.5-thetah_int)**17
		P_keep = (thetas_dom+0.5)**16 * (1.5-thetas_dom)**11 * (thetas_int+0.5)**16 * (1.5-thetas_int)**17 / denom
		indices = rand(N) <= P_keep
		keepers_int = thetas_int[indices]
		keepers_dom = thetas_dom[indices]
		bools = keepers_dom >= keepers_int
		num_keepers += sum(bools)
		num_total += bools.size
	P_result = num_keepers*1.0 / num_total
	print "Probability of theta_dom >= theta_int is %f" % P_result
	return P_result

def main():
	alpha_int, beta_int = 1.2, 1.5 # Best std Beta fit is alpha = 1.200000, beta = 1.500000
	alpha_dom, beta_dom = 0.9, 1.4 # Best std Beta fit is alpha = 0.900000, beta = 1.400000
	find_probability_posterior(alpha_int, beta_int, alpha_dom, beta_dom) # Probability of theta_dom >= theta_int is 0.754662
	find_probability_hb() # Probability of theta_dom >= theta_int is 0.696833


if __name__ == '__main__':
	main()