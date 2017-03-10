from pylab import *
from scipy.special import gamma
from scipy.special import factorial

def get_gamma_pdf(x, alpha, beta):
	'''
		Return the pdf of gamma function on given points.
	'''
	return 1.0 / (beta**alpha * gamma(alpha)) * x**(alpha-1) * exp(-x/beta)

def get_possion_pmf(x, theta):
	'''
		Return the pdf of possion function on given points.
	'''
	return exp(-theta) * theta**x / factorial(x)

def find_best_gamma_fit(hist_data):
	'''
		Find the best alpha and beta for the gamma prior.
	'''
	hist_data = hist_data / sum(hist_data)
	hist_cdf = cumsum(hist_data)
	best_alpha, best_beta = 0, 0
	best_err = 99
	X = arange(1,8)
	for alpha in arange(0.1,10,0.01):
		for beta in arange(0.1,10,0.01):
			P = 1.0 / factorial(X) * factorial(alpha+X-1) / factorial(alpha-1) * beta**X * (1+beta)**(alpha+X)
			max_err = max( abs(P[3:5] - hist_data[3:5]) )
			# max_err = max( abs(P - hist_data) )
			if max_err < best_err:
				best_err = max_err
				best_alpha, best_beta = alpha, beta
				best_P = copy(P)
	print "Best gamma fit is alpha = %f, beta = %f" % (best_alpha, best_beta) 
	plot(range(1,8), hist_data, label='data')
	plot(range(1,8), best_P, label='fit')
	legend(loc=1)
	show()
	return best_alpha, best_beta

def find_bayes_credible_region_gamma(alpha, beta, perc):
	'''
		Find the Bayes credible region by knocking off the head and tail using Gamma.
	'''
	r = arange(1,7.01,0.01)
	pdf = get_gamma_pdf(r, alpha, beta)
	pdf = pdf / sum(pdf)
	cdf = cumsum( pdf )
	for a in arange(len(cdf)):
		b = len(cdf)-a-1
		if cdf[b]-cdf[a] < perc:
			print "Bayes credible region is (%f, %f)" % (r[a], r[b]) 
			return a, b

def find_bayes_credible_region_hb(perc):
	'''
		Find the Bayes credible region using A/R and hyper priors.
	'''
	N = 100	
	keepers = array([])
	for i in range(1000):
		# Generate prior distribution
		alpha = rand()*5 + 1 # uniform on [1,6]
		beta = rand()*(6/alpha - 2/alpha) + 2/alpha # uniform on [2/alpha, 6/alpha]
		r = arange(1,7.01,0.01)
		pdf = get_gamma_pdf(r, alpha, beta)
		pdf = pdf / sum(pdf)
		# Generate thetas
		thetas = choice(r, size=N, p=pdf)
		M = max(get_possion_pmf(4, thetas))
		P = get_possion_pmf(4, thetas) / M
		keepers = append(keepers, thetas[rand(N) <= P])
	keepers = sort(keepers)
	# Knocking off the head and tail 2.5%
	offset = int((1 - perc) * len(keepers))
	print "Bayes credible region is (%f, %f)" % (keepers[offset], keepers[-offset]) 
	return keepers[offset], keepers[-offset]

def main():
	data = array([4,3,5,5,1,3,1,4,3,5,6,2,4,4,3,4,3,4,4,3,4,2,4,4,6,4,5,3,2,3,7,7,2,5,1,3,6,5,7,4,4,3,2,3,6,3,1,2,1,3])
	hist_data, _ = histogram(data, bins=7)
	hist_data = array(hist_data, dtype=float)
	alpha, beta = find_best_gamma_fit(hist_data) # Best gamma fit is alpha = 2.200000, beta = 0.300000
	find_bayes_credible_region_gamma(alpha+4, beta/(1+beta), 0.95) # Bayes credible region is (1.050000, 6.950000)
	find_bayes_credible_region_hb(0.95) # Bayes credible region is (1.810000, 6.180000)

if __name__ == '__main__':
	main()
