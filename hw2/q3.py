from pylab import *
from scipy.interpolate import interp1d
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
	thetas = linspace(0.1,10,100)
	for alpha in linspace(0.1,10,100):
		for beta in linspace(0.1,10,100):
			P = zeros(7)
			for x in arange(1,8):
				P[x-1] = sum( get_possion_pmf(x, thetas) * get_gamma_pdf(x, alpha, beta) ) * 0.1
			max_err = max( abs(P - hist_data) )
			if max_err < best_err:
				best_err = max_err
				best_alpha, best_beta = alpha, beta
				best_P = copy(P)
	print "Best gamma fit is alpha = %f, beta = %f" % (best_alpha, best_beta) 
	plot(hist_data, label='data')
	plot(best_P, label='fit')
	legend(loc=4)
	show()
	return best_alpha, best_beta

def find_bayes_credible_region(alpha, beta, perc):
	'''
		Find the Bayes credible region by knocking off the head and tail using Gamma.
	'''
	cdf = cumsum( get_gamma_pdf(linspace(0.01,10.0,1000), alpha+beta, beta/(beta+1)) ) * 0.001
	for a in arange(len(cdf)):
		b = len(cdf)-a-1
		if cdf[b]-cdf[a] < perc:
			print "Bayes credible region is (%f, %f)" % (a*10.0/len(cdf), b*10.0/len(cdf)) 
			return a, b

def main():
	data = array([4,3,5,5,1,3,1,4,3,5,6,2,4,4,3,4,3,4,4,3,4,2,4,4,6,4,5,3,2,3,7,7,2,5,1,3,6,5,7,4,4,3,2,3,6,3,1,2,1,3])
	hist_data, _ = histogram(data, bins=7)
	hist_data = array(hist_data, dtype=float)
	alpha, beta = find_best_gamma_fit(hist_data)
	find_bayes_credible_region(alpha, beta, 0.95)

if __name__ == '__main__':
	main()
