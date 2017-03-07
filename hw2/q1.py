from pylab import *
from scipy.stats import beta as std_beta
from scipy.special import gamma

def get_gen_beta_pdf(x, alpha, beta, a, b):
	'''
		Return the pdf of a generalized beta function on given points.
	'''
	return gamma(alpha+beta)*1.0 / (gamma(alpha)*gamma(beta)) * \
			( (x-a)**(alpha-1) * (b-x)**(beta-1) *1.0 / (b-a)**(alpha+beta-1) )

def find_best_gen_beta_fit(hist_data):
	'''
		Find the (alpha, beta) of the closest generalized (-0.5 to 1.5) Beta distribution.
	'''
	hist_data = hist_data / sum(hist_data)
	hist_cdf = cumsum(hist_data)
	best_alpha, best_beta = 0, 0
	best_err = 99
	for alpha in linspace(0.1,10,100):
		for beta in linspace(0.1,10,100):
			cdf = cumsum( get_gen_beta_pdf(linspace(-0.5,1.5,20), alpha, beta, -0.5, 1.5) )
			# Restrict on [0,1]
			cdf = cdf[4:14] / (cdf[14] - cdf[4])
			max_err = max(abs(cdf - hist_cdf))
			if max_err < best_err:
				best_err = max_err
				best_alpha, best_beta = alpha, beta
	print "Best gen Beta fit is alpha = %f, beta = %f" % (best_alpha, best_beta) 
	plot(hist_cdf, label='data')
	cdf = cumsum( get_gen_beta_pdf(linspace(-0.5,1.5,20), best_alpha, best_beta, -0.5, 1.5) )
	cdf = cdf[4:14] / (cdf[14] - cdf[4])
	plot( cdf , label='fit')
	legend(loc=4)
	show()
	return best_alpha, best_beta

def find_bayes_interval(alpha, beta):
	'''
		Find the interval that gives the minimal expected posterior loss using gen Beta.
	'''
	cdf = cumsum( get_gen_beta_pdf(linspace(-0.5,1.5,200), alpha, beta, -0.5, 1.5) )
	best_a, best_b = 0, 0
	best_loss = 999
	for a in linspace(0.01,1.00,100):
		for b in arange(a,1.0,0.01):
			loss = 1 + (b-a-1) * (cdf[int((b+0.5) / 0.01)] - cdf[int((a+0.5) / 0.01)]) / (cdf[149] - cdf[49])
			if loss < best_loss:
				best_loss = loss
				best_a, best_b = a, b
	print "Best Bayes interval is (%f, %f)" % (best_a, best_b) 
	return best_a, best_b

def find_bayes_credible_region(alpha, beta, perc):
	'''
		Find the Bayes credible region by knocking off the head and tail using gen Beta.
	'''
	cdf = cumsum( get_gen_beta_pdf(linspace(-0.5,1.5,200), alpha, beta, -0.5, 1.5) )
	cdf = cdf[49:149] / (cdf[149] - cdf[49])
	for a in arange(len(cdf)):
		b = len(cdf)-a-1
		if cdf[b]-cdf[a] < perc:
			print "Bayes credible region is (%f, %f)" % (a*1.0/len(cdf), b*1.0/len(cdf)) 
			return a, b

def main():
	hist_int = array([0.45,0.74,1.07,1.22,1.47,1.50,1.48,1.09,0.61,0.36])
	hist_dom = array([0.55,0.90,1.54,1.59,1.52,1.30,1.04,0.82,0.47,0.25])
	alpha, beta = find_best_gen_beta_fit(hist_int) # Best gen Beta fit is alpha = 5.800000, beta = 6.700000
	find_bayes_interval(alpha, beta) # Best Bayes interval is (0.180000, 0.630000)
	find_bayes_credible_region(alpha, beta, 0.95) # Bayes credible region is (0.060000, 0.930000)

if __name__ == '__main__':
	main()