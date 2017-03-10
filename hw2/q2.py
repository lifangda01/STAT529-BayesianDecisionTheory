from pylab import *
from scipy.stats import beta as std_beta
from scipy.interpolate import interp1d

def find_best_std_beta_fit(hist_data):
	'''
		Find the (alpha, beta) of the closest standard Beta distribution.
	'''
	hist_data = hist_data / sum(hist_data)
	hist_cdf = cumsum(hist_data)
	best_alpha, best_beta = 0, 0
	best_err = 99
	for alpha in linspace(0.1,10,100):
		for beta in linspace(0.1,10,100):
			max_err = max(abs(std_beta.cdf(linspace(0.0,1.0,10), alpha, beta) - hist_cdf))
			if max_err < best_err:
				best_err = max_err
				best_alpha, best_beta = alpha, beta
	print "Best std Beta fit is alpha = %f, beta = %f" % (best_alpha, best_beta) 
	plot(hist_cdf, label='data')
	plot( std_beta.cdf(linspace(0.0,1.0,10), best_alpha, best_beta) , label='fit')
	legend(loc=4)
	show()
	return best_alpha, best_beta

def find_bayes_credible_region_std_beta(alpha, beta, perc, y=16, n=33, thetah=0.47):
	'''
		Find the Bayes credible region by knocking off the head and tail using std Beta.
	'''
	N = 10000
	thetas = std_beta.rvs(alpha, beta, size=N)
	denom = (thetah+0.5)**y * (1.5-thetah)**(n-y)
	P = (thetas+0.5)**y * (1.5-thetas)**(n-y) / denom
	keepers = thetas[rand(N) <= P]
	keepers = sort(keepers)
	offset = int((1 - perc) * len(keepers))
	print "Bayes credible region is (%f, %f)" % (keepers[offset], keepers[-offset]) 
	return keepers[offset], keepers[-offset]

def find_bayes_credible_region_jeffery(perc, y=16, n=33, thetah=0.47):
	'''
		Find the Bayes credible region by knocking off the head and tail using Jeffery's prior.
	'''
	return find_bayes_credible_region_std_beta(0.5, 0.5, perc, y, n, thetah)

def find_bayes_credible_region_hist(hist_data, perc, y=16, n=33, thetah=0.47):
	'''
		Find the Bayes credible region by knocking off the head and tail using histogram data prior.
	'''
	hist_data = hist_data / sum(hist_data)
	hist_cdf = cumsum(hist_data)
	hist_cdf = hstack((zeros(1),hist_cdf))
	# Use interpolation to generate continuous cdf
	F = interp1d(linspace(0.0, 1.0, 11), hist_cdf)
	cdf = F(linspace(0.0,1.0,1000))
	for a in arange(len(cdf)):
		b = len(cdf)-a-1
		if cdf[b]-cdf[a] < perc:
			print "Bayes credible region is (%f, %f)" % (a*1.0/len(cdf), b*1.0/len(cdf)) 
			return a, b

def main():
	hist_int = array([0.45,0.74,1.07,1.22,1.47,1.50,1.48,1.09,0.61,0.36])
	hist_dom = array([0.55,0.90,1.54,1.59,1.52,1.30,1.04,0.82,0.47,0.25])
	# alpha, beta = find_best_std_beta_fit(hist_int) # Best fit is a = 1.200000, b = 1.500000
	find_best_std_beta_fit(hist_int) # Best std Beta fit is alpha = 1.200000, beta = 1.500000
	find_best_std_beta_fit(hist_dom) # Best std Beta fit is alpha = 0.900000, beta = 1.400000
	# find_bayes_credible_region_std_beta(alpha, beta, 0.95) # Bayes credible region is (0.202, 0.717)
	# find_bayes_credible_region_jeffery(0.95) # Bayes credible region is (0.165, 0.776)
	# find_bayes_credible_region_hist(hist_int, 0.95) # Bayes credible region is (0.062000, 0.937000)


if __name__ == '__main__':
	main()