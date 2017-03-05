from pylab import *
from scipy.stats import beta

# Find the (alpha, beta) of the closest Beta distribution
def find_best_std_beta_fit(hist_data):
	hist_data = hist_data / sum(hist_data)
	accum_int = cumsum(hist_data)
	best_a, best_b = 0, 0
	best_err = 99
	for a in linspace(0.1,10,100):
		for b in linspace(0.1,10,100):
			max_err = max(abs(beta.cdf(linspace(0.0,1.0,10), a, b) - accum_int))
			if max_err < best_err:
				best_err = max_err
				best_a, best_b = a, b
	print "Best fit is a = %f, b = %f" % (best_a, best_b) 
	plot(accum_int)
	plot( beta.cdf(linspace(0.0,1.0,10), best_a, best_b) )
	show()

def find_best_gen_beta_fit(hist_data):
	hist_data = hist_data / sum(hist_data)
	accum_int = cumsum(hist_data)
	best_a, best_b = 0, 0
	best_err = 99
	for a in linspace(0.1,10,100):
		for b in linspace(0.1,10,100):
			max_err = max(abs(beta.cdf(linspace(0.0,1.0,10), a, b) - accum_int))
			if max_err < best_err:
				best_err = max_err
				best_a, best_b = a, b
	print "Best fit is a = %f, b = %f" % (best_a, best_b) 
	plot(accum_int)
	plot( beta.cdf(linspace(0.0,1.0,10), best_a, best_b) )
	show()

def main():
	hist_int = array([0.45,0.74,1.07,1.22,1.47,1.50,1.48,1.09,0.61,0.36])
	hist_dom = array([0.55,0.90,1.54,1.59,1.52,1.30,1.04,0.82,0.47,0.25])
	find_best_beta_fit(hist_int) # Best fit is a = 1.200000, b = 1.500000
	find_best_beta_fit(hist_dom) # Best fit is a = 0.900000, b = 1.400000

if __name__ == '__main__':
	main()