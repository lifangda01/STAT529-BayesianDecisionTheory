from pylab import *

def get_optimal_rules_backward_induction(c):
	'''
		Obtain the optimal stopping rules for each stage.
	'''
	# Optimal thresholds from stage 0 to 14
	# Stage 15 doesn't have a rule
	thresh = zeros(15)
	E_beta = 25 - 15 * c
	thresh[-1] = 25 - c
	for n in range(13, -1, -1):
		y_n = arange(51.) - n * c
		E_beta = (sum(y_n[y_n > E_beta]) + E_beta * sum(y_n <= E_beta)) / 50.
		thresh[n] = E_beta + n * c
	return thresh

def main():
	set_printoptions(precision=3)
	c = 5
	E_beta = get_optimal_rules_backward_induction(c)
	print E_beta

if __name__ == '__main__':
	main()