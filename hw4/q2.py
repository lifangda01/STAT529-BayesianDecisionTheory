from pylab import *

def get_optimal_rules_backward_induction(c):
	'''
		Obtain the optimal stopping thresholds for each stage.
		Thresholds are designated for x.
	'''
	# Optimal thresholds from stage 0 to 14
	# Stage 15 doesn't have a rule
	thresh = zeros(15)
	thresh[14] = 25 - c # Stage 14
	E_beta = zeros(15)
	E_beta[14] = 25 - 15 * c # beta_14^15
	for n in range(13, 0, -1):
		thresh[n] = E_beta[n + 1] + n * c
		E_beta[n] = (E_beta[n + 1] + n * c)**2 / 100. + 25 - n * c
	return thresh[1:], E_beta[1] - c

def get_rules_osa(c):
	'''
		Get OSA x thresholds, but no optimality guaranteed.
	'''
	return ones(15) * (25 - c)

def get_expected_payoff_osa(c, thresh):
	'''
		Repeat the game and return the average payoff.
	'''
	N = 1000
	Y = zeros(N)
	for i in range(N):
		j = 0
		x = rand() * 50
		while j < 15 and x < thresh[j]:
			j += 1
			x = rand() * 50
		Y[i] = x - j * c
	return mean(Y)

def main():
	set_printoptions(precision=3, suppress=True)
	C = [0, 1, 2, 5, 10, 25, 50]
	Y_opt = zeros(len(C))
	Y_osa = zeros((len(C)))
	figure()
	for i, c in enumerate(C):
		t_opt, Y_opt[i] = get_optimal_rules_backward_induction(c)
		t_osa = get_rules_osa(c)
		Y_osa[i] = get_expected_payoff_osa(c, t_osa)
		print "c =", c
		print t_opt
		print "Optimal expected payoff =", Y_opt[i]
		print "OSA expected payoff =", Y_osa[i]
		step(arange(14) + 1, t_opt, label="c = %d" % c)
	title('Optimal Stopping Thresholds for x_i')
	xlabel('Stage')
	ylabel('Threshold')
	legend()
	# Plot the expected payoffs
	figure()
	bar_width = 0.35
	index = arange(len(C))
	bar(index, Y_opt, bar_width, label='OPT')
	bar(index + bar_width, Y_osa, bar_width, label='OSA')
	title('Expected Payoffs: Optimal v.s. OSA')
	xlabel('c')
	ylabel('Payoff')
	xticks(index, [str(c) for c in C])
	legend()
	tight_layout()
	show()

if __name__ == '__main__':
	main()