from pylab import *
from scipy.stats import beta
fact = math.factorial

def q3():
	loss = array([[0, 10], [10, 0]]) # action * theta
	X = array([[0.5, 0.5], [0.5, 0.5]]) # x * theta
	prior = array([0.4, 0.6])
	r = zeros((4,2))
	R = zeros(4)

	for i, t in enumerate([(0,0), (0,1), (1,0), (1,1)]):
		a1, a2 = t
		for theta in arange(0,2):
			risk = loss[a1,theta] * X[a1,theta] + loss[a2,theta] * X[a2,theta]
			r[i,theta] = risk
			R[i] += risk * prior[theta]

	print 'Frequentist risks:'
	print r
	print 'Min:', argmin(r,axis=0)

	print 'Bayesian risks:'
	print R
	print 'Min:', argmin(R,axis=0)

def q4a():
	X = range(21)
	a_prior, b_prior = 10, 10
	for x in X:
		a_posterior, b_posterior = x+10, 30-x
		beta_posterior = beta.cdf(0.5, a_posterior, b_posterior, loc=0, scale=1)
		print 'x = %d, cdf(0.5) = %f, diff = %f' % (x, beta_posterior, beta_posterior-2./3)

def q4b():
	C = range(21)
	f = zeros(21)
	for i,c in enumerate(C):
		f[i] = fact(20) / (fact(c) * fact(20-c)) * (0.5)**c * (1-0.5)**(20-c)
	plot(C,f,'o')
	show()
	print C
	print f

def q4c():
	x = 14
	diffs = zeros((20,20))
	for a_prior in arange(1,21):
		for b_prior in arange(1,21):
			a_posterior, b_posterior = x+a_prior, 20-x+b_prior
			beta_posterior = beta.cdf(0.5, a_posterior, b_posterior, loc=0, scale=1)
			diff = beta_posterior-2./3
			diffs[a_prior-1, b_prior-1] = diff
			# print 'x = %d, cdf(0.5) = %f, diff = %f' % (x, beta_posterior, diff)
	best_index = squeeze(argwhere(diffs == amin(abs(diffs))))
	print amin(abs(diffs))
	print 'Best a_prior = %d, b_prior = %d' % (best_index[0], best_index[1])

	# xx = linspace(0,1,1000)
	# beta_prior = beta.cdf(xx, a_prior, b_prior, loc=0, scale=1)
	# beta_posterior = beta.cdf(xx, a_posterior, b_posterior, loc=0, scale=1)
	# plot(xx, beta_prior)
	# title('prior, a = %d, b = %d'%(a_prior, b_prior))
	# figure()
	# plot(xx, beta_posterior)
	# title('posterior, a = %d, b = %d'%(a_posterior, b_posterior))
	# show()

def main():
	q4a()

if __name__ == '__main__':
	main()