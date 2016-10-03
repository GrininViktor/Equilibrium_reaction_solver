# vim: noet
import numpy as np
np.seterr(all='raise')

e = 1e-20
# 0.00000000000000000001

# reactionDepth1, reactionDepth2, reactionDepth0 = np.array([-1.9975799,-5.8705*(10**(-8)),0.941224, 0.99758]), np.array([-1.9975799,-5.8705*(10**(-8)),0.941224, 0.99758]), np.array([-1.9975799,-5.8705*(10**(-8)),0.941224, 0.99758])
reactionDepth1 = np.array([-1.9975799,-5.8705e-8, 0.941224, 0.99758])

#reactionDepth2 = reactionDepth1.copy()
#reactionDepth0 = reactionDepth1.copy()

#[  1.99757990e+00   9.88482786e-08   5.87759607e-02   2.41999997e-03
#   0.00000000e+00   2.93880394e+00   2.41999997e-03   9.97580000e-01]
#reactionDepth1, reactionDepth2 = np.array([-1,-0.5,0.5,0.7]),np.array([-1,-0.5,0.5,0.7])

species = ['OH-', 'HCO3+', 'CO2g', 'CaCO3', 'H2O', 'CO2l', 'H+', 'Ca2+']

V = np.array([[-1,0,0,0],
[0,-1,0,0],
[0,0,-1,0],
[0,0,0,-1],
[1,1,0,1],
[0,1,1,1],
[-1,-1,0,-2],
[0,0,0,1]])

for j in range(V.shape[1]):
	s = '0 = '
	for i in range(V.shape[0]):
		if V[i, j] != 0:
			s += str(V[i,j]) + ' ' + species[i] + ' + '
	print s

K = np.array([10**(-14), 10**(-5.928), 50, 10**(-8.094)])

N0 = np.array([0,0,1,1,1,1,0,0])

def test_a(N0, V, reactionDepth, C):
	a = 1
	N = N0 + V.dot(reactionDepth)
	dN = V.dot(C)
	#print '\n',N
	#print '\n', dN
	for i in range(len(N)):
		if (N[i] + dN[i] < 0):
			if (-N[i]/dN[i] < a):
				a = -N[i]/dN[i]
	return a

def function_F(K, V, N0, reactionDepth):
	N = N0 + V.dot(reactionDepth)
	print 'N =', list(zip(species, N))
	return np.log(K) - V.T.dot(np.log(N))

def function_dF(V, N0, reactionDepth):
	#return V.T.dot(np.linalg.solve(np.diag(N0 + V.dot(reactionDepth)),V))
	N = N0 + V.dot(reactionDepth)
	iN = 1. / N
	return V.T.dot(np.diag(iN)).dot(V)

#print function_F(K, V, N0, reactionDepth1)
def newtons_method(reactionDepth, e):
	l = 0
	while True:
		print 
		l += 1
		reactionDepthOld = reactionDepth.copy()

		A = function_F(K,V,N0,reactionDepth)
		print 'residual =', A
#		B = np.linalg.inv(function_dF(V,N0,reactionDepth))
#		C = B.dot(A)

		B = function_dF(V, N0, reactionDepth)
		C = np.linalg.solve(B + 100 * np.eye(4), A)
		print 'dx =', C

		a = test_a(N0, V, reactionDepth, C)
		print 'it =', l, '||A|| =', np.linalg.norm(A), 'a =', a

		reactionDepth += 0.5 * a * C

		print 'xi =', reactionDepth

		if np.linalg.norm(reactionDepth - reactionDepthOld) < e:
			return reactionDepth, l

print newtons_method(reactionDepth1, e)
