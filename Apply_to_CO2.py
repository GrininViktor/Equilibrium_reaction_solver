import numpy as np

e = 0.00000000000000000001

reactionDepth1, reactionDepth2,reactionDepth0 = np.array([-1.9975799,-5.8705*(10**(-8)),0.941224, 0.99758]),np.array([-1.9975799,-5.8705*(10**(-8)),0.941224, 0.99758]),np.array([-1.9975799,-5.8705*(10**(-8)),0.941224, 0.99758])
#[  1.99757990e+00   9.88482786e-08   5.87759607e-02   2.41999997e-03
#   0.00000000e+00   2.93880394e+00   2.41999997e-03   9.97580000e-01]
#reactionDepth1, reactionDepth2 = np.array([-1,-0.5,0.5,0.7]),np.array([-1,-0.5,0.5,0.7])
V = np.array([[-1,0,0,0],
[0,-1,0,0],
[0,0,-1,0],
[0,0,0,-1],
[1,1,0,1],
[0,1,1,1],
[-1,-1,0,-2],
[0,0,0,1]])

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
	return np.log(K) - V.T.dot(np.log(N))
def function_dF(V, N0, reactionDepth):
	#return V.T.dot(np.linalg.solve(np.diag(N0 + V.dot(reactionDepth)),V))
	return np.dot(np.dot(V.transpose(), np.diag(N0 + np.dot(V,reactionDepth1))),V)
#print function_F(K, V, N0, reactionDepth1)
def newtons_method(reactionDepth1, reactionDepth2, e):
	while True:
		A = function_F(K,V,N0,reactionDepth0)
		#print '\n',np.linalg.norm(A), A
		B = np.linalg.inv(function_dF(V,N0,reactionDepth0))
		C = B.dot(A)
		a = test_a(N0, V, reactionDepth1, C)
		#print a
		reactionDepth2 = reactionDepth1 + 0.5*a*C
		print reactionDepth1, reactionDepth2
		reactionDepth1[0] = reactionDepth1[0] - 1;
		print reactionDepth1, reactionDepth2
		#if np.linalg.norm( reactionDepth2 - reactionDepth1) < e:
			#return reactionDepth2, l
		reactionDepth1 = reactionDepth2 	
print newtons_method(reactionDepth1, reactionDepth2, e)
