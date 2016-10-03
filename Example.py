import numpy as np

e = 0.0001

reactionDepth1, reactionDepth2 = np.array([0,0,0]),np.array([0,0,0])

V = np.array([[1,0,0],
[0,0,0],
[0,0,0],
[0,1,0],
[0,0,0],
[0,0,1],
[0,0,0]])

V1 = np.array([[1,0,0],
[-1,0,0],
[-1,-1,-2],
[0,1,0],
[0,-1,0],
[0,0,1],
[0,0,-1]])

K = np.array([1, 10**0.7, 10**3.4576])

N0 = np.array([1,0,0,1,0,1,0])

def test_a(N0, V, reactionDepth, C):
	a = 1
	N = N0 + V.dot(reactionDepth)
	dN = V.dot(C)
	for i in range(len(N)):
		if (N[i] + dN[i] < 0):
			if (abs(N[i]/dN[i]) < a):
				a = abs(N[i]/dN[i]) 
	return 0.99*a

def function_F(K, V, N0, reactionDepth):
	N = N0 + V.dot(reactionDepth)
	print N
	for i in range(len(N)):
		if (N[i] == 0):
			N[i] = 1
	return np.log(K) - V.T.dot(np.log(N))

def function_dF(V, N0, reactionDepth):
	#return V.T.dot(np.linalg.solve(np.diag(N0 + V.dot(reactionDepth)),V))
	return np.dot(np.dot(V.transpose(), np.diag(N0 + np.dot(V,reactionDepth1))),V)

def newtons_method(reactionDepth1, reactionDepth2, e):
	while True:
		A = function_F(K,V,N0,reactionDepth1)
		B = np.linalg.inv(function_dF(V,N0,reactionDepth1))
		C = B.dot(A)
		a = test_a(N0, V, reactionDepth1, C)
		reactionDepth2 = reactionDepth1 + a*C
		if np.linalg.norm( reactionDepth2 - reactionDepth1) < e:
			return reactionDepth2
		reactionDepth1 = reactionDepth2 
				
print(newtons_method(reactionDepth1, reactionDepth2, e))

#My_Repository
