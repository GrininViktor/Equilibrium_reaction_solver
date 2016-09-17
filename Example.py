import numpy as np

K = 3
a = 1
e = 0.1

reactionDepth1, reactionDepth2 = np.array([0,0,0]),np.array([0,0,0])

V = np.array([[1,0,0],
[-1,0,0],
[-1,-1,-2],
[0,1,0],
[0,-1,0],
[0,0,1],
[0,0,-1]])

K = np.array([1, 10**0.7, 10**3.4576])

N0 = np.array([1,0.1,0.1,1,0.1,1,0.1])

def function_F(K, V, N0, reactionDepth1):
	return np.log(K) - np.dot(V.transpose(),np.log(N0 + np.dot(V,reactionDepth1)))

def function_dF(V, N0, reactionDepth1):
	return np.dot(np.dot(V.transpose(), np.diag(N0 + np.dot(V,reactionDepth1))),V)
a = 1/float(2)
print V.T
