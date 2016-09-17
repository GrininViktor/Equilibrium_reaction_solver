import numpy as np

e = 0.00001

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

def test_a(N0,V,reactionDepth):
	if (np.min(N0 + np.dot(V,reactionDepth)) >= 0):
		return 1
	else:
		return 0
	

def function_F(K, V, N0, reactionDepth):
	return np.log(K) - np.dot(V.transpose(),np.log(N0 + np.dot(V,reactionDepth)))

def function_dF(V, N0, reactionDepth):
	return np.dot(np.dot(V.transpose(), np.linalg.inv(np.diag(N0 + np.dot(V,reactionDepth)))),V)

def newtons_method(reactionDepth1, reactionDepth2, e):
	while True:
		print reactionDepth2
		a = 1
		i = 0
		while (i == 0) :
			reactionDepth2 = reactionDepth1 + a*np.dot(np.linalg.inv(function_dF(V,N0,reactionDepth1)),function_F(K,V,N0,reactionDepth1))
			if (test_a(N0,V,reactionDepth2) == 1):
				i = 1	
			a = a/float(2)
		if np.linalg.norm( reactionDepth2 - reactionDepth1) < e:
			return reactionDepth2
		reactionDepth1 = reactionDepth2 
		
print(newtons_method(reactionDepth1, reactionDepth2, e))



#My_Repository
