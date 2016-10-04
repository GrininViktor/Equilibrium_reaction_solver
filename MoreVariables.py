import numpy as np
np.seterr(all='raise')
e = 1e-20

variablesSet = np.array([-0.5,-0.7,0.5,0.5,-0.6931471805599453,-0.3566749439387324,-0.6931471805599453,-0.6931471805599453,-1.2039728043259359,0.2623642644674911,-1.6094379124341005,-0.6931471805599453])

VOld = np.array([[-1.0,0.0,0.0,0.0],
              [0.0,-1.0,0.0,0.0],
              [0.0,0.0,-1.0,0.0],
              [0.0,0.0,0.0,-1.0],
              [1.0,1.0,0.0,1.0],
              [0.0,1.0,1.0,1.0],
              [-1.0,-1.0,0.0,-2.0],
              [0.0,0.0,0.0,1.0]])

V = np.vstack((np.hstack((np.zeros((4,4)), -VOld.T)), np.hstack((-VOld, np.zeros((8,8))))))

K = np.array([10**(-14), 10**(-5.928), 50, 10**(-8.094)])

N0 = np.array([0.0,0.0,1.0,1.0,1.0,1.0,0.0,0.0])



def function_F(K, V, N0, variablesSet):
    return np.hstack((-np.log(K), -N0)) + np.diag([0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]).dot(np.exp(variablesSet)) + V.dot(variablesSet)

def function_dF(V, variablesSet):
    dF = V.copy()
    for i in range(8):
        dF[i+4][i+4] = np.exp(variablesSet[4+i])
    return dF 
print (function_dF(V, variablesSet))
def newtons_method(variablesSet, e):
    l = 0
    a = 0.5
    while True:
        l += 1
        variablesSetOld = variablesSet.copy()
        A = function_F(K, V, N0, variablesSet)
        print ('residual =', A)
        B = function_dF(V, variablesSet)
        C = np.linalg.solve(B, A)
        print ('dx =', C)
        print ('it =', l, '||A|| =', np.linalg.norm(A), 'a =', a)
        variablesSet -= a * C
        print ('xi =', variablesSet)
        if np.linalg.norm(variablesSet - variablesSetOld) < e:
            return variablesSet, l
 
print (newtons_method(variablesSet, e))
