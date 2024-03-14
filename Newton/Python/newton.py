import math
import numpy as np

def f1(x):
    y=np.zeros(2)
    y[0]=2*x[0]-math.cos(x[1])
    y[1]=2*x[1]+math.sin(x[0])
    return y

def Jacobiana(x):
    J=np.zeros((2,2))
    J[0,0]=2
    J[0,1]=math.sin(x[1])
    J[1,0]=math.cos(x[0])
    J[1,1]=1
    return J

def newton(atol, rtol, max_it, x0, J, f):
    n_it = 0
    errore = True
    while errore:
        dx = np.linalg.solve(J(x0), f(x0))
        x1 = x0 - dx
        n_it += 1
        errore = np.linalg.norm(np.abs(x1 - x0) / (atol + rtol * np.abs(x1)), np.inf) > 1 and n_it < max_it
        x0 = x1
    return x1, n_it


x=np.array([0,0])
print(x)
x,n_it=newton(1e-5,1e-5,1000,x,Jacobiana,f1)
print(x)
print(n_it)