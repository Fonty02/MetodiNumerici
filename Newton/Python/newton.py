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
    J[1,1]=2
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

def quasiNewton(atol, rtol, max_it, x0, J, f):
    n_it = 0
    errore = True
    A=J(x0)
    while errore:
        dx = np.linalg.solve(A, f(x0))
        x1 = x0 - dx
        n_it += 1
        errore = np.linalg.norm(np.abs(x1 - x0) / (atol + rtol * np.abs(x1)), np.inf) > 1 and n_it < max_it
        x0 = x1
    return x1, n_it

def quasiNewtonApprossimato(atol, rtol, max_it, x0, J, f):
    n_it = 0
    errore = True
    A=J(x0)
    a=np.linalg.norm(A)
    while errore:
        dx =f(x0)/a
        x1 = x0 - dx
        n_it += 1
        errore = np.linalg.norm(np.abs(x1 - x0) / (atol + rtol * np.abs(x1)), np.inf) > 1 and n_it < max_it
        df = f(x1) - f(x0)
        dx = x1 - x0
        A = A + np.outer((df - A @ dx), dx) / np.dot(dx, dx)
        x0 = x1
    return x1, n_it





x=np.array([0,0])
x,n_it=newton(1e-5,1e-5,1000,x,Jacobiana,f1)
print("NEWTON")
print(x)
print(n_it)
x=np.array([0,0])
x,n_it=quasiNewton(1e-5,1e-5,1000,x,Jacobiana,f1)
print("QNEWTON")
print(x)
print(n_it)
x=np.array([0,0])
x,n_it=quasiNewtonApprossimato(1e-5,1e-5,1000,x,Jacobiana,f1)
print("QANEWTON")
print(x)
print(n_it)