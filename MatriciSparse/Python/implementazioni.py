import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.sparse.linalg as spla


def split(A):
    D=sparse.dia_matrix((A.diagonal(),[0]),shape=(A.shape[0],A.shape[1])).tocsr()
    E=sparse.tril(A,k=-1)
    F=sparse.triu(A,k=1)
    return D,E,F


def Jacobi(A:sparse._matrix,b,x0,tol=1e-15,max_iter=5000):
    D,E,F=split(A)
    M=D
    N=-(E+F)
    it=0
    stop=False
    while it<max_iter and not stop:
        x1=spla.spsolve(M,N@x0+b)
        if np.linalg.norm(x1-x0)<tol:
            stop=True
            break
        x0=x1
        it+=1
    if stop:
        return x1
    else:
        print('Il metodo non converge')
        return None
    

def SORForward(A:sparse._matrix,b,x0,tol=1e-15,max_iter=5000,omega=1.0):
    D,E,F=split(A)
    M=D-omega*E
    N=-(omega*F+(1-omega)*D)
    b=omega*b
    it=0
    stop=False
    while it<max_iter and not stop:
        x1=spla.spsolve(M,N @ x0 + b)
        if np.linalg.norm(x1-x0)<tol:
            stop=True
            break
        x0=x1
        it+=1
    if stop:
        return x1
    else:
        print('Il metodo non converge')
        return None


def SORBackward(A:sparse._matrix,b,x0,tol=1e-15,max_iter=5000):
    D,E,F=split(A)
    M=D-F
    N=-E
    it=0
    stop=False
    while it<max_iter and not stop:
        x1=spla.spsolve(M,N@x0+b)
        if np.linalg.norm(x1-x0)<tol:
            stop=True
            break
        x0=x1
        it+=1
    if stop:
        return x1
    else:
        print('Il metodo non converge')
        return None

def SORSymmetric(A:sparse._matrix,b,x0,tol=1e-15,max_iter=5000,omega=1.0):
    D,E,F=split(A)
    M=D-omega*E
    N=-(omega*F+(1-omega)*D)
    b=omega*b
    it=0
    stop=False
    while it<max_iter and not stop:
        x1=spla.spsolve(M,N @ x0 + b)
        x2=spla.spsolve(M,N @ x1 + b)
        if np.linalg.norm(x1-x0)<tol:
            stop=True
            break
        x0=x2
        it+=1
    if stop:
        return x1
    else:
        print('Il metodo non converge')
        return None
    


A=sparse.csr_matrix(np.array([[4,0,0],[0,4,0],[0,0,4]]))
b=np.array([1,2,3])
x0=np.random.rand(3)
print("Jacobi",Jacobi(A,b,x0))
print("SORForward",SORForward(A,b,x0))
print("SORBackward",SORBackward(A,b,x0))
print("SORSymmetric"SORSymmetric(A,b,x0))

print("\n\n",spla.spsolve(A,b))