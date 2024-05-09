import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.sparse.linalg as spla


def Jacobi(A,b,x0,tol=1e-10,max_iter=1000):
    M=sparse.diags(A.diagonal())
    N=(sparse.tril(A,-1)+sparse.triu(A,1))
    it=0
    stop=False
    x0=sparse.csr_matrix(x0)
    while it<max_iter and not stop:
        x1=spla.spsolve(M,N@x0+b)
        if spla.norm(x1-x0)<tol:
            stop=True
            break
        x0=x1
        it+=1
    if stop:
        return x1
    else:
        print('Il metodo non converge')
        return None
    

    


#create a 5x5 sparse matrix
A=sparse.csr_matrix(np.random.rand(5,5))
b=np.array([1,2,3,4,5])
x0=np.array([0,0,0,0,0])
print(Jacobi(A,b,x0))


print("\n\n",np.linalg.solve(A,b))