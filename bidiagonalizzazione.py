import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.sparse.linalg as spla


def bidiagonalize_LGK(A, b):
    m, n = A.shape
    P = np.zeros((m, n+1))
    B = np.zeros((n+1, n))
    Z = np.zeros((n, n))
    betas=np.zeros(n+1)
    alphas=np.zeros(n)
    beta1p1=b
    betas[0]=np.linalg.norm(beta1p1)
    zeta=np.zeros(n)
    P[:,0]=beta1p1/betas[0]
    for i in range(n):
        w=A.T@P[:,i]+betas[i]*zeta
        alphas[i]=np.linalg.norm(w)
        zeta=w/alphas[i]
        Z[:,i]=zeta
        y=A@zeta-alphas[i]*P[:,i]
        betas[i+1]=np.linalg.norm(y)
        P[:,i+1]=y/betas[i+1]
    B[0,0]=alphas[0]
    for i in range(1,n):
        B[i,i]=alphas[i]
        B[i,i-1]=betas[i]
    B[n,n-1]=betas[n]
    return P, B, Z



A = sparse.rand(10, 4, density=0.8, format="csr", random_state=42)
b = np.random.rand(10)
P,B,Z=bidiagonalize_LGK(A,b)
A_ricostruita=P@B@Z.T
print(np.linalg.norm(A-A_ricostruita))