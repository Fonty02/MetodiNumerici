import numpy as np

def randomizedSVD(A,l):
    m,n = A.shape
    Omega = np.random.randn(n,l)
    Y = A @ Omega
    Q, _ = np.linalg.qr(Y, mode='complete')
    B = Q @ A
    Sigma=np.zeros((m,n))
    U_tilde, S, V = np.linalg.svd(B, full_matrices=True, compute_uv=True)
    U = Q @ U_tilde
    Sigma[:m,:m] = np.diag(S)
    return U, Sigma, V

def randomizedSVD_v2(A,l,q):
   m,n = A.shape
   Omega = np.random.randn(n,l)
   Y_0=A @ Omega
   Q0,R0 = np.linalg.qr(Y_0, mode='complete')
   for i in range(1,q):
    Yi=A.T @ Q0
    Qi,Ri = np.linalg.qr(Yi, mode='complete')
    Y0=A @ Qi
    Q0,R0 = np.linalg.qr(Y0, mode='complete')
   B = Q0 @ A
   Sigma=np.zeros((m,n))
   U_tilde, S, V = np.linalg.svd(B, full_matrices=True, compute_uv=True)
   U = Q0 @ U_tilde
   Sigma[:m,:m] = np.diag(S)
   return U, Sigma, V



A=np.random.randn(100,200)
l=5
U, S, V = randomizedSVD(A,l)
A_reconstructed = U @ S @ V
print(np.linalg.norm(A-A_reconstructed,ord='fro'))
U, S, V = randomizedSVD_v2(A,l,1000)
A_reconstructed = U @ S @ V
print(np.linalg.norm(A-A_reconstructed,ord='fro'))
