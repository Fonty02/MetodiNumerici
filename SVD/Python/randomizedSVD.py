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



A=np.random.randn(100,200)
l=5
U, S, V = randomizedSVD(A,l)
print(U.shape, S.shape, V.shape)
A_reconstructed = U @ S @ V
print(np.linalg.norm(A-A_reconstructed,ord='fro'))