import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.sparse.linalg as spla
import matplotlib.pylab as plt
import time


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
    ME=D-omega*E
    NF=-(omega*F+(1-omega)*D)
    MF=D-omega*F
    NE=-(omega*E+(1-omega)*D)
    b=omega*b
    it=0
    stop=False
    while it<max_iter and not stop:
        x1=spla.spsolve(ME,NF @ x0 + b)
        x2=spla.spsolve(MF,NE @ x1 + b)
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



def Met_PotenzeNorm(u0,A,tol=1e-15,it_max=10000):
    n_it = 0
    u1 = A@u0  #sarebbe xk
    u1 = u1*(1/( np.linalg.norm((u1).todense()))) #Normalizzazione (sarebbe zk)    #ci potrebbe mettere un IF per vedere se abbiamo sparse o dense e fare la norma giusta
    lam0 = u1.T@(A@u1)/(u1.T@u1)
    u0 = u1
    #Collezioniamo in una lista le varie approssimazioni di lambda
    #approx = []   
    #approx.append(lam0)
    approx=lam0
    #Usiamo anche un'altra lista per memorizzare le varie stime dell'errore
    #err = []
    #err.append(1)
    err=1
    while((err>tol) & (n_it <it_max) ):
        u1 = A@u0
        u1 = u1*(1/( np.linalg.norm((u1).todense())))
      
        lam = u1.T@(A@u1)/(u1.T@u1)
        #approx.append(np.asarray(lam))
        approx=lam
        err=abs(lam-lam0)/(1+abs(lam))
        #err.append(abs(lam-lam0)/(1+abs(lam))) #Errore relativo, crea una matrice 1x1 (questo è fatto a scopo didattico)
        lam0 = lam
        u0 = u1
        n_it = n_it+1
        
    return lam,u0,n_it,err, approx




def Met_PotenzeGoogle(u0,A,tol=1e-15,it_max=10000,alfa=0.85):
    n_it = 0
    n=A.shape[0]
    #u1 = alfa*A@u0 + (((1-alfa)/n)*u0.sum()*np.ones(n)).reshape(n,1) SPARSO -> più costoso perchè alla fine tanto diventa pieno
    u1 = alfa*A@u0 + (((1-alfa)/n)*u0.sum()) # PIENO
    u1 = u1*(1/( np.linalg.norm(u1)))
    lam0 = u1.T@(A@u1)/(u1.T@u1)
    u0 = u1
    err=1
    while((err>tol) & (n_it <it_max) ):
        #u1 = alfa*A@u0 + (((1-alfa)/n)*u0.sum()*np.ones(n)).reshape(n,1)
        u1 = alfa*A@u0 + (((1-alfa)/n)*u0.sum())
        u1 = u1*(1/( np.linalg.norm(u1)))
        lam = u1.T@(A@u1)/(u1.T@u1)
        err=abs(lam-lam0)/(1+abs(lam))
        lam0 = lam
        u0 = u1
        n_it = n_it+1
        
    return lam,u0,n_it,err


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


def LGKbidiag(A,b,k): #quello della prof
    #k è il numero di componenti (serve poi per il calcolo dei minimi quadrati)
    m,n=A.shape
    tolb=1e-12 #tolleranza per il calcolo degli elementi diagonali di B (potrebbero essere nulli)
    beta=np.zeros((k+2,1))
    alfa=np.zeros((k+1,1))
    P=sparse.csc_matrix((m,k+2))
    Z=sparse.csc_matrix((n,k+1))
    beta[0]=np.linalg.norm(b,2)
    P[:,0]=b/beta[0] #MATRICE C (b A) -> prima colonna B
    Z[:,0]=A.T@P[:,0]
    alfa[0]=np.linalg.norm((Z[:,0]).todense(),2)
    Z[:,0]=Z[:,0]/alfa[0]
    P[:,1]=A@Z[:,0]-alfa[0][0]*P[:,0]
    beta[1]=np.linalg.norm((P[:,1]).todense(),2)
    P[:,1]=P[:,1]/beta[1]
    for i in range(1,k+1):
        Z[:,i]=A.T@P[:,i]-beta[i][0]*Z[:,i-1]
        alfa[i]=np.linalg.norm((Z[:,i]).todense(),2)
        Z[:,i]=Z[:,i]/alfa[i]
        P[:,i+1]=A@Z[:,i]-alfa[i][0]*P[:,i]
        beta[i+1]=np.linalg.norm((P[:,i+1]).todense(),2)
        P[:,i+1]=P[:,i+1]/beta[i+1]
        if (abs(alfa[i])<tolb or abs(beta[i+1])<tolb):
            break
    km=i-1
    diags=np.zeros((2,km+1))
    diags[0,:]=beta[1:km+2].reshape(km+1)
    diags[1,:]=alfa[0:km+1].reshape(km+1)
    ioff=np.array([-1,0])
    B=sparse.dia_matrix((diags,ioff),shape=(km+2,km+1))
    return P[:,0:km+2],Z[:,0:km+1],B,beta[0:km+2],alfa[0:km+2]



A = sparse.rand(10, 4, density=0.8, format="csr", random_state=42)
b = np.random.rand(10)
P,B,Z=bidiagonalize_LGK(A,b)
A_ricostruita=P@B@Z.T
print(np.linalg.norm(A-A_ricostruita))

P,Z,B,beta,alfa=LGKbidiag(A,b,10)
A_ricostruita=P@B@Z.T
print(np.linalg.norm(A-A_ricostruita.todense()))







'''
A=sparse.csr_matrix(np.array([[4,0,0],[0,4,0],[0,0,4]]))
b=np.array([1,2,3])
x0=np.random.rand(3)
print("Jacobi",Jacobi(A,b,x0))
print("SORForward",SORForward(A,b,x0))
print("SORBackward",SORBackward(A,b,x0))
print("SORSymmetric",SORSymmetric(A,b,x0))

print("\n\n",spla.spsolve(A,b))


n=1000
A = sparse.rand(n, n, density=0.15, format="csr", random_state=42)
#u0 = sparse.rand(n, 1, density=0.35, format="csr", random_state=42)
u0=np.random.rand(n)
t1=time.time()
lam,u0,n_it,err = Met_PotenzeGoogle(u0,A)
t2=time.time()
print('tempo',t2-t1)
print('autovalore',lam)
print('iterate',n_it)
print('errore',err)


err = [err[i].item() for i in range(1,len(err))]
plt.semilogy(err)
plt.show()

A = A # + scipy.sparse.identity(n)
A=A/np.max(A.sum(1),1)
A.sum(1) #creata matrice di adiacenza. Ha ancora i dead end però'''