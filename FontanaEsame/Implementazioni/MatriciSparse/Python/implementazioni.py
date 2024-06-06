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


def cerchioGershgorin(M,N):
    T=spla.inv(M)@N
    r,c=T.shape
    # Calcola il valore assoluto di ogni elemento della matrice
    T_abs = abs(T)
    diag_T = T_abs.diagonal()
    # Somma i valori assoluti lungo le righe
    somme_righe = T_abs.sum(axis=1)
    #ad ogni riga sottraggo il valore sulla diagonale
    raggio_r = somme_righe - diag_T
    # Trova il massimo valore
    raggio_r_max = np.max(raggio_r)

    # Calcolo il raggio di Gershgorin per ogni colonna
    somme_colonne = T_abs.sum(axis=0)
    #ad ogni colonna sottraggo il valore sulla diagonale
    raggio_c = somme_colonne - diag_T
    # Trova il massimo valore
    raggio_c_max = np.max(raggio_c)
    # Restituisce intersezione dei due raggi
    return min(raggio_r_max,raggio_c_max)




def Jacobi(A:sparse._matrix,b,x0,tol=1e-15,max_iter=5000):
    D,E,F=split(A)
    M=D
    N=-(E+F)
    if (cerchioGershgorin(M,N)>=1):
        print('Il metodo di Jacobi non converge')
        return None
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
        print('Il metodo non converge in {} iterazioni'.format(max_iter))
        return None
    

def SORForward(A:sparse._matrix,b,x0,tol=1e-15,max_iter=5000,omega=1.0):
    D,E,F=split(A)
    M=D-omega*E
    N=-(omega*F+(1-omega)*D)
    if (cerchioGershgorin(M,N)>=1):
        print('Il metodo di SOR non converge')
        return None
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
        print('Il metodo non converge in {} iterazioni'.format(max_iter))
        return None


def SORBackward(A:sparse._matrix,b,x0,tol=1e-15,max_iter=5000):
    D,E,F=split(A)
    M=D-F
    N=-E
    if (cerchioGershgorin(M,N)>=1):
        print('Il metodo di SOR non converge')
        return None
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
        print('Il metodo non converge in {} iterazioni'.format(max_iter))
        return None

def SORSymmetric(A:sparse._matrix,b,x0,tol=1e-15,max_iter=5000,omega=1.0):
    D,E,F=split(A)
    ME=D-omega*E
    NF=-(omega*F+(1-omega)*D)
    MF=D-omega*F
    NE=-(omega*E+(1-omega)*D)
    if (cerchioGershgorin(ME,NF)>=1 or cerchioGershgorin(MF,NE)>=1):
        print('Il metodo di SOR non converge')
        return None
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
        print('Il metodo non converge in {} iterazioni'.format(max_iter))
        return None



def Met_PotenzeNorm(u0,A,tol=1e-15,it_max=100):
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




def Met_PotenzeGoogle(u0,A,tol=1e-15,it_max=100,alfa=0.85):
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
    tolb = 1e-12
    
    P = np.zeros((m, n+1))
    B = np.zeros((n+1, n))
    Z = np.zeros((n, n))
    beta = np.zeros(n+2)
    alfa = np.zeros(n)
        
    # PRIMO CICLO
    beta[0] = np.linalg.norm(b, 2)
    P[:, 0] = b / beta[0]
    Z[:, 0] = A.T @ P[:, 0]
    alfa[0]=np.linalg.norm(Z[:,0],2)
    Z[:,0]=Z[:,0]/alfa[0]
    P[:, 1] = A @ Z[:, 0] - alfa[0] * P[:, 0]
    beta[1]=np.linalg.norm(P[:,1],2)
    P[:,1]=P[:,1]/beta[1]
    
    
    for i in range(1,n):
        Z[:,i]=A.T@P[:,i]-beta[i]*Z[:,i-1]
        alfa[i]=np.linalg.norm(Z[:,i],2)
        Z[:,i]=Z[:,i]/alfa[i]
        P[:,i+1]=A@Z[:,i]-alfa[i]*P[:,i]
        beta[i+1]=np.linalg.norm(P[:,i+1],2)
        P[:,i+1]=P[:,i+1]/beta[i+1]
        if (abs(alfa[i])<tolb or abs(beta[i+1])<tolb):
            break
    km=i-1
    diags=np.zeros((2,km+1))
    diags[0,:]=beta[1:km+2].reshape(km+1)
    diags[1,:]=alfa[0:km+1].reshape(km+1)
    ioff=np.array([-1,0])
    B=sparse.dia_matrix((diags,ioff),shape=(km+2,km+1))
    return P[:,0:km+2],B,Z[:,0:km+1]



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

# Esempio di utilizzo
A = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]], dtype=float)
print("A")
print(A)

b = np.array([1, 2, 3], dtype=float)
x0 = np.random.rand(3)

print("A COME MATRICE SPARSA")
A_sparse = sparse.csr_matrix(A)
tic = time.time()
soluzione = Jacobi(A_sparse, b, x0)
tempo = time.time() - tic
print('Jacobi')
print(soluzione)
print('tempo')
print(tempo)

tic = time.time()
soluzione = SORForward(A_sparse, b, x0)
tempo = time.time() - tic
print('SORForward')
print(soluzione)
print('tempo')
print(tempo)

tic = time.time()
soluzione = SORBackward(A_sparse, b, x0)
tempo = time.time() - tic
print('SORBackward')
print(soluzione)
print('tempo')
print(tempo)

tic = time.time()
soluzione = SORSymmetric(A_sparse, b, x0)
tempo = time.time() - tic
print('SORSymmetric')
print(soluzione)
print('tempo')
print(tempo)

tic = time.time()
soluzione = spla.spsolve(A_sparse, b)
tempo = time.time() - tic
print('Solver classico')
print(soluzione)
print('tempo')
print(tempo)

A = A_sparse.toarray()

tic = time.time()
soluzione = Jacobi(A, b, x0)
tempo = time.time() - tic
print('Jacobi')
print(soluzione)
print('tempo')
print(tempo)

tic = time.time()
soluzione = SORForward(A, b, x0)
tempo = time.time() - tic
print('SORForward')
print(soluzione)
print('tempo')
print(tempo)

tic = time.time()
soluzione = SORBackward(A, b, x0)
tempo = time.time() - tic
print('SORBackward')
print(soluzione)
print('tempo')
print(tempo)

tic = time.time()
soluzione = SORSymmetric(A, b, x0)
tempo = time.time() - tic
print('SORSymmetric')
print(soluzione)
print('tempo')
print(tempo)

tic = time.time()
soluzione = np.linalg.solve(A, b)
tempo = time.time() - tic
print('Solver classico')
print(soluzione)
print('tempo')
print(tempo)

print("MATRICE RANDOMICA A SPARSA GRANDE 100x100 con 90% di zeri")
A = sparse.random(100, 100, density=0.1, format='csr')
b = np.random.rand(100)
x0 = np.random.rand(100)

tic = time.time()
soluzione = spla.spsolve(A, b)
tempo = time.time() - tic
print('Solver classico')
print('tempo')
print(tempo)

print("STESSA MATRICE A MA PIENA")
A = A.toarray()

tic = time.time()
soluzione = np.linalg.solve(A, b)
tempo = time.time() - tic
print('Solver classico')
print('tempo')
print(tempo)

print("MATRICE RANDOMICA A SPARSA GRANDE 100x100 con 40% di zeri")
A = sparse.random(100, 100, density=0.6, format='csr')
b = np.random.rand(100)
x0 = np.random.rand(100)

tic = time.time()
soluzione = spla.spsolve(A, b)
tempo = time.time() - tic
print('Solver classico')
print('tempo')
print(tempo)

print("STESSA MATRICE A MA PIENA")
A = A.toarray()

tic = time.time()
soluzione = np.linalg.solve(A, b)
tempo = time.time() - tic
print('Solver classico')
print('tempo')
print(tempo)

print("METODO DELLE POTENZE PER MATRICE DI GOOGLE SPARSA")
n = 1000
A = sparse.random(n, n, density=0.85)
u0 = np.random.rand(n)
t1 = time.time()
lam, u , it, err= Met_PotenzeGoogle(u0, A)
t2 = time.time() - t1
print('tempo:', t2)
print('autovalore:', lam)
lam, u = spla.eigs(A, k=1)
print('autovalore:', lam[0])