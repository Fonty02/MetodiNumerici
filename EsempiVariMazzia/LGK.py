#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 18:46:42 2024

@author: frances
"""
# my implementation of iterative LGK bidiagonalization algorithm
import numpy as np
import scipy.sparse.linalg
import time
import scipy as sp
import scipy.sparse as sparse

def LGK(A,b,k):
    (m,n)= A.shape # find the shape of the matrix
    tolb = 1e-12
    beta = np.zeros((k+2,1))
    alfa = np.zeros((k+1,1))
    alfag = np.zeros((k+1,1))
    cg = np.zeros((k+1,1))
    da = np.zeros((k+1,1))
    du = np.zeros((k+1,1))
    sg = np.zeros((k+1,1))
    gam = np.zeros((k+1,1))
    P = np.zeros((m,k+2)) 
    Z = np.zeros((n,k+1))
    beta[0]=np.linalg.norm(b,2)
    # initialization
    gam[0]=beta[0]
    P[:,0]= b[:,0]/beta[0]
    # alfa and beta are the elements of the bidiagonal matrix
    Z[:,0]= A.T@P[:,0]
    alfa[0] = np.linalg.norm(Z[:,0],2)
    Z[:,0]  = Z[:,0]/alfa[0]
    P[:,1]= A@Z[:,0] - alfa[0]*P[:,0]
    beta[1]= np.linalg.norm(P[:,1],2)
    P[:,1]  = P[:,1]/beta[1]
    alfag[0] = alfa[0]
    for i in range(1,k+1):
        Z[:,i]= A.T@P[:,i]-beta[i]*Z[:,i-1] 
        alfa[i] = np.linalg.norm(Z[:,i],2)
        Z[:,i]  = Z[:,i]/alfa[i]
        P[:,i+1]= A@Z[:,i] -alfa[i]*P[:,i]
        beta[i+1]= np.linalg.norm(P[:,i+1],2)
        P[:,i+1]  = P[:,i+1]/beta[i+1]
        #  Givens orthogonal matrices for the QR factorization with elements: 
        #    cg e sg  at the step  i-1 not optimized
        cg[i-1] = alfag[i-1]/np.sqrt( alfag[i-1]**2 + beta[i]**2)
        sg[i-1] = beta[i]/np.sqrt( alfag[i-1]**2 + beta[i]**2)
        alfag[i] = cg[i-1]*alfa[i]
        # da and  du  elements on the main and upper diagonal of R 
        da[i-1] = cg[i-1]*alfag[i-1]+sg[i-1]*beta[i]
        du[i-1] = sg[i-1]*alfa[i]
        #  residual and  scaled known term  (scaled by cg) of the  bidiagonal least squares problem  
        gam[i]=-gam[i-1]*sg[i-1]
        if (abs(alfa[i])<tolb or abs(beta[i+1])<tolb):
           break
    #  solution of the least squares problem with bidiagonal matrix upper triangular
    km=i-1
    y = np.zeros((km+1,1) )  #last element
    y[km]=gam[km]*cg[km]/da[km]
  
    for i in range(km-1,-1,-1):
        y[i]=(gam[i]*cg[i]-du[i]*y[i+1])/da[i]
        
    # solution of the original least square problem using km steps 
    x=np.dot(Z[:,0:(km+1)],y[0:(km+1)])
   
    return( (x, P[:,0:km+2],Z[:,0:km+1],beta[0:km+2],alfa[0:km+2],da[0:km+1],du[0:km+1],cg[0:km+1],sg[0:km+1],gam[0:km+2]))    
  



def LGKbidiag(A,b,k):
    (m,n)= A.shape # find the shape of the matrix
    tolb = 1e-12
    beta = np.zeros((k+2,1))
    alfa = np.zeros((k+1,1))
    
    P =  sparse.csr_matrix((m,k+2)) 
    Z =  sparse.csr_matrix((n,k+1))
    beta[0]=np.linalg.norm(b,2)
    # initialization
   # gam[0]=beta[0]
    P[:,0]= b[:,0]/beta[0]
    # alfa and beta are the elements of the bidiagonal matrix
    Z[:,0]= A.T@P[:,0]
    alfa[0] =  (np.linalg.norm( (Z[:,0]).todense() ,2))
    Z[:,0]  = Z[:,0]/alfa[0]
    P[:,1]= A@Z[:,0] - alfa[0][0]*P[:,0]
    beta[1]= np.linalg.norm(P[:,1].todense(),2)
    P[:,1]  = P[:,1]/beta[1]

    for i in range(1,k+1):
        Z[:,i]= A.T@P[:,i]-beta[i][0]*Z[:,i-1]
        alfa[i] = np.linalg.norm(Z[:,i].todense(),2)
        Z[:,i]  = Z[:,i]/alfa[i]
        P[:,i+1]= A@Z[:,i] - alfa[i][0]*P[:,i]
        beta[i+1]= np.linalg.norm(P[:,i+1].todense(),2)
        P[:,i+1]  = P[:,i+1]/beta[i+1]
        
        if (abs(alfa[i])<tolb or abs(beta[i+1])<tolb):
           break
    km=i-1
    diags = np.zeros((2,(km+1)))
    diags[0,:] = beta[1:(km+2)].reshape(km+1)
    diags[1,:] = alfa[0:(km+1)].reshape(km+1)
    ioff = np.array([-1, 0])
    print(diags)
    B = sparse.dia_matrix((diags, ioff), shape=(km+2, km+1))
    
    return( (P[:,0:km+2],Z[:,0:km+1],B,beta[0:km+2],alfa[0:km+2]))    
    



A = np.array([
[0, 0, 0, 1., 0],
[0, 0, 0, 0, 1.],
[0, 0, 0, 0, 1.],
[1., 0, 1., 0, 0],
[1., 0, 0, 0, 0],
[0, 1., 0, 0, 0],
[1., 0, 1., 1., 0],
[0, 1., 1., 0, 0],
[0, 0, 1., 1., 1.],
[0, 1., 1., 0, 0]])

# query vectors
q1 = np.array([ 
[0],
[0],
[0],
[0],
[0],
[0],
[0],
[1],
[1],
[1]])

q2 = np.array([ 
[0],
[1],
[1],
[0],
[0],
[0],
[0],
[0],
[0],
[0]])
nc = 3
startT1 = time.time()
(x1,P1,Z1,beta1,alfa1,da1,du1,cg1,sg1,gam1) = LGK(A,q1,nc)
endT1 = time.time()
startT = time.time()
As = sparse.lil_matrix(A)
(P1,Z1,B,beta1,alfa1) = LGKbidiag(As,q1,nc)
endT = time.time()
rh = np.zeros((nc+1,1))
rh[0]=beta1[0][0]
(x, istop, itn, r1norm, r2norm,anorm,acond,arnorm,xnorm,var) = sparse.linalg.lsqr(B, rh)
x = Z1*x
print('time q1',endT1-startT1)
print('time q1',endT-startT)
print('bidiagonal matrix for  q1')
print( np.round((P1.T@A)@Z1-B, decimals=3) )
print('B')
print(B)
print(sparse.linalg.norm(B-(P1.T@As)@Z1))
print('beta')
print(beta1)
print('alfa')
print(alfa1)
print( np.round((P1.T@A)@Z1, decimals=3) )
startT = time.time()
(x2,P2,Z2,beta2,alfa2,da2,du2,cg2,sg2,gam2) = LGK(A,q2,nc)
endT = time.time()
startT = time.time()
(P1,Z1,B,beta1,alfa1) = LGKbidiag(As,q2,nc)
endT = time.time()
print('time q2',endT-startT)
print('bidiagonal matrix   for q2')
print( np.round((P2.T@A)@Z2,decimals=3) )
