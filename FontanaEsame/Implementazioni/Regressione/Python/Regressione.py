import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import time

data = pd.read_csv('winequality-red.csv', delimiter=';')


X = data.iloc[:, :-1]
y = data.iloc[:, -1]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=5)

#Ai fini della regressione lineare, aggiungiamo una colonna di 1 per il termine noto
X_train = np.hstack((np.ones((X_train.shape[0], 1)), X_train.values))
X_test = np.hstack((np.ones((X_test.shape[0], 1)), X_test.values))

# QR
start_time = time.time()
Q, R = np.linalg.qr(X_train)
sol_qr = np.linalg.solve(R, Q.T @ y_train)
qr_time = time.time() - start_time
print("Soluzione con QR:")
print(sol_qr)
print(f'Tempo per QR: {qr_time:.6f} secondi')

# SVD
start_time = time.time()
U, S, VT = np.linalg.svd(X_train, full_matrices=False)
sol_svd = VT.T @ np.linalg.solve(np.diag(S), U.T @ y_train)
svd_time = time.time() - start_time
print('Soluzione con SVD:')
print(sol_svd)
print(f'Tempo per SVD: {svd_time:.6f} secondi')

# Equazioni normali
start_time = time.time()
sol_ne = np.linalg.solve(X_train.T @ X_train, X_train.T @ y_train)
ne_time = time.time() - start_time
print('Soluzione con equazioni normali:')
print(sol_ne)
print(f'Tempo per equazioni normali: {ne_time:.6f} secondi')

# RMSE
y_pred = X_test @ sol_ne
rmse = np.sqrt(np.mean((y_pred - y_test) ** 2))
print(f'RMSE: {rmse:.6f}')
