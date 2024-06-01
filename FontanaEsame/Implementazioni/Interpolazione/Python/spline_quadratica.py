import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import time

def spline_quadratica(x, y, xx):
    n = len(x)
    matrix = np.zeros((n + 1, n + 1))
    coefficienti = np.zeros(n + 1)
    termini_noti = np.zeros(n + 1)
    b = np.zeros(n - 1)
    c = np.zeros(n - 1)

    # Inserisco nella matrice i termini noti delle condizioni di continuità della funzione
    for i in range(n - 1):
        hi = x[i + 1] - x[i]
        termini_noti[i] = y[i + 1] - y[i]
        matrix[i, i] = hi
        matrix[i, i + 2] = hi ** 2
    matrix[0, 2] = 0


    # Inserisco nella matrice i termini noti delle condizioni di continuità della derivata prima
    j = 0
    for i in range(n-1, n+1):
        matrix[i, j] = 1
        matrix[i, j + 1] = -1
        matrix[i, j + 2] = 2 * (x[j + 1] - x[j])
        j += 1
    matrix[n-1, 2] = 0
 

    coefficienti = np.linalg.solve(matrix, termini_noti)
    for i in range(n - 1):
        b[i] = coefficienti[i]
    j = 1
    for i in range(len(coefficienti)-2,len(coefficienti)):
        c[j] = coefficienti[i]
        j += 1
    yy = np.zeros_like(xx)
    for i in range(n - 1):
        index = np.where(np.logical_and(xx >= x[i], xx <= x[i + 1]))
        yy[index] = y[i] + b[i] * (xx[index] - x[i]) + c[i] * (xx[index] - x[i]) ** 2
    return yy



x = np.array([1, 2, 3, 4])
y = np.array([2, 3, 5, 2])
xx = np.linspace(1, 4, 100)
tempo=time.time()
yy = spline_quadratica(x, y, xx)
tempo_finale=time.time()-tempo
plt.plot(x, y, 'ro', label='Punti dati')
plt.plot(xx, yy, 'b-', label='Interpolazione spline quadratica')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Interpolazione spline quadratica con tempo di esecuzione: {tempo_finale:.10f}')
plt.legend()
plt.grid(True)
plt.show()
