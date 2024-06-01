import numpy as np
import matplotlib.pyplot as plt
from skimage import io

# Carica l'immagine e la converte in scala di grigi
ref = io.imread('peppers.png', as_gray=True)
ref = np.double(ref)

# Dimensioni dell'immagine
m, n = ref.shape

# Creazione di una maschera per rimuovere alcune informazioni
P = np.random.rand(m, n)
P[P <= 0.5] = 1
P[P != 1] = 0
Rimossa = ref * P

# Trova gli indici dei pixel rimossi
indici = np.argwhere(Rimossa == 0)

# Applica la media dei pixel circostanti ai pixel rimossi
for row, col in indici:
    if row == 0 and col == 0:  # angolo in alto a sinistra
        Rimossa[row, col] = (Rimossa[row + 1, col] + Rimossa[row, col + 1]) / 2
    elif row == 0 and col == n - 1:  # angolo in alto a destra
        Rimossa[row, col] = (Rimossa[row + 1, col] + Rimossa[row, col - 1]) / 2
    elif row == m - 1 and col == 0:  # angolo in basso a sinistra
        Rimossa[row, col] = (Rimossa[row - 1, col] + Rimossa[row, col + 1]) / 2
    elif row == m - 1 and col == n - 1:  # angolo in basso a destra
        Rimossa[row, col] = (Rimossa[row - 1, col] + Rimossa[row, col - 1]) / 2
    elif row == 0 or row == m - 1:  # prima riga o ultima riga
        Rimossa[row, col] = (Rimossa[row, col - 1] + Rimossa[row, col + 1]) / 2
    elif col == 0 or col == n - 1:  # prima colonna o ultima colonna
        Rimossa[row, col] = (Rimossa[row - 1, col] + Rimossa[row + 1, col]) / 2
    else:  # caso generale
        Rimossa[row, col] = (Rimossa[row - 1, col] + Rimossa[row + 1, col] +
                             Rimossa[row, col - 1] + Rimossa[row, col + 1]) / 4

# Mostra l'immagine originale e quella con le informazioni rimosse
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.imshow(ref, cmap='gray')
plt.title('Originale')
plt.subplot(1, 2, 2)
plt.imshow(Rimossa, cmap='gray')
plt.title('Laplace')
plt.show()