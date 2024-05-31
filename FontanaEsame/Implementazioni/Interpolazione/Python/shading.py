import numpy as np
import matplotlib.pyplot as plt

# Definizione dei vettori
v1 = np.array([96, 147, 141])
v2 = np.array([126, 207, 242])
v3 = np.array([38, 245, 202])
v4 = np.array([194, 18, 89])

# Creazione di un'immagine vuota
img = np.zeros((100, 100, 3), dtype=np.float64)

# Assegnazione dei colori agli angoli
img[0, 0, :] = v1
img[0, -1, :] = v2
img[-1, 0, :] = v3
img[-1, -1, :] = v4

# Calcolo dei vettori interpolati per le colonne
x = np.arange(100)  # Array da 0 a 99
a = v1 * (1 - x / 99)[:, np.newaxis] + v3 * (x / 99)[:, np.newaxis] #creo nuova dimensione per poter effettuare il broadcasting
b = v2 * (1 - x / 99)[:, np.newaxis] + v4 * (x / 99)[:, np.newaxis]

# Interpolazione lungo le righe
y = np.arange(100)  # Array da 0 a 99
img = a * (1 - y / 99)[:, np.newaxis, np.newaxis] + b * (y / 99)[:, np.newaxis, np.newaxis]

# Visualizzazione dell'immagine
plt.imshow(img.astype(np.uint8))
plt.axis('off')
plt.show()
