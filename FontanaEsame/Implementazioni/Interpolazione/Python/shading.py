import numpy as np
import matplotlib.pyplot as plt
import time

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
tempo_iniziale = time.time()
x = np.arange(100)  # Array da 0 a 99
a = img[0, 0, :] * ((x - 99) / -99)[:, np.newaxis] + img[-1, 0, :] * (x / 99)[:, np.newaxis]
b = img[0, -1, :] * ((x - 99) / -99)[:, np.newaxis] + img[-1, -1, :] * (x / 99)[:, np.newaxis]

# Interpolazione lungo le righe
y = np.arange(100)  # Array da 0 a 99
img = a * ((y - 99) / -99)[:, np.newaxis, np.newaxis] + b * (y / 99)[:, np.newaxis, np.newaxis]
#Al contrario di MatLab non è possibile fare il broadcast di un array 1D con uno 3D, quindi è necessario aggiungere due dimensioni
tempo_finale = time.time()-tempo_iniziale
# Visualizzazione dell'immagine
plt.imshow(img.astype(np.uint8))
plt.axis('off')
plt.title("Tempo di esecuzione: "+str(tempo_finale))
plt.show()