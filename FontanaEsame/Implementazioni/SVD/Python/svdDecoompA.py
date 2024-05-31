import os
import numpy as np
import cv2
import matplotlib.pyplot as plt

def svdDecompA():
    myFolder = 'C:/Users/fonta/Desktop/MetodiNumerici/SVD/MatLab/faces94/faces94/female/9336923'
    
    # Cerco tutti e soli i files con estensione .jpg
    theFiles = [f for f in os.listdir(myFolder) if f.endswith('.jpg')]
    N = len(theFiles) 

    S = []
    for i, baseFileName in enumerate(theFiles):
        fullFileName = os.path.join(myFolder, baseFileName)
        print(f'Now reading {fullFileName}')
        
        # Trasformare l'immagine in '2D' passando alla scala di grigio per avere una MATRICE invece di un TENSORE
        fi = cv2.imread(fullFileName, cv2.IMREAD_GRAYSCALE)
        
        # Visualizzare nella stessa finestra tutte le immagini presenti nel database
        plt.subplot(int(np.ceil(np.sqrt(N))), int(np.ceil(np.sqrt(N))), i + 1)
        plt.imshow(fi, cmap='gray')
        
        # Numero di righe e colonne di fi
        ss1, ss2 = fi.shape  # ss1=righe, ss2=colonne
        
        M = ss1 * ss2  # DIMENSIONE
        
        fi = fi.reshape(M, 1).astype(np.float64)  # Ogni faccia diventa un vettore colonna
        
        S.append(fi)  # Creazione matrice delle facce
    
    S = np.hstack(S)
    
    fbar = np.mean(S, axis=1, keepdims=True)  # Immagine media
    
    A = S - fbar  # Matrice A -> ogni immagine sottratta la media
    
    U, _, _ = np.linalg.svd(A, full_matrices=False)  # Siccome mi serve solo la matrice U per il range di A
    rank_A = np.linalg.matrix_rank(A)
    RangeA = U[:, :rank_A]  # Range di A
    X = RangeA.T @ A
    
    # Prendo una faccia della stessa persona (ovvero una colonna di S)
    path_faccia_nota = 'C:/Users/fonta/Desktop/MetodiNumerici/SVD/MatLab/faces94/faces94/female/9336923/9336923.1.jpg'
    faccia_nota = cv2.imread(path_faccia_nota, cv2.IMREAD_GRAYSCALE).reshape(M, 1).astype(np.float64)
    proiezione = RangeA.T @ (faccia_nota - fbar)
    
    # Calcolo la distanza tra la proiezione e ogni colonna di X e prendo la più piccola
    distanza_min = np.inf
    for i in range(N):
        distanza = np.linalg.norm(proiezione - X[:, i:i+1])
        if distanza < distanza_min:
            distanza_min = distanza
    
    plt.figure()
    plt.imshow(faccia_nota.reshape(ss1, ss2).astype(np.uint8), cmap='gray')
    plt.title(f'Distanza minima: {distanza_min}')
    
    # Prendo una faccia di una persona DIVERSA
    path_new_face = 'C:/Users/fonta/Desktop/MetodiNumerici/SVD/MatLab/faces94/faces94/female/anpage/anpage.3.jpg'
    new_face = cv2.imread(path_new_face, cv2.IMREAD_GRAYSCALE).reshape(M, 1).astype(np.float64)
    proiezione_new_face = RangeA.T @ (new_face - fbar)
    distanza_min_new_face = np.inf
    for i in range(N):
        distanza = np.linalg.norm(proiezione_new_face - X[:, i:i+1])
        if distanza < distanza_min_new_face:
            distanza_min_new_face = distanza
    
    plt.figure()
    plt.imshow(new_face.reshape(ss1, ss2).astype(np.uint8), cmap='gray')
    plt.title(f'Distanza minima (nuova faccia): {distanza_min_new_face}')
    
    plt.show()
    
    return A, U, fbar, M, ss1, ss2, N, RangeA

# Eseguire la funzione
A, U, fbar, M, ss1, ss2, N, RangeA = svdDecompA()
