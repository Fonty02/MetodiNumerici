import numpy as np
import matplotlib.pyplot as plt


def interpolazioneDifferenzeDivise(x,f):
    diff_div=np.copy(f).astype(float)
    a=np.copy(f).astype(float)
    righe=np.shape(x)[0]
    for i in range(1,righe):
        diff_div=(diff_div[1:]-diff_div[0:-1])/(x[i:]-x[:-(i)])
        a[i]=diff_div[0]
    return a




def pol(coef, x):
    # Controlla che la lunghezza di coef e x sia coerente
    assert len(coef) == len(x), "Il numero di coefficienti deve essere uguale al numero di elementi in x."

    # Numero di coefficienti
    n = len(coef)

    # Inizializza la stringa della funzione polinomiale
    poly_string = f"{coef[0]}"

    # Aggiungi i termini del polinomio
    for i in range(1, n):
        poly_string += f" + {coef[i]}"
        for j in range(i):
            poly_string += f" * (t - {x[j]})"

    # Valuta la stringa come una funzione lambda
    polynomial = eval(f"lambda t: {poly_string}")
    return polynomial



x = np.array([1, 2, 3, 4])
f = np.array([96, 126, 38, 194])

a = interpolazioneDifferenzeDivise(x, f)

t = np.arange(1, 4, 0.01)
polynomial = pol(a, x)
y = polynomial(t)

plt.plot(t, y, label='Polinomio interpolante')
plt.scatter(x, f, label='Punti di interpolazione')
#nel titolo inserisci i valori di a
plt.title('Interpolazione di Newton\n'+ str(a))
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)
plt.show()
