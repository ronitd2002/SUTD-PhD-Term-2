import numpy as np
import matplotlib.pyplot as plt

a = 1.0
alpha_values = [1,5,9]

lam = np.linspace(1e-6, 10, 10000)

def even_eq(lam, alpha, a):
    return lam - alpha * (1 + np.exp(-2 * lam * a))

def odd_eq(lam, alpha, a):
    return lam - alpha * (1 - np.exp(-2 * lam * a))

def find_roots(f_vals, lam):
    roots = []
    for i in range(len(lam) - 1):
        if f_vals[i] * f_vals[i+1] < 0:
            root = lam[i] - f_vals[i] * (lam[i+1] - lam[i]) / (f_vals[i+1] - f_vals[i])
            roots.append(root)
    return roots

for alpha in alpha_values:

    f_even = even_eq(lam, alpha, a)
    f_odd = odd_eq(lam, alpha, a)

    even_roots = find_roots(f_even, lam)
    odd_roots = find_roots(f_odd, lam)

    # Energies: E = -lambda^2 (in scaled units)
    E_even = [-r**2 for r in even_roots]
    E_odd = [-r**2 for r in odd_roots]

    print(f"\nalpha = {alpha}")
    print(f"Even states: {len(even_roots)}, Energies: {E_even}")
    print(f"Odd states : {len(odd_roots)}, Energies: {E_odd}")

    # Plot transcendental equations
    plt.figure(figsize=(7,5))
    plt.plot(lam, lam, label=r"$\lambda$")
    plt.plot(lam, alpha*(1 + np.exp(-2*lam*a)), label="Even RHS")
    plt.plot(lam, alpha*(1 - np.exp(-2*lam*a)), label="Odd RHS")

    for r in even_roots:
        plt.plot(r, r, 'o')
    for r in odd_roots:
        plt.plot(r, r, 's')

    plt.xlabel(r"$\lambda$")
    plt.ylabel("Function value")
    plt.title(rf'alpha = {alpha}' + '\n' + 
              rf'Even states: {len(even_roots)}, Energies: {E_even}' + '\n' + rf'Odd states : {len(odd_roots)}, Energies: {E_odd}')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()