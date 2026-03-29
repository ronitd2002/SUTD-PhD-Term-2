import numpy as np
import matplotlib.pyplot as plt

# Range for |gamma|a
x = np.linspace(1e-4, 10.0, 4000)

# Scattering length ratio
a_s_over_a = 1.0 - np.tanh(x) / x

# Hard-sphere limit
hard_sphere_limit = np.ones_like(x)

plt.figure(figsize=(8, 5))
plt.plot(x, a_s_over_a, label=r'$\dfrac{a_s}{a}=1-\dfrac{\tanh(|\gamma|a)}{|\gamma|a}$')
plt.plot(x, hard_sphere_limit, '--', label=r'Hard-sphere limit: $a_s/a = 1$')

plt.xlabel(r'$|\gamma|a$')
plt.ylabel(r'$a_s/a$')
plt.title(r'Repulsive spherical well: scattering length')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()