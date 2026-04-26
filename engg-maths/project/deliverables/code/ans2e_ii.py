import numpy as np
import matplotlib.pyplot as plt

# Range for gamma*a
x = np.linspace(1e-3, 3*np.pi, 5000)  # avoid zero

# Scattering length expression
a_s_over_a = 1 - np.tan(x)/x

# Mask divergences (where tan blows up)
a_s_over_a[np.abs(a_s_over_a) > 10] = np.nan

plt.figure(figsize=(8,5))

plt.plot(x, a_s_over_a, label=r'$a_s/a = 1 - \frac{\tan(\gamma a)}{\gamma a}$')

# Mark asymptotes (divergences)
for n in range(1, 4):
    asymptote = (n - 0.5)*np.pi
    plt.axvline(asymptote, linestyle='--', alpha=0.5)

plt.axhline(0, linewidth=1)

plt.xlabel(r'$\gamma a$')
plt.ylabel(r'$a_s/a$')
plt.title(r'$s$-wave scattering length for attractive spherical well')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.show()