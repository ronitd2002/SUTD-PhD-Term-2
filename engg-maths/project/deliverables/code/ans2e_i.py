import numpy as np
import matplotlib.pyplot as plt
from scipy.special import spherical_jn, spherical_yn

# Range of ka
ka = np.linspace(0.1, 10.0, 2000)

# Partial waves to plot
ells = [0, 1, 2, 3]

plt.figure(figsize=(8, 5))

for ell in ells:
    jl = spherical_jn(ell, ka)
    nl = spherical_yn(ell, ka)

    # Principal-value phase shift from tan(delta_l) = j_l / n_l
    delta_l = np.arctan2(jl, nl)

    plt.plot(ka, delta_l, label=rf'$\ell={ell}$')

plt.xlabel(r'$ka$')
plt.ylabel(r'$\delta_\ell$')
plt.title(r'Hard-sphere phase shifts $\delta_\ell$ vs $ka$')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()