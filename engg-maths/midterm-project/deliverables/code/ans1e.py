import numpy as np
import matplotlib.pyplot as plt

U0 = 1.0
a_values = [1, 5, 100]
q_max = np.sqrt(U0)
eps = 1e-8

def rhs(q, U0=1.0):
    return np.sqrt(U0 / q**2 - 1.0)

def f_even(q, a, U0=1.0):
    return np.tan(a * q) - rhs(q, U0)

def f_odd(q, a, U0=1.0):
    return -1.0 / np.tan(a * q) - rhs(q, U0)

def bisect_root(func, left, right, args=(), tol=1e-12, maxiter=200):
    fL = func(left, *args)
    fR = func(right, *args)

    if not np.isfinite(fL) or not np.isfinite(fR):
        return None
    if fL == 0:
        return left
    if fR == 0:
        return right
    if fL * fR > 0:
        return None

    for _ in range(maxiter):
        mid = 0.5 * (left + right)
        fM = func(mid, *args)

        if not np.isfinite(fM):
            return None

        if abs(fM) < tol or 0.5 * (right - left) < tol:
            return mid

        if fL * fM < 0:
            right = mid
            fR = fM
        else:
            left = mid
            fL = fM

    return 0.5 * (left + right)

def find_even_roots(a, U0=1.0):
    roots = []
    n = 0
    while True:
        left = n * np.pi / a + eps
        right = (n + 0.5) * np.pi / a - eps

        if left >= q_max:
            break
        right = min(right, q_max - eps)

        if left < right:
            root = bisect_root(f_even, left, right, args=(a, U0))
            if root is not None and 0 < root < q_max:
                roots.append(root)

        n += 1
    return roots

def find_odd_roots(a, U0=1.0):
    roots = []
    n = 0
    while True:
        left = (n + 0.5) * np.pi / a + eps
        right = (n + 1.0) * np.pi / a - eps

        if left >= q_max:
            break
        right = min(right, q_max - eps)

        if left < right:
            root = bisect_root(f_odd, left, right, args=(a, U0))
            if root is not None and 0 < root < q_max:
                roots.append(root)

        n += 1
    return roots

def masked_for_plot(y, ycap=10):
    y = np.array(y, dtype=float)
    y[(~np.isfinite(y)) | (y < 0) | (y > ycap)] = np.nan
    return y

for a in a_values:
    even_roots = find_even_roots(a, U0)
    odd_roots = find_odd_roots(a, U0)

    n_even = len(even_roots)
    n_odd = len(odd_roots)
    n_total = n_even + n_odd

    print(f"a = {a}")
    print(f"  even states : {n_even}")
    print(f"  odd states  : {n_odd}")
    print(f"  total states: {n_total}")
    if n_even:
        print("  even q-values:", [f"{r:.8f}" for r in even_roots])
    if n_odd:
        print("  odd  q-values:", [f"{r:.8f}" for r in odd_roots])
    print()

    # Plotting grid
    q = np.linspace(1e-4, q_max - 1e-4, 50000)

    y_rhs = rhs(q, U0)
    y_even = np.tan(a * q)
    y_odd = -1.0 / np.tan(a * q)

    plt.figure(figsize=(9, 6))
    plt.plot(q, masked_for_plot(y_rhs, ycap=10), label=r'$\sqrt{U_0/q^2-1}$')
    plt.plot(q, masked_for_plot(y_even, ycap=10), label=rf'$\tan({a}q)$')
    plt.plot(q, masked_for_plot(y_odd, ycap=10), label=rf'$-\cot({a}q)$')

    # Mark even roots
    for i, r in enumerate(even_roots, start=1):
        yr = rhs(r, U0)
        plt.plot(r, yr, 'o')
        plt.annotate(f'E{i}', (r, yr), xytext=(5, 5), textcoords='offset points')

    # Mark odd roots
    for i, r in enumerate(odd_roots, start=1):
        yr = rhs(r, U0)
        plt.plot(r, yr, 's')
        plt.annotate(f'O{i}', (r, yr), xytext=(5, -12), textcoords='offset points')

    plt.xlim(0, q_max)
    plt.ylim(0, 10)
    plt.xlabel(r'$q$')
    plt.ylabel('function value')
    plt.title(
        rf'$a={a},\ U_0={U0}$'
        + '\n'
        + rf'even states = {n_even}, odd states = {n_odd}, total = {n_total}'
    )
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

plt.show()