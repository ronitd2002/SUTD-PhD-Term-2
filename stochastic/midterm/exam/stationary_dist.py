import numpy as np

# your matrix
P = np.array([
    [0.2776781881, 0.7223218119, 0],
    [0.1078767679185, 0.344921232320815, 1-0.1078767679185-0.344921232320815],
    [0.1, 0.89, 0.01]
])

# ---------------------------
# 1. Compute large power P^n
# ---------------------------

n = 200
A = np.eye(3)

for _ in range(n):
    A = np.dot(A, P)

print("P^n ≈")
print(A)

print("\nfirst row (stationary distribution approx):")
print(A[0])


# ---------------------------
# 2. Compute stationary vector
#    using correct nullspace
#    (P^T - I)π = 0
# ---------------------------

M = P.T - np.eye(3)

U, S, Vt = np.linalg.svd(M)

pi = Vt[-1]
pi = pi / np.sum(pi)

print("\nstationary distribution from nullspace:")
print(pi)


# ---------------------------
# 3. Verify πP = π
# ---------------------------

print("\nπP =")
print(np.dot(pi, P))