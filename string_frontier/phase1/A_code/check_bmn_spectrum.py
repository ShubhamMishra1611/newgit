"""Quadratic fluctuation spectrum of the BMN (plane-wave) matrix model about the single fuzzy-sphere vacuum.
Hamiltonian (DSV hep-th/0205185 eq.(1.3), R=1, mu=1):
  H = Tr( 1/2 Pi^2 - 1/4 [X_A,X_B]^2 + 1/2 (1/3)^2 X_i^2 + 1/2 (1/6)^2 X_a^2 + i (1/3) eps_ijk X_i X_j X_k )
      + Tr( -1/2 Psi^T gamma^A [X_A, Psi] + i/8 Psi^T gamma^123 Psi )   [the (R/2) prefactor in DSV (1.3) multiplies i mu/(4R)]
Compare with DSV Table 1 (masses in units of mu):
  SO(6): 1/6 + j/3 (0<=j<=N-1, mult 6(2j+1));  SO(3): 1/3 + j/3 (0<=j<=N-2), j/3 (1<=j<=N), mult (2j+1);
  fermions: 1/4 + j/3 (1/2<=j<=N-3/2), 1/12 + j/3 (1/2<=j<=N-1/2), mult 4(2j+1);
  sum rule: sum(boson masses) = sum(fermion masses) = (2/9) N (8 N^2 + 1).
"""
import numpy as np, itertools
from collections import Counter

def su2(N):
    j = (N-1)/2; m = np.arange(j, -j-1, -1)
    Jp = np.zeros((N, N), complex)
    for k in range(1, N):
        Jp[k-1, k] = np.sqrt(j*(j+1) - m[k]*(m[k]+1))
    Jx = (Jp + Jp.T.conj())/2; Jy = (Jp - Jp.T.conj())/(2j); Jz = np.diag(m).astype(complex)
    return [Jx, Jy, Jz]

def herm_basis(N):
    B = []
    for a in range(N):
        E = np.zeros((N, N), complex); E[a, a] = 1; B.append(E)
    for a in range(N):
        for b in range(a+1, N):
            E = np.zeros((N, N), complex); E[a, b] = E[b, a] = 1/np.sqrt(2); B.append(E)
            E = np.zeros((N, N), complex); E[a, b] = -1j/np.sqrt(2); E[b, a] = 1j/np.sqrt(2); B.append(E)
    return B   # orthonormal: Tr(T^A T^B) = delta

eps = np.zeros((3, 3, 3))
for (i, j, k), s in {(0,1,2):1, (1,2,0):1, (2,0,1):1, (0,2,1):-1, (2,1,0):-1, (1,0,2):-1}.items():
    eps[i, j, k] = s

def V(X):
    v = 0.0
    for A in range(9):
        for B in range(9):
            C = X[A] @ X[B] - X[B] @ X[A]
            v += -0.25*np.trace(C @ C).real
    v += 0.5*(1/3)**2*sum(np.trace(X[i] @ X[i]).real for i in range(3))
    v += 0.5*(1/6)**2*sum(np.trace(X[a] @ X[a]).real for a in range(3, 9))
    v += sum((1j*(1/3)*eps[i, j, k]*np.trace(X[i] @ X[j] @ X[k])).real
             for i in range(3) for j in range(3) for k in range(3) if eps[i, j, k] != 0)
    return v

def gammas9():
    s1 = np.array([[0,1],[1,0]]); s3 = np.array([[1,0],[0,-1]]); e = np.array([[0,1],[-1,0]]); one = np.eye(2)
    mats = []
    for combo in itertools.product([one, s1, s3, e], repeat=4):
        M = combo[0]
        for c in combo[1:]: M = np.kron(M, c)
        if np.allclose(M, M.T) and not np.allclose(M, np.eye(16)): mats.append(M)
    def search(chosen, start):
        if len(chosen) == 9: return chosen
        for k in range(start, len(mats)):
            if all(np.allclose(mats[k] @ c + c @ mats[k], 0) for c in chosen):
                r = search(chosen + [mats[k]], k+1)
                if r: return r
        return None
    return search([], 0)

def spectrum(N, sgn_vac=+1, sgn_f=+1):
    J = su2(N); T = herm_basis(N); n = len(T)
    X0 = [sgn_vac*J[i]/3 for i in range(3)] + [np.zeros((N, N), complex)]*6
    # numerical Hessian in orthonormal real coordinates (central differences)
    dim = 9*n; h = 1e-4
    def Vc(xvec):
        X = [X0[A] + sum(xvec[A*n+k]*T[k] for k in range(n)) for A in range(9)]
        return V(X)
    H = np.zeros((dim, dim)); z = np.zeros(dim); V0 = Vc(z)
    for p in range(dim):
        for q in range(p, dim):
            ep = np.zeros(dim); ep[p] = h; eq = np.zeros(dim); eq[q] = h
            H[p, q] = H[q, p] = (Vc(ep+eq) - Vc(ep-eq) - Vc(-ep+eq) + Vc(-ep-eq))/(4*h*h)
    w = np.linalg.eigvalsh(H); bos = np.sqrt(np.clip(w, 0, None))
    # fermions: H_F = (i/2) Psi^T A Psi ; frequencies = |eig(A)|
    g = gammas9(); g123 = g[0] @ g[1] @ g[2]
    K = []
    for i in range(3):
        Ki = np.zeros((n, n))
        for a in range(n):
            for bb in range(n):
                c = np.trace(T[a] @ (X0[i] @ T[bb] - T[bb] @ X0[i]))   # = i K_ab
                Ki[a, bb] = (c/1j).real
        K.append(Ki)
    A = -sum(np.kron(g[i], K[i]) for i in range(3)) + sgn_f*0.25*np.kron(g123, np.eye(n))
    ev = np.linalg.eigvals(A); fer = np.sort(np.abs(ev.imag))[::2]
    return V0, bos, fer

for N in [2, 3, 4]:
    for sv, sf in [(+1, +1), (+1, -1), (-1, +1), (-1, -1)]:
        V0, bos, fer = spectrum(N, sv, sf)
        if abs(V0) > 1e-9:
            continue   # not the supersymmetric vacuum branch
        nz = bos[bos > 1e-4]
        target = 2/9*N*(8*N**2 + 1)
        print(f"N={N} vac_sign={sv} fermion_sign={sf}: V(vac)={V0:.2e}  #zero modes={np.sum(bos<=1e-4)} (expect N^2-1={N*N-1})")
        print("   boson masses (x12):", sorted(Counter(np.round(12*nz).astype(int)).items()))
        print("   fermion masses (x12):", sorted(Counter(np.round(12*fer).astype(int)).items()))
        print(f"   sum bos = {nz.sum():.6f}, sum fer = {fer.sum():.6f}, DSV (2/9)N(8N^2+1) = {target:.6f}")
