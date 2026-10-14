"""
R1 independent check of the BMN fuzzy-sphere spectrum (Dasgupta-Sheikh-Jabbari-Van Raamsdonk hep-th/0205185)
and of the reading of the fermion mass term in DSV eq. (1.3):

 H = R Tr( 1/2 Pi^2 - 1/4 [X_A,X_B]^2 - 1/2 Psi^T gamma^A [X_A,Psi] )                      (1.1)
   + R/2 Tr( (mu/3R)^2 X_i^2 + (mu/6R)^2 X_a^2 + i mu/(4R) Psi^T g123 Psi + i 2mu/(3R) eps X_i X_j X_k )   (1.3)

Literal reading: fermion bilinear coefficient (R/2)(i mu/4R) = i mu/8.   Alternative ("no R/2"): i mu/4.
We build real symmetric SO(9) gamma matrices, expand around X_i = (mu/3R) J_i, and diagonalise.
Units R = 1, mu = 1 (masses in units of mu).
Fermion: L = (i/2) psi^T psi_dot - (1/2) psi^T A psi  => frequencies = |eig(A)| (pairs +-).
DSV Table 1 / sec. 5 (as read): X^a: (2j+1)/6 (deg 6(2j+1), j=0..N-1);
 Y: J/3 (deg 2J+1, J=1..N), (J+1)/3 (deg 2J+1, J=0..N-2), 0 (N^2-1 gauge modes);
 fermions (x4 SU(4)): (J+1/4)/3, J=1/2..N-1/2 ; (J+3/4)/3, J=1/2..N-3/2   (deg 2J+1 each).
"""
import itertools
import numpy as np

s0 = np.eye(2); s1 = np.array([[0, 1], [1, 0]]); s3 = np.array([[1, 0], [0, -1]]); ep = np.array([[0, 1], [-1, 0]])
cands = []
for combo in itertools.product([s0, s1, s3, ep], repeat=4):
    m = combo[0]
    for c in combo[1:]:
        m = np.kron(m, c)
    if np.allclose(m, m.T) and np.allclose(m @ m, np.eye(16)) and not np.allclose(m, np.eye(16)):
        cands.append(m)

def search(chosen, start):
    if len(chosen) == 9:
        return chosen
    for k in range(start, len(cands)):
        c = cands[k]
        if all(np.allclose(c @ g + g @ c, 0) for g in chosen):
            r = search(chosen + [c], k + 1)
            if r:
                return r
    return None

gam = search([], 0)
prod = np.eye(16)
for g in gam:
    prod = prod @ g
print("found 9 real symmetric anticommuting 16x16 gammas; gamma^{1..9} = +-1 :", np.round(prod[0, 0], 6), np.allclose(prod, prod[0, 0]*np.eye(16)))
g123 = gam[0] @ gam[1] @ gam[2]

def su2(N):
    j = (N - 1)/2
    m = np.arange(j, -j - 1, -1)
    Jp = np.zeros((N, N), complex)
    for k in range(1, N):
        Jp[k - 1, k] = np.sqrt(j*(j + 1) - m[k]*(m[k] + 1))
    Jx = (Jp + Jp.conj().T)/2; Jy = (Jp - Jp.conj().T)/(2j); Jz = np.diag(m)
    return [Jx, Jy, Jz]

def herm_basis(N):
    B = []
    for a in range(N):
        E = np.zeros((N, N), complex); E[a, a] = 1; B.append(E)
    for a in range(N):
        for b in range(a + 1, N):
            E = np.zeros((N, N), complex); E[a, b] = E[b, a] = 1/np.sqrt(2); B.append(E)
            E = np.zeros((N, N), complex); E[a, b] = -1j/np.sqrt(2); E[b, a] = 1j/np.sqrt(2); B.append(E)
    return B

def spectrum(N, mu=1.0, ferm_coeff=1/8, sign123=+1):
    J = su2(N); B = herm_basis(N); n = len(B)
    Xv = [mu/3*J[0], mu/3*J[1], mu/3*J[2]] + [np.zeros((N, N))]*6
    # bosonic Hessian by finite differences of the potential (exact quartic: use central differences)
    def V(X):
        v = 0
        for A in range(9):
            for Bb in range(9):
                C = X[A] @ X[Bb] - X[Bb] @ X[A]
                v += -0.25*np.trace(C @ C)
        for i in range(3):
            v += 0.5*(mu/3)**2*np.trace(X[i] @ X[i])
        for a in range(3, 9):
            v += 0.5*(mu/6)**2*np.trace(X[a] @ X[a])
        for (i, j, k), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
            v += 1j*(mu/3)*s*np.trace(X[i] @ X[j] @ X[k])
        return v.real
    dim = 9*n
    def vec2X(x):
        return [Xv[A] + sum(x[A*n + a]*B[a] for a in range(n)) for A in range(9)]
    h = 1e-3
    Hs = np.zeros((dim, dim))
    V0 = V(Xv)
    E = np.eye(dim)
    Vp = [V(vec2X(h*E[p])) for p in range(dim)]
    Vm = [V(vec2X(-h*E[p])) for p in range(dim)]
    for p in range(dim):
        Hs[p, p] = (Vp[p] - 2*V0 + Vm[p])/h**2
        for q in range(p + 1, dim):
            Hs[p, q] = Hs[q, p] = (V(vec2X(h*(E[p] + E[q]))) - V(vec2X(h*(E[p] - E[q])))
                                   - V(vec2X(h*(-E[p] + E[q]))) + V(vec2X(-h*(E[p] + E[q]))))/(4*h**2)
    wB = np.sqrt(np.clip(np.linalg.eigvalsh(Hs), 0, None))
    # fermions: H_F = -1/2 sum_A Psi^T g^A [X_A,Psi] + c * i Psi^T g123 Psi  ;  Psi = sum psi^a T^a
    ad = [np.array([[np.trace(B[a] @ (X @ B[b] - B[b] @ X)) for b in range(n)] for a in range(n)]) for X in Xv[:3]]
    A = sum(-np.kron(gam[i], ad[i]) for i in range(3)) + 2*ferm_coeff*mu*1j*sign123*np.kron(g123, np.eye(n))
    ev = np.linalg.eigvalsh(A)
    wF = np.sort(np.abs(ev))[::2]          # +- pairs -> one oscillator each
    return np.sort(wB), np.sort(wF)

def dsv_table(N):
    b = []
    for j in range(N):
        b += [(2*j + 1)/6]*(6*(2*j + 1))
    for Jt in range(1, N + 1):
        b += [Jt/3]*(2*Jt + 1)
    for Jt in range(0, N - 1):
        b += [(Jt + 1)/3]*(2*Jt + 1)
    b += [0.0]*(N*N - 1)
    f = []
    for k in range(1, 2*N, 2):        # J = k/2 = 1/2..N-1/2
        f += [(k/2 + 0.25)/3]*(4*(k + 1))
    for k in range(1, 2*N - 2, 2):    # J = 1/2..N-3/2
        f += [(k/2 + 0.75)/3]*(4*(k + 1))
    # each complex oscillator = 2 real-fermion pairs?  real 16 N^2 components -> 8 N^2 oscillators
    return np.sort(b), np.sort(f)

for N in [2, 3, 4]:
    bT, fT = dsv_table(N)
    print(f"\nN={N}: DSV-table counts: bosons {len(bT)} (expect 9N^2={9*N*N}), fermion oscillators {len(fT)} (expect 8N^2={8*N*N})")
    print(f"  DSV sum m_B = {bT.sum():.6f}, sum m_F = {fT.sum():.6f}, (2/9)N(8N^2+1) = {2/9*N*(8*N*N+1):.6f}")
    for label, coeff, sg in [("literal (R/2)(i mu/4R)=i mu/8, +g123", 1/8, +1), ("literal, -g123", 1/8, -1),
                             ("no-R/2 reading i mu/4, +g123", 1/4, +1), ("no-R/2 reading, -g123", 1/4, -1)]:
        wB, wF = spectrum(N, ferm_coeff=coeff, sign123=sg)
        okB = np.allclose(wB, bT, atol=2e-3)
        okF = len(wF) == len(fT) and np.allclose(wF, fT, atol=1e-6)
        print(f"  {label:40s}: bosons match DSV (atol 2e-3, zero modes ~sqrt(FD error)): {okB} sum m_B={wB.sum():.4f};  fermions match DSV: {okF};  sum m_F = {wF.sum():.6f}; "
              f"distinct m_F = {sorted(set(np.round(wF, 6)))[:6]}")
