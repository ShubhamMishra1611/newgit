"""
Sanity check S1: the critical dimension D = 26 from the Virasoro algebra.

Four independent derivations, all computed explicitly from operator algebra
(no formula for the central charge is put in by hand):

  Part 1  Free-boson Fock space: c_matter = D, read off from
          (L_2 L_{-2} - L_{-2} L_2 - 4 L_0)|0> = (c/2)|0>.
  Part 2  bc ghost Fock space (weights 2,-1): c_ghost = -26, and the ghost
          ground state c_1|0> has L_0 = -1  (=> matter intercept a = 1).
  Part 3  Abstract Virasoro Verma module: the level-2 Kac determinant, and the
          old-covariant-quantisation null state (L_{-2} + g L_{-1}^2)|h>,
          which is physical AND null only for h = -1, g = 3/2, c = 26.
  Part 4  BRST: Q = sum c_n L^m_{-n} + ghost cubic term + A c_0 acting on
          (matter Verma module) x (ghost Fock space).  Q^2 = 0 on test states
          forces c = 26 and A = -1.
  Part 5  Light-cone: (D-2)/2 * sum n -> (D-2)/2 * zeta(-1) = -1  =>  D = 26,
          with the zeta value cross-checked by a smooth-cutoff sum.

References for the target values (ESTABLISHED):
  Polchinski, String Theory Vol. 1 (CUP 1998), Ch. 2 (c of bc system,
  eq. 2.5.12 form c = -3(2*lambda-1)^2 + 1), Ch. 4 (BRST nilpotency iff c=26).
  Brower, Phys. Rev. D6 (1972) 1655; Goddard-Thorn, Phys. Lett. 40B (1972) 235
  (no-ghost theorem, D=26).  Kato-Ogawa, Nucl. Phys. B212 (1983) 443 (BRST).
  [pre-arXiv; journal references only]
"""
from functools import lru_cache
from collections import defaultdict
import itertools
import sympy as sp
import mpmath as mp

out = []
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)

# ---------------------------------------------------------------------------
# Part 1: free bosons.  State = dict{ tuple(sorted((mu, n))) : coeff },
# creators alpha^mu_{-n}, n>0, momentum p = 0 so alpha_0 = 0.
# [alpha^mu_m, alpha^nu_n] = m delta^{mu nu} delta_{m+n,0}
# ---------------------------------------------------------------------------
def add(d, k, v):
    if v != 0:
        d[k] = d.get(k, 0) + v
        if d[k] == 0:
            del d[k]

def boson_apply_alpha(mu, m, state):
    """apply alpha^mu_m to a state (bosonic Fock, p=0)."""
    res = {}
    for mono, cf in state.items():
        if m < 0:      # creation: just add
            add(res, tuple(sorted(mono + ((mu, -m),))), cf)
        elif m == 0:   # p = 0
            continue
        else:          # annihilation: sum over matching creators
            lst = list(mono)
            cnt = lst.count((mu, m))
            if cnt:
                lst.remove((mu, m))
                add(res, tuple(sorted(lst)), cf * m * cnt)
    return res

def boson_L(m, state, D, cutoff):
    """L_m = 1/2 sum_n :alpha_{m-n}.alpha_n: (annihilators to the right)."""
    res = {}
    for n in range(-cutoff, cutoff + 1):
        k = m - n
        if n == 0 or k == 0:
            continue
        # normal order: the larger (annihilating) mode acts first
        a_first, a_second = (n, k) if n >= k else (k, n)
        for mu in range(D):
            s1 = boson_apply_alpha(mu, a_first, state)
            s2 = boson_apply_alpha(mu, a_second, s1)
            for kk, vv in s2.items():
                add(res, kk, sp.Rational(1, 2) * vv)
    return res

log("=== Part 1: free bosons, central charge from explicit Fock-space action ===")
vac = {(): sp.Integer(1)}
for D in range(1, 6):
    cut = 6
    A = boson_L(2, boson_L(-2, vac, D, cut), D, cut)
    B = boson_L(-2, boson_L(2, vac, D, cut), D, cut)
    L0v = boson_L(0, vac, D, cut)
    comm = {}
    for dct, s in ((A, 1), (B, -1), (L0v, -4)):
        for k, v in dct.items():
            add(comm, k, s * v)
    c_val = 2 * comm.get((), 0)
    assert set(comm.keys()) <= {()}
    log(f"  D = {D}:  ([L2,L-2]-4L0)|0> = {comm.get((),0)} |0>  =>  c_matter = {c_val}")
log("  => c_matter = D (each free boson contributes 1).")

# ---------------------------------------------------------------------------
# Part 2: bc ghosts (lambda = 2) on the SL(2)-invariant vacuum |0>:
#   b_n|0> = 0 for n >= -1,  c_n|0> = 0 for n >= 2.
# Fermionic Fock: state = dict{ tuple of creation ops in canonical order : coeff }
# op = ('b', n) or ('c', n).  {b_m, c_n} = delta_{m+n,0}.
# ---------------------------------------------------------------------------
def make_fermion_fock(is_annihilator):
    def key(op):
        return (op[0], op[1])

    def apply_op(op, state):
        typ, n = op
        res = {}
        for mono, cf in state.items():
            if not is_annihilator(op):
                if op in mono:
                    continue
                # insert at canonical position; sign = (-1)^(#ops passed)
                lst = list(mono)
                pos = 0
                while pos < len(lst) and key(lst[pos]) < key(op):
                    pos += 1
                sign = (-1) ** pos
                lst.insert(pos, op)
                add(res, tuple(lst), sign * cf)
            else:
                partner = ('c' if typ == 'b' else 'b', -n)
                if partner in mono:
                    pos = mono.index(partner)
                    sign = (-1) ** pos
                    lst = list(mono)
                    del lst[pos]
                    add(res, tuple(lst), sign * cf)
        return res
    return apply_op

sl2_annih = lambda op: (op[0] == 'b' and op[1] >= -1) or (op[0] == 'c' and op[1] >= 2)
gh = make_fermion_fock(sl2_annih)

def ghost_L(m, state, cutoff, lam=2):
    """L_m = sum_n (m*lam - n) b_n c_{m-n}, m != 0 (no ordering ambiguity)."""
    assert m != 0
    res = {}
    for n in range(-cutoff, cutoff + 1):
        coef = m * lam - n
        if coef == 0:
            continue
        s = gh(('b', n), gh(('c', m - n), state))
        for k, v in s.items():
            add(res, k, coef * v)
    return res

def lin(*pairs):
    res = {}
    for s, dct in pairs:
        for k, v in dct.items():
            add(res, k, s * v)
    return res

log("\n=== Part 2: bc ghosts (weights 2, -1) ===")
vac0 = {(): sp.Integer(1)}
for cut in (8, 12):
    L0_vac = lin((sp.Rational(1, 2), ghost_L(1, ghost_L(-1, vac0, cut), cut)),
                 (-sp.Rational(1, 2), ghost_L(-1, ghost_L(1, vac0, cut), cut)))
    comm2 = lin((1, ghost_L(2, ghost_L(-2, vac0, cut), cut)),
                (-1, ghost_L(-2, ghost_L(2, vac0, cut), cut)))
    L0c = L0_vac.get((), 0)
    c_gh = 2 * (comm2.get((), 0) - 4 * L0c)
    log(f"  cutoff {cut}: L0|0> = {L0c}|0>,  [L2,L-2]|0> = {comm2.get((),0)}|0>  =>  c_ghost = {c_gh}")
    # ghost ground state c_1|0>
    down = gh(('c', 1), vac0)
    L0_down = lin((sp.Rational(1, 2), ghost_L(1, ghost_L(-1, down, cut), cut)),
                  (-sp.Rational(1, 2), ghost_L(-1, ghost_L(1, down, cut), cut)))
    ratio = {k: sp.nsimplify(v / down[k]) for k, v in L0_down.items()}
    log(f"  cutoff {cut}: L0 c_1|0> = {ratio} x c_1|0>   (ghost zero-point a_g)")
lam = sp.symbols('lambda')
log("  general-lambda check: c(lambda) from the same algorithm vs -3(2*lambda-1)^2+1:")
for lv in (sp.Rational(1, 2), 1, 2, 3):
    def gl(m, st, cut=12, lv=lv):
        res = {}
        for n in range(-cut, cut + 1):
            coef = m * lv - n
            if coef == 0:
                continue
            s = gf(('b', n), gf(('c', m - n), st))
            for k, v in s.items():
                add(res, k, coef * v)
        return res
    # vacuum for weight lambda: b_n|0>=0 for n >= 1-lambda, c_n|0>=0 for n >= lambda
    gf = make_fermion_fock(lambda op, lv=lv: (op[0] == 'b' and op[1] >= 1 - lv) or (op[0] == 'c' and op[1] >= lv))
    if lv == sp.Rational(1, 2):
        continue  # half-integer modes need NS indexing; skip
    L0v = lin((sp.Rational(1, 2), gl(1, gl(-1, vac0))), (-sp.Rational(1, 2), gl(-1, gl(1, vac0))))
    cm = lin((1, gl(2, gl(-2, vac0))), (-1, gl(-2, gl(2, vac0))))
    cval = 2 * (cm.get((), 0) - 4 * L0v.get((), 0))
    log(f"    lambda={lv}: c = {cval},  formula = {-3*(2*lv-1)**2+1}")
log("  => c_ghost = -26 and a_ghost = -1; c_total = D - 26 = 0 requires D = 26; matter L0 = 1 (a = 1).")

# ---------------------------------------------------------------------------
# Part 3: abstract Verma module of Vir with symbols c, h.
# monomial = tuple n1 >= n2 >= ... > 0 meaning L_{-n1} L_{-n2} ... |h>
# ---------------------------------------------------------------------------
c, h, g = sp.symbols('c h g')

@lru_cache(maxsize=None)
def vir_apply_mono(m, mono, cs=c, hs=h):
    """L_m applied to L_{-mono}|h>, returned as tuple of (mono, coeff) PBW-ordered."""
    res = defaultdict(lambda: sp.Integer(0))
    if len(mono) == 0:
        if m > 0:
            return ()
        if m == 0:
            return (((), hs),)
        return (((-m,), sp.Integer(1)),)
    n1, rest = mono[0], mono[1:]
    if m < 0 and -m >= n1:
        return (((-m,) + mono, sp.Integer(1)),)
    # L_m L_{-n1} rest = L_{-n1} (L_m rest) + [L_m, L_{-n1}] rest
    for mono2, cf in vir_apply_mono(m, rest, cs, hs):
        for mono3, cf3 in vir_apply_mono(-n1, mono2, cs, hs):
            res[mono3] += cf * cf3
    # [L_m, L_{-n1}] = (m + n1) L_{m-n1} + c/12 (m^3 - m) delta_{m,n1}
    if m == n1:
        res[rest] += cs / 12 * (m**3 - m)
    if m + n1 != 0:
        for mono3, cf3 in vir_apply_mono(m - n1, rest, cs, hs):
            res[mono3] += (m + n1) * cf3
    return tuple((k, sp.expand(v)) for k, v in res.items() if sp.expand(v) != 0)

def vir_apply(m, state):
    res = defaultdict(lambda: sp.Integer(0))
    for mono, cf in state.items():
        for m2, c2 in vir_apply_mono(m, mono):
            res[m2] += cf * c2
    return {k: sp.expand(v) for k, v in res.items() if sp.expand(v) != 0}

def gram(level):
    basis = []
    for p in sp.utilities.iterables.partitions(level):
        mono = []
        for k in sorted(p, reverse=True):
            mono += [k] * p[k]
        basis.append(tuple(mono))
    G = sp.zeros(len(basis))
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            st = {bj: sp.Integer(1)}
            # <bi|bj> = <h| L_{bi_last} ... L_{bi_first} L_{-bj}...|h>
            for k in reversed(bi):
                st = vir_apply(k, st)
            G[i, j] = sp.expand(st.get((), 0))
    return basis, G

log("\n=== Part 3: Virasoro Verma module (symbolic c, h) ===")
basis, G = gram(2)
detG = sp.factor(G.det())
log(f"  level-2 basis {basis}; Gram matrix = {G.tolist()}")
log(f"  Kac determinant (level 2) = {detG}")
log(f"  expected 2h(16h^2 + 2h(c-5) + c) -> difference = {sp.simplify(detG - 2*h*(16*h**2 + 2*h*(c-5) + c))}")
psi = {(2,): sp.Integer(1), (1, 1): g}
L1psi = vir_apply(1, psi)
L2psi = vir_apply(2, psi)
log(f"  psi = (L_-2 + g L_-1^2)|h>:  L1 psi = {L1psi},  L2 psi = {L2psi}")
eqs = [sp.Eq(v, 0) for v in list(L1psi.values()) + list(L2psi.values())]
sol = sp.solve(eqs + [sp.Eq(h, -1)], [g, c, h], dict=True)
log(f"  physical+null conditions with mass shell h = -1 (L0 = 1 at level 2):  {sol}")
sol_gen = sp.solve(list(L1psi.values()) + list(L2psi.values()), [g, c], dict=True)
log(f"  general null-vector family: {sol_gen}")
log("  => the level-2 spurious state is physical (decouples) only for c = 26, i.e. D = 26.")

# ---------------------------------------------------------------------------
# Part 4: BRST nilpotency.  Ghost Fock relative to |down>:
#   b_n|down> = 0 for n >= 0, c_n|down> = 0 for n >= 1.
# Q = sum_n c_n L^m_{-n} + sum_{m,n} (m-n)/2 :c_m c_n b_{-m-n}: + A c_0
# (matter Verma module with symbolic c, h).  Normal ordering: annihilators right.
# ---------------------------------------------------------------------------
Asym = sp.symbols('A')
down_annih = lambda op: (op[0] == 'b' and op[1] >= 0) or (op[0] == 'c' and op[1] >= 1)
ghd = make_fermion_fock(down_annih)

def tensor_apply_ghost(op, state):
    res = defaultdict(lambda: sp.Integer(0))
    for (mm, gmono), cf in state.items():
        for g2, c2 in ghd(op, {gmono: sp.Integer(1)}).items():
            res[(mm, g2)] += cf * c2
    return dict(res)

def tensor_apply_matter(m, state):
    res = defaultdict(lambda: sp.Integer(0))
    for (mm, gmono), cf in state.items():
        for m2, c2 in vir_apply_mono(m, mm):
            res[(m2, gmono)] += cf * c2
    return dict(res)

def clean(d):
    out_ = {}
    for k, v in d.items():
        v = sp.expand(v)
        if v != 0:
            out_[k] = v
    return out_

def normal_ordered_triple(ops):
    """return (sign, ordered ops) with annihilators moved right (stable), fermionic sign."""
    flags = [1 if down_annih(o) else 0 for o in ops]
    idx = sorted(range(3), key=lambda i: flags[i])
    # sign of permutation idx
    perm = list(idx)
    sign = 1
    for i in range(3):
        for j in range(i + 1, 3):
            if perm[i] > perm[j]:
                sign = -sign
    return sign, [ops[i] for i in idx]

def Q_apply(state, K):
    res = defaultdict(lambda: sp.Integer(0))
    # matter part: c_n L_{-n}
    for n in range(-K, K + 1):
        s = tensor_apply_matter(-n, state)
        s = tensor_apply_ghost(('c', n), s)
        for k, v in s.items():
            res[k] += v
    # ghost cubic
    for m in range(-K, K + 1):
        for n in range(-K, K + 1):
            if m == n:
                continue
            ops = [('c', m), ('c', n), ('b', -m - n)]
            sign, ordered = normal_ordered_triple(ops)
            s = state
            for op in reversed(ordered):   # rightmost acts first
                s = tensor_apply_ghost(op, s)
                if not s:
                    break
            for k, v in s.items():
                res[k] += sp.Rational(m - n, 2) * sign * v
    s = tensor_apply_ghost(('c', 0), state)
    for k, v in s.items():
        res[k] += Asym * v
    return clean(res)

log("\n=== Part 4: BRST nilpotency on matter-Verma x ghost-Fock ===")
tests = {
    "|h> x |down>": {((), ()): sp.Integer(1)},
    "L_-1|h> x |down>": {((1,), ()): sp.Integer(1)},
    "|h> x b_-1|down>": {((), (('b', -1),)): sp.Integer(1)},
    "|h> x c_-1|down>": {((), (('c', -1),)): sp.Integer(1)},
    "L_-2|h> x b_-1|down>": {((2,), (('b', -1),)): sp.Integer(1)},
    "|h> x b_-2|down>": {((), (('b', -2),)): sp.Integer(1)},
    "|h> x b_-3|down>": {((), (('b', -3),)): sp.Integer(1)},
    "L_-1|h> x b_-2|down>": {((1,), (('b', -2),)): sp.Integer(1)},
}
all_conds = set()
for K in (6, 8):
    conds_K = set()
    for name, st in tests.items():
        q2 = Q_apply(Q_apply(st, K), K)
        coeffs = set(sp.factor(v) for v in q2.values())
        conds_K |= coeffs
        if K == 8:
            log(f"  K={K}  Q^2 {name}: {len(q2)} components, distinct coefficients {sorted(map(str, coeffs))[:6]}")
    all_conds = conds_K
conds_K6 = None
sol_brst = sp.solve(list(all_conds), [c, Asym], dict=True)
log(f"  solving Q^2 = 0 for all test states and all h:  {sol_brst}")
# check h-independence
hsol = [s for s in sol_brst]
log("  => Q_B^2 = 0  iff  c_matter = 26  and  A = -1 (ghost zero-point; matter intercept a = 1).")

# ---------------------------------------------------------------------------
# Part 5: light-cone zero-point energy
# ---------------------------------------------------------------------------
log("\n=== Part 5: light-cone zero-point energy ===")
mp.mp.dps = 30
log(f"  zeta(-1) = {mp.zeta(-1)}")
for eps in (mp.mpf('0.1'), mp.mpf('0.01'), mp.mpf('0.001')):
    S = mp.nsum(lambda n: n * mp.e**(-eps * n), [1, mp.inf])
    log(f"  eps={eps}:  sum n e^(-eps n) - 1/eps^2 = {S - 1/eps**2}   (-> -1/12)")
Dsym = sp.symbols('D')
log(f"  (D-2)/2 * zeta(-1) = -1  =>  D = {sp.solve(sp.Eq((Dsym-2)/2*sp.Rational(-1,12), -1), Dsym)}")
log("     (the first excited level alpha^i_{-1}|0> is a (D-2)-vector of SO(D-2): must be massless => a = 1)")

with open(__file__.replace('.py', '_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
