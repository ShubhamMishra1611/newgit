"""
R2-T1 step (0): baseline reproduction for the infinite-spin-tower (IST) amplitude of
arXiv:2607.14230 (Berman, Caron-Huot, Chandra, Elvang, Herderschee, Lin, Morales).

Checks (all with kappa^2 = 1, m^2 = 1 unless stated):
  1. Low-energy coefficients of F_IST = k2/(stu) + lam2/((m2-s)(m2-t)(m2-u))   [2607.14230 eq. 5.12]
     reproduce eq. (5.13).
  2. The nonlinear relations listed in 2607.14230 Sec. 5.1 (O(s^7) ... O(s^13)) hold iff lam2 = 2 k2.
     Same relations checked for Virasoro-Shapiro (5.8)/(5.9) as a control, and shown to FAIL for the
     massive-supermultiplet amplitude (5.21) as a negative control.
  3. Residue at s=m^2 of M(zz zbar zbar)=s^4 F equals 8 k2 m^4/(9-z^2) [2610.02124 eq. 106] and its
     Legendre coefficients are (2J+1) Q_J(3)/3 * 8 for even J (positive), reproducing 2610.02124 eq. (107).
  4. Forward dispersive sum rule g_{2k} = 2^{1-k} k2 int dmu/s'^{2k+3} with dmu = delta(s'-m^2)
     [2610.02124 eq. 98, 102] reproduces g_0, g_2, g_4, ...
"""
import sympy as sp
import mpmath as mp

s, t, eps = sp.symbols('s t epsilon')
u = -s - t
ORDER = 13

def basis_decompose(expr_series, kmax=ORDER):
    """Decompose symmetric polynomial (in s,t with u=-s-t), given order by order in eps, into
    sigma3^n1 sigma2^n2.  Returns dict k -> list of ((n1,n2), coeff) ordered by increasing n1."""
    sig2 = sp.expand(s**2 + t**2 + u**2)
    sig3 = sp.expand(s*t*u)
    out = {}
    for k in range(0, kmax + 1):
        poly_k = sp.expand(expr_series.coeff(eps, k))
        pairs = [(n1, (k - 3*n1)//2) for n1 in range(0, k//3 + 1) if (k - 3*n1) % 2 == 0]
        if not pairs:
            assert poly_k == 0, (k, poly_k)
            continue
        cs = sp.symbols('c0:%d' % len(pairs))
        ans = sum(c * sig3**n1 * sig2**n2 for c, (n1, n2) in zip(cs, pairs))
        eqs = sp.Poly(sp.expand(ans - poly_k), s, t).coeffs()
        sol = sp.solve(eqs, cs, dict=True)
        assert len(sol) == 1, (k, sol)
        out[k] = [((n1, n2), sp.nsimplify(sol[0][c])) for c, (n1, n2) in zip(cs, pairs)]
    return out

def label(dec):
    """Map to the paper's labels g_k, g_k', g_k'' (ordered by increasing n1)."""
    g = {}
    for k, lst in dec.items():
        for i, (_, c) in enumerate(lst):
            g[(k, i)] = c
    return g

def check_relations(g, k2):
    G = lambda k, i=0: g.get((k, i), 0)
    rels = {
        'O(s7):  k2 g3 = g0^2/2':                 k2*G(3) - G(0)**2/2,
        'O(s9):  k2 g5 = g2 g0':                   k2*G(5) - G(2)*G(0),
        "O(s10): k4 g6' = 8/3 k4 g6 + g0^3/6":     k2**2*G(6, 1) - sp.Rational(8, 3)*k2**2*G(6) - G(0)**3/6,
        'O(s11): k2 g7 = g4 g0 + g2^2/2':          k2*G(7) - G(4)*G(0) - G(2)**2/2,
        "O(s12): k4 g8' = 8 k4 g8 + g2 g0^2/2":    k2**2*G(8, 1) - 8*k2**2*G(8) - G(2)*G(0)**2/2,
        'O(s13): k2 g9 = g6 g0 + g4 g2':           k2*G(9) - G(6)*G(0) - G(4)*G(2),
        "O(s13): k6 g9' = 8/3 k4 g6 g0 + g0^4/24": k2**3*G(9, 1) - sp.Rational(8, 3)*k2**2*G(6)*G(0) - G(0)**4/24,
    }
    return {k: sp.simplify(v) for k, v in rels.items()}

def series_of(fun):
    e = fun(eps*s, eps*t)
    return sp.series(e, eps, 0, ORDER + 1).removeO()

lam2 = sp.symbols('lambda2', positive=True)
print('=== 1. IST low-energy coefficients (kappa^2=1, m^2=1, general lambda^2) ===')
f_IST = lambda S, T: lam2/((1 - S)*(1 - T)*(1 + S + T))   # F minus graviton pole
gI = label(basis_decompose(series_of(f_IST)))
names = {(0,0):'g0',(2,0):'g2',(3,0):'g3',(4,0):'g4',(5,0):'g5',(6,0):'g6',(6,1):"g6'",(7,0):'g7',(8,0):'g8',(8,1):"g8'",(9,0):'g9',(9,1):"g9'"}
for key, nm in names.items():
    print(f'  {nm:4s} = {gI[key]}')
print('  paper eq.(5.13): g0=l2, g2=l2/2, g3=l2, g4=l2/4, g5=l2, g6=l2/8, g6\'=l2, g7=3 l2/4')

print('\n=== 2. nonlinear relations for IST: residual as function of lambda^2 (kappa^2=1) ===')
for k, v in check_relations(gI, 1).items():
    print(f'  {k:42s} residual = {sp.factor(v)}')
print('  -> every residual vanishes at lambda^2 = 2 (=2 kappa^2), cf. 2607.14230 text after (5.13)')
for k, v in check_relations({kk: vv.subs(lam2, 2) for kk, vv in gI.items()}, 1).items():
    assert v == 0, k

print('\n=== 2b. control: Virasoro-Shapiro via exponential form (5.11), alpha\'=1 ===')
z = sp.symbols('zeta3 zeta5 zeta7 zeta9 zeta11 zeta13')
zeta = dict(zip([3, 5, 7, 9, 11, 13], z))
S_, T_ = eps*s, eps*t
U_ = -S_ - T_
expo = 2*sum(zeta[n]/n*(S_**n + T_**n + U_**n) for n in [3, 5, 7, 9, 11, 13])
# F = (1/stu) exp(expo); f = F - 1/stu.  expo is divisible by stu, so divide order by order exactly.
num = sp.expand(sp.series(sp.exp(expo) - 1, eps, 0, ORDER + 4).removeO())
ser = 0
for k in range(0, ORDER + 1):
    ck = sp.expand(num.coeff(eps, k + 3))
    q, rmd = sp.div(sp.Poly(ck, s, t), sp.Poly(sp.expand(s*t*(-s - t)), s, t))
    assert rmd.is_zero
    ser += eps**k*q.as_expr()
gV = label(basis_decompose(sp.expand(ser)))
print('  VS g0,g2,g3 =', gV[(0,0)], gV[(2,0)], gV[(3,0)], ' (paper 5.9: 2z3, z5, 2z3^2)')
for k, v in check_relations(gV, 1).items():
    print(f'  {k:42s} residual = {sp.simplify(v)}')

print('\n=== 2c. negative control: massive supermultiplet exchange (5.21), lambda^2=1 ===')
f_SM = lambda S, T: 1/(1 - S) + 1/(1 - T) + 1/(1 + S + T)
gS = label(basis_decompose(series_of(f_SM)))
print('  g0,g2,g3,g4,g5 =', [gS[(k,0)] for k in [0,2,3,4,5]], ' (paper 5.22: 3,1,3,1/2,5/2)')
k2 = sp.symbols('k2', positive=True)
r = check_relations(gS, k2)
first_two = [list(r.values())[0], list(r.values())[1]]
print('  first two relations as equations in kappa^2:', first_two,
      ' -> solutions:', sp.solve(first_two[0], k2), sp.solve(first_two[1], k2), '(incompatible)')

print('\n=== 3. residue and Legendre positivity (D=4) ===')
mp.mp.dps = 30
R = lambda tt: 2/((1 - tt)*(2 + tt))        # -Res_{s=1} of F_IST*s^4 with kappa^2=1, lambda^2=2
for zz in [-0.9, -0.3, 0.0, 0.5, 0.99]:
    tt = (zz - 1)/2
    assert abs(R(tt) - 8/(9 - zz**2)) < 1e-12
print('  -Res_{s=m^2} s^4 F_IST = 8 k2 m^4/(9-z^2)  [matches 2610.02124 eq. 106]')
def legendre_coeff(J):
    return (2*J + 1)/mp.mpf(2)*mp.quad(lambda x: 8/(9 - x**2)*mp.legendre(J, x), [-1, 1])
def closed(J):
    # 1/(z-x) = sum (2J+1) Q_J(z) P_J(x); 8/(9-x^2) = (4/3)[1/(3-x)+1/(3+x)]
    Q = mp.legenq(J, 0, 3, type=3).real
    return 8*(2*J + 1)*Q/3 if J % 2 == 0 else 0
print('   J   numeric a_J          closed (8/3)(2J+1)Q_J(3)   ratio a_{J+2}/a_J')
prev = None
for J in range(0, 21, 2):
    a = legendre_coeff(J); c = closed(J)
    print(f'  {J:2d}  {mp.nstr(a, 12):>20s}  {mp.nstr(c, 12):>20s}   {mp.nstr(a/prev, 6) if prev else "":>10s}')
    assert a > 0 and abs(a - c) < 1e-15*abs(c) + 1e-25
    prev = a
print('  asymptotic ratio (3-2sqrt2)^2 =', mp.nstr((3 - 2*mp.sqrt(2))**2, 8), '(exponential decay of couplings)')
for J in [1, 3, 5]:
    assert abs(legendre_coeff(J)) < 1e-20
print('  odd-J coefficients vanish (only even spins exchanged).')

print('\n=== 4. forward sum rule: g_{2k} = 2^{1-k} k2 / m^{2(2k+3)} with one unit of measure at m^2 ===')
for k, key in [(0, (0,0)), (1, (2,0)), (2, (4,0)), (3, (6,0)), (4, (8,0))]:
    print(f'  g_{2*k}: sum rule {sp.Rational(2)**(1-k)}  vs expansion {gI[key].subs(lam2, 2)}')
    assert sp.Rational(2)**(1-k) == gI[key].subs(lam2, 2)
print('\nALL BASELINE CHECKS PASSED')
