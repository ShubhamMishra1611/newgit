"""R1: Firat 2102.03936 eqs. (3a),(5),(6),(7),(8): t^3 coefficient (1/3) exp[(6 v(l,l,l) + 3 pi)/l]  (l~ = -1).
v = (1/2i) Log[ Gamma(-il)^2/Gamma(il)^2 * gam(a1)gam(a2)/(gam(a3)gam(a4)) ],  gam(x)=Gamma(1-x)/Gamma(x).
Test (a) principal-branch Log (as printed), (b) continuous log-gamma branch, at l = arcsinh(1)/pi and l -> infinity.
Minimal-area reference: (1/3) r^-6 with r = 4/(3 sqrt 3): 6561/4096 = 1.6018."""
import mpmath as mp
mp.mp.dps = 30
def gam(x): return mp.gamma(1 - x)/mp.gamma(x)
def v_principal(l1, l2, l3):
    i = 1j
    arg = (mp.gamma(-i*l1)**2/mp.gamma(i*l1)**2 * gam((1 + i*l1 + i*l2 + i*l3)/2)*gam((1 + i*l1 - i*l2 + i*l3)/2)
           / (gam((1 - i*l1 - i*l2 + i*l3)/2)*gam((1 - i*l1 + i*l2 + i*l3)/2)))
    return mp.log(arg)/(2*i)
def lg(x): return mp.loggamma(x)
def lgam(x): return lg(1 - x) - lg(x)
def v_cont(l1, l2, l3):
    i = 1j
    s = (2*lg(-i*l1) - 2*lg(i*l1) + lgam((1 + i*l1 + i*l2 + i*l3)/2) + lgam((1 + i*l1 - i*l2 + i*l3)/2)
         - lgam((1 - i*l1 - i*l2 + i*l3)/2) - lgam((1 - i*l1 + i*l2 + i*l3)/2))
    return s/(2*i)
def coeff(v, l): return mp.exp((6*v + 3*mp.pi)/l)/3
lam = mp.asinh(1)/mp.pi
vp, vc = v_principal(lam, lam, lam), v_cont(lam, lam, lam)
print("lambda* =", mp.nstr(lam, 10))
print("principal Log : v =", mp.nstr(vp, 12), " coeff =", mp.nstr(coeff(mp.re(vp), lam), 8))
print("continuous    : v =", mp.nstr(vc, 12), " coeff =", mp.nstr(coeff(mp.re(vc), lam), 8))
print("v_cont - v_principal in units of pi:", mp.nstr((vc - vp)/mp.pi, 8))
print("Firat eq.(8): 1.62187e8")
print("\nlarge lambda (coefficient with each branch; 6561/4096 =", mp.nstr(mp.mpf(6561)/4096, 8), ")")
for l in [1, 3, 10, 30, 100, 1000, 10000]:
    vp, vc = v_principal(l, l, l), v_cont(l, l, l)
    # branch choice that keeps v continuous AND bounded: shift by the integer multiple of pi nearest to continuous v - principal
    print(f"  l={l:6d}: principal coeff={mp.nstr(coeff(mp.re(vp), l), 8):>12s}   continuous coeff={mp.nstr(coeff(mp.re(vc), l), 8):>12s}"
          f"   v_cont/l={mp.nstr(mp.re(vc)/l, 6)}")
