"""
R2 independent check of E2: optimal lower bound on g3~ = g3 M^2/g2 using the three heavy averages of
Caron-Huot & Van Duong 2011.02957 (eqs. 2.22, 2.23, functional ansatz 3.7):
   g2 = <1/m^4>,  g3 = <(3 - 4 J2/(d-2))/m^6>,  0 = <J2 (2 J2 - (5d-4))/m^8>,  J2 = J(J+d-3), J even, m >= M = 1.
We need min alpha s.t. exists beta with  alpha x^2 + (3-4J2/(d-2)) x^3 + beta J2(2J2-(5d-4)) x^4 >= 0
for x = 1/m^2 in (0,1], all even J.  Dividing by x^2: alpha >= max_{x,J} [ c1 x + c2 x^2 ],
c1 = 4J2/(d-2) - 3,  c2 = -beta J2 (2J2-(5d-4)).  alpha(beta) is convex (max of affine in beta) -> 1D minimisation.
Not a re-run of the authors' LP: exact per-J maximisation of a quadratic on [0,1].
"""
import mpmath as mp
mp.mp.dps = 40

def alpha_of_beta(beta, d, Jmax=2000):
    best = mp.mpf(0)   # J=0 gives c1=-3 -> max 0 at x->0
    for J in range(2, Jmax+1, 2):
        J2 = J*(J+d-3)
        c1 = mp.mpf(4)*J2/(d-2) - 3
        c2 = -beta*J2*(2*J2-(5*d-4))
        # max of c1 x + c2 x^2 on [0,1]
        cands = [mp.mpf(0), c1 + c2]
        if c2 < 0:
            xs = -c1/(2*c2)
            if 0 < xs < 1:
                cands.append(c1*xs + c2*xs**2)
        best = max(best, max(cands))
    return best

def optimal(d):
    # golden-section on convex function
    f = lambda b: alpha_of_beta(b, d)
    a, b = mp.mpf('1e-4'), mp.mpf('1')
    g = (mp.sqrt(5)-1)/2
    c, e = b - g*(b-a), a + g*(b-a)
    fc, fe = f(c), f(e)
    for _ in range(160):
        if fc < fe:
            b, e, fe = e, c, fc
            c = b - g*(b-a); fc = f(c)
        else:
            a, c, fc = c, e, fe
            e = a + g*(b-a); fe = f(e)
    bstar = (a+b)/2
    return f(bstar), bstar

if __name__ == "__main__":
    for d in [4, 5, 10]:
        al, be = optimal(d)
        an = (5*d-2 + mp.sqrt(mp.mpf((d+3)*(319*d**3+76*d**2-292*d+32))/(6*(d+1)*(d+4))))/(2*(d-2))
        print(f"d={d}: alpha*={mp.nstr(al,12)} beta*={mp.nstr(be,10)}  | CHVD analytic eq.3.12 alpha*={mp.nstr(an,12)}"
              f" beta*={mp.nstr((an-mp.mpf(5*d-2)/(d-2))/(2*d*(d-1)),10)}")
    print("d=4 closed form 9/2 + 7/4 sqrt(61/5) =", mp.nstr(mp.mpf(9)/2+mp.mpf(7)/4*mp.sqrt(mp.mpf(61)/5), 12))
    print("Cauchy-Schwarz bound 2(5d-4)/(d-2) at d=4:", 2*(5*4-4)/(4-2))
