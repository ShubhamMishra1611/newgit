# Check Strominger 2105.14346: eq.(2.5) soft-graviton current algebra + rescaling (2.6) => wedge w_{1+inf} (2.7),
# and Jacobi identity of (2.7) on the wedge 1-p<=m<=p-1.
from sympy import Rational as R, factorial, simplify, Symbol
import itertools
kappa=Symbol('kappa')
def Hcomm(k,m,l,n):
    # coefficient c such that [H^k_m, H^l_n] = c * H^{k+l}_{m+n}   (2.5)
    a=R(2-k,2); b=R(2-l,2)
    return -kappa/2*(n*(2-k)-m*(2-l))*factorial(a-m+b-n-1)/(factorial(a-m)*factorial(b-n))*factorial(a+m+b+n-1)/(factorial(a+m)*factorial(b+n))
def N(p,n):  # w^p_n = N(p,n) H^{-2p+4}_n   (2.6)
    return factorial(p-n-1)*factorial(p+n-1)/kappa
ok=True; cnt=0
ps=[R(x,2) for x in range(2,11)]  # p=1,3/2,...,5
for p,qq in itertools.product(ps,ps):
    for m in [p-1-j for j in range(int(2*p-1))]:
        for n in [qq-1-j for j in range(int(2*qq-1))]:
            r=p+qq-2
            if r<1 or abs(m+n)>r-1:   # target outside wedge: skip (coefficient should vanish or not defined)
                continue
            k=-2*p+4; l=-2*qq+4
            lhs=N(p,m)*N(qq,n)*Hcomm(k,m,l,n)/N(r,m+n)
            rhs=(m*(qq-1)-n*(p-1))
            if simplify(lhs-rhs)!=0: ok=False; print('mismatch',p,qq,m,n,lhs,rhs)
            cnt+=1
print('(2.5)+(2.6) -> (2.7) checked on',cnt,'cases; all match:',ok)
# Jacobi for (2.7)
def br(A,B):
    (p,m),(qq,n)=A,B
    return {(p+qq-2,m+n): m*(qq-1)-n*(p-1)}
def br_lin(X,B):
    out={}
    for key,c in X.items():
        for k2,c2 in br(key,B).items(): out[k2]=out.get(k2,0)+c*c2
    return out
bad=0;tot=0
gens=[(p,m) for p in ps[:6] for m in [p-1-j for j in range(int(2*p-1))]]
for A,B,C in itertools.product(gens,repeat=3):
    s={}
    for X,Y,Z in [(A,B,C),(B,C,A),(C,A,B)]:
        for k2,c in br_lin(br(X,Y),Z).items(): s[k2]=s.get(k2,0)+c
    tot+=1
    if any(v!=0 for v in s.values()): bad+=1
print('Jacobi triples tested',tot,'violations',bad)
