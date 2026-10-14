import sys; sys.path.insert(0,'.')
from rt1core import *
from fractions import Fraction as Fr
for n in [2,3,4,5,8,20]:
    f,k=open_facs(n,1,-1,0); print('bos',n,float(Dcrit_factors(f,k)[0]), Dcrit_factors(f,k)[1])
for n in [3,5,7,100]:
    f,k=open_facs(n,0,0,1); d=Dcrit_factors(f,k); print('typeI',n,float(d[0]),d[1])
    f,k=open_facs(n,0,0,0); d=Dcrit_factors(f,k); print('c0 a=0',n,float(d[0]),d[1])
