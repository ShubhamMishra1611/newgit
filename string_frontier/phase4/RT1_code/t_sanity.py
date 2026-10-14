import sys; sys.path.insert(0,'.')
from rt1core import *
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 80
# bosonic a=1, m2=-1
for n in [1,2,3,4,5,6,8,10,20]:
    f,k = open_factors(n, 1, -1, 0)
    d,j = Dcrit_level(f,k)
    print('bos n',n, mp.nstr(d,12), j)
for n in [2,3,4,5,6,8,10,20]:
    f,k = open_factors(n, 0, 0, 1)
    d,j = Dcrit_level(f,k)
    print('typeI n',n, mp.nstr(d,12), j)
# exact check: bosonic level 2,3 at D=26 exactly
for n in [2,3]:
    f,k = open_factors(n,1,-1,0)
    ok,r = level_ok(f,k,26)
    print(n, ok, r)
