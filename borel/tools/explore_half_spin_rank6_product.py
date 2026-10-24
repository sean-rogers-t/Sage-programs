"""Independent Spin_12 fiber from the product character, not a weight list."""
from sage.all import GF, PolynomialRing, prod
from monomial_staircase import staircase_length
from time import monotonic
import sys
p=int(sys.argv[1]) if len(sys.argv)>1 else 32003
sigma=int(sys.argv[2]) if len(sys.argv)>2 else 1
r=int(sys.argv[3]) if len(sys.argv)>3 else 6
R=PolynomialRing(GF(p),['z'+str(i+1) for i in range(r)]+['b','v'],order='degrevlex')
z=R.gens()[:r]; b,v=R.gens()[r:]; D=prod(z)
F=(prod(1+x for x in z)+sigma*(-1)**r*prod(1-x for x in z))/2
target=[1,2,3,5,7,11,17,23][:r]
eq=[2*x*F.derivative(x)-F-2**(r-2)*b*y for x,y in zip(z,target)]
eq += [b*b-D,v*b*D-1]
print('Starting product character',p,sigma,r,flush=True); start=monotonic()
gb=R.ideal(eq).groebner_basis(algorithm='libsingular:slimgb')
length=staircase_length([tuple(f.lm().exponents()[0]) for f in gb])
print('Length',length,'seconds',monotonic()-start,flush=True)
