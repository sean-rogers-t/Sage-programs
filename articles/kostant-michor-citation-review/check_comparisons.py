"""Exact checks of comparisons in the literature report, not the E8 proof.

Run with Python and SymPy; output is written next to this script.
"""
import json
from pathlib import Path
import sympy as s

checks = {}
t = s.symbols('t')
I3 = s.eye(3)
X = s.diag(t, t, -2*t)
C = (I3 + X) * (I3 - X).inv()
expected = (1+t)**2*(1-2*t)/((1-t)**2*(1+2*t))
assert s.cancel(C.det()-expected) == 0
assert s.cancel(C.det()-1) != 0
checks['sl3_classical_cayley_determinant'] = str(s.factor(C.det()))

J = s.zeros(4)
J[:2, 2:] = s.eye(2)
J[2:, :2] = -s.eye(2)
g = s.diag(2, -s.Rational(1, 2), s.Rational(1, 2), -2)
assert g.T*J*g == J
theta = (g-g.inv())/2
assert theta == s.diag(s.Rational(3, 4), s.Rational(3, 4), -s.Rational(3, 4), -s.Rational(3, 4))
assert theta.T*J + J*theta == s.zeros(4)
a,b = g[0,0],g[1,1]
positive_root_values = [a/b,a*b,a*a,b*b]
assert all(v != 1 for v in positive_root_values)
positive_root_differentials = [theta[0,0]-theta[1,1],theta[0,0]+theta[1,1],2*theta[0,0],2*theta[1,1]]
assert positive_root_differentials[0] == 0
checks['sp4_source_positive_root_values'] = list(map(str,positive_root_values))
checks['sp4_target_positive_root_values'] = list(map(str,positive_root_differentials))

S = s.diag(2,2,s.Rational(1,4))
U = s.eye(3)
U[0,1] = 1
def theta_sl3(A): return A-s.trace(A)*s.eye(3)/3
assert S.det() == U.det() == 1 and S*U == U*S
N = s.zeros(3); N[0,1]=1
assert theta_sl3(S*U)-theta_sl3(S) == 2*N
assert theta_sl3(U) == N
assert N*N == s.zeros(3) and theta_sl3(S)*N == N*theta_sl3(S)
checks['jordan_nilpotent_factor'] = {'translated':2,'at_identity':1}

# Independent exact identity check for the harmonic-family scalar on a torus.
I4 = s.eye(4)
gt = s.diag(2,3,s.Rational(1,2),s.Rational(1,3))
D = (I4+gt).det()
Xt = (gt-I4)*(gt+I4).inv()
delta = (I4-Xt).det()
assert D*delta == 16
for m in [1,2,3]:
    left = s.Rational(2,4**(2*m))*D**(m-1)*(gt-I4)*(gt+I4).adjugate()
    right = 2*delta**(-m)*Xt
    assert left == right
checks['symplectic_scalar_identity'] = {'n':2,'m_checked':[1,2,3],'D':str(D),'delta':str(delta)}

result = {'status':'passed','arithmetic':'exact rational and symbolic characteristic zero',
          'scope':'Displayed matrix comparisons only; no claim of full proof re-verification',
          'sympy_version':s.__version__,'checks':checks}
path = Path(__file__).with_name('comparison-checks.json')
path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
