"""Characteristic-zero checks for the singular-block radial attempt.

The obstruction is an abstract potential, NOT an adjoint character.
No inverse of the singular quadratic matrix is used.
"""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def differential_checks():
    q, y, z = s.symbols('q y z')
    k, c = s.symbols('k c', nonzero=True)
    h, l, Q = [s.Function(n)(y, z) for n in ('h', 'l', 'Q')]
    a, d, b = [s.Function(n)(q, y, z) for n in ('a', 'd', 'b')]
    f = h*q*q/2+l*q+Q
    F, E, mu = s.diff(f, q), s.diff(f, y), z*s.diff(f, z)
    defect = a*F+d*E+b*mu-k*f-c
    # Identity BEFORE setting either F or E to zero. Hence its restriction
    # is valid in the actual incidence ring, including its nilpotents.
    rhs = (h*a+s.diff(F,y)*d+z*s.diff(F,z)*b
           +s.diff(a,q)*F+s.diff(d,q)*E+s.diff(b,q)*mu-k*F)
    assert s.expand(s.diff(defect,q)-rhs) == 0
    # Full matrix identity, ordinary derived and logarithmic central columns.
    xs = (q,y)
    A = s.hessian(f,xs)
    B = s.Matrix([s.diff(mu,x) for x in xs])
    av = s.Matrix([a,d])
    derivative_a = s.Matrix([[s.diff(v,x) for v in (a,d)] for x in xs])
    derivative_b = s.Matrix([s.diff(b,x) for x in xs])
    identity = A*av+B*b+derivative_a*s.Matrix([F,E])+derivative_b*mu-k*s.Matrix([F,E])
    assert all(s.expand(v-s.diff(defect,x)) == 0 for v,x in zip(identity,xs))
    # A root pair's discriminant logarithmic derivative. No logarithm choice.
    t, n = s.symbols('t n')
    disc = (1-t)*(1-1/t)
    assert s.cancel(n*t*s.diff(disc,t)/disc-n*(t+1)/(t-1)) == 0
    return dict(unrestricted_block_derivative_identity=True,
                unrestricted_full_matrix_identity=True,
                external_root_pair_logarithmic_derivative=True)


def control_checks():
    q,y,z,k,c=s.symbols('q y z k c')
    S=z+1/z; U=z-1/z
    f=S*(1+y*q*q)-c/k
    b=k*S/U
    D=lambda v:z*s.diff(v,z)
    assert s.cancel(b*D(f)-k*f-c)==0
    assert s.cancel(f.subs(z,1/z)-f)==0
    # For integral even k this is grad(log(U**k)); use k=6 to check
    # the Laurent inversion-invariant unit without symbolic powers.
    disc=U**6
    assert s.cancel(D(disc)/disc-b.subs(k,6))==0
    assert s.cancel(disc.subs(z,1/z)-disc)==0
    gradients=[s.cancel(s.diff(f,q)/(2*S)),s.cancel(s.diff(f,y)/S)]
    I=s.groebner(gradients,q,y)
    assert list(I)==[q*q,q*y]
    assert I.reduce(q)[1]==q and I.reduce(q*q)[1]==0
    # Independent monomial intersection: (q) intersect (q^2,y).
    left=[(1,0)];right=[(2,0),(0,1)]
    intersection=sorted({tuple(max(a,b) for a,b in zip(u,v)) for u in left for v in right})
    assert intersection==[(1,1),(2,0)]
    mu=D(f)
    assert s.cancel(mu-U-U*y*q*q)==0
    H=s.diff(f,q,2)
    assert H.subs(y,0)==0
    J=s.Matrix([s.diff(f,q),s.diff(f,y)]).jacobian((q,y,z))
    assert J.subs({q:0,y:0})==s.zeros(2,3)
    assert s.cancel(D(U)-S)==0
    # Evaluate the differential identity on the bad component, retaining the
    # coefficient of dlog(z). This is NONzero on the chosen unit open.
    obstruction=s.cancel(U*(k-D(b)))
    assert s.cancel(obstruction-k*S*S/U)==0
    assert s.cancel(b*D(U)-obstruction)==0
    assert s.diff(b,q)==s.diff(b,y)==0
    # Every nonzero target h with h^2+4 != 0 has preimages: z^2-h*z-1=0.
    target=s.symbols('target')
    assert s.discriminant(z*z-target*z-1,z)==target*target+4
    # Never use the determinant stratum to discard normal derivative equations.
    assert s.diff(f.subs(y,0),q)==0
    assert s.cancel(s.diff(f,y).subs(y,0)-S*q*q)==0
    return dict(radial_affine_eigenidentity=True,
        inversion_symmetric_Laurent_potential=True,
        logarithmic_gradient_coefficient=True,
        actual_incidence_ideal=['q^2','q*y'],
        primary_decomposition=['(q)','(q^2,y)'],
        nonzero_nilpotent='q',singular_support='q=y=0',
        quadratic_block_rank_on_bad_support=0,
        target_on_bad_support='z-z^(-1)',
        target_logarithmic_derivative='z+z^(-1)',
        bad_image_contains_open='h*(h^2+4) != 0',
        differentiated_radial_rhs='k*(z+z^(-1))^2/(z-z^(-1))',
        bad_image_dominant=True,normal_equations_must_be_retained=True,
        actual_adjoint_character=False)


def main():
    report=dict(result='passed',engine='SymPy rational identities and monomial intersection',
        differential=differential_checks(),control=control_checks(),
        abstract_radial_identity_suffices=False,
        uniform_adjoint_conjecture_proved=False,adjoint_counterexample=False,
        new_Levi_coverage=False,geometric_argument_verified_by_computer=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'adjoint-singular-radial-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
