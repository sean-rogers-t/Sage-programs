"""Independent Clifford trace audit of the fundamental-spin fiber excess."""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json

from sage.all import QQ, QuadraticField, PolynomialRing, WeylCharacterRing, matrix, prod
from sage.env import SAGE_VERSION
from verify_spin_cayley_geometry import gamma_matrices

ROOT = Path(__file__).resolve().parents[1]


def original_matrix_checks():
    A = QuadraticField(-1, 'ii'); ii = A.gen(); rows = []
    for n in range(2, 6):
        N = 2*n+1; D = 2**n
        gam = gamma_matrices(N, A); one = matrix.identity(A, D)
        basis = [gam[i]*gam[j]/2 for i,j in combinations(range(N), 2)]
        entries = [B.dict() for B in basis]; d = len(basis)
        def weighted_pair(a, b, weights):
            return sum(weights[i]*v*entries[b].get((j,i), 0)
                       for (i,j),v in entries[a].items())
        gram = matrix(A, d, d, lambda a,b: weighted_pair(a,b,[1]*D))
        assert gram == -QQ(D)/4*matrix.identity(A,d)
        Z = -ii*gam[0]*gam[1]
        assert Z*Z == one and (Z*Z).trace() == D
        for seed in (2, 5):
            us = [ii]+[A(seed+j) for j in range(n-1)]
            g = one
            for j,u in enumerate(us):
                g *= (u+1/u)/2*one+(u-1/u)/(2*ii)*gam[2*j]*gam[2*j+1]
            y = prod(u+1/u for u in us[1:])
            lam = 2*ii*y/D; h = lam*Z
            assert g.is_diagonal() and g.trace() == 0 and lam != 0
            gd = g.diagonal(); hd = h.diagonal()
            for B in basis:
                assert (g*B).trace() == (h*B).trace()
            differential = matrix(A,d,d,lambda a,b: weighted_pair(a,b,gd))
            assert differential.rank() == 4*n
            def centralizer_dimension(diag):
                coeffs = {(i*D+j,b):(diag[i]-diag[j])*v
                          for b,entry in enumerate(entries) for (i,j),v in entry.items()
                          if (diag[i]-diag[j])*v}
                return d-matrix(A,D*D,d,coeffs,sparse=True).rank()
            assert centralizer_dimension(gd) == n
            dimK = (n-1)*(2*n-1)+1
            assert centralizer_dimension(hd) == dimK
            expected = dimK-n
            assert d-differential.rank() == expected+n-2 == n*(2*n-3)
            rows.append(dict(rank=n,spin_dimension=D,seed=seed,target_scalar=str(lam),
                             trace_pairings=d,differential_pairings=d*d,
                             differential_rank=4*n,source_centralizer_dimension=n,
                             target_centralizer_dimension=dimK,
                             proposed_fiber_dimension=expected,
                             actual_tangent_dimension=n*(2*n-3),excess=n-2))
        print('Full Clifford traces and differentials passed at rank',n,flush=True)
    return rows


def symbolic_critical_scheme():
    # Clear only powers of z, which is a unit on the central torus.
    C = PolynomialRing(QQ,'c'); c = C.gen(); F = C.fraction_field()
    R = PolynomialRing(F,['z','y']); z,y = R.gens()
    ideal = R.ideal(z*z+1,(z*z-1)*y-c*z)
    assert ideal == R.ideal(z*z+1,2*y+c*z)
    assert ideal == ideal.radical()
    # Ordinary coordinate Hessian, evaluated using z^2=-1, y=-c*z/2.
    mixed = 1-1/z**2
    yy = R(0); zz = 2*y/z**3+c/z**2
    determinant = yy*zz-mixed**2
    assert ideal.reduce(R(determinant.numerator())) == -4
    # The numerator above is - (z^2-1)^2; denominator z^4 is 1 on the locus.
    assert ideal.reduce(R(determinant.denominator())) == 1
    return dict(critical_ideal='(z^2+1, 2*y+c*z)',reduced=True,
                normal_hessian_determinant=-4,unused_coordinates='n-2',
                generic_target_field='QQ(c)',all_nonzero_complex_targets_proved_in_note=True)


def independent_weights():
    rows = []
    for n in range(2,8):
        W = WeylCharacterRing('B'+str(n),style='coroots')
        weights = W(*([0]*(n-1)+[1])).weight_multiplicities()
        assert len(weights)==2**n and all(m==1 for m in weights.values())
        for seed in (2,5):
            us = [QQ(seed+j) for j in range(n)]
            value = sum(m*prod(u**int(2*w[j]) for j,u in enumerate(us))
                        for w,m in weights.items())
            assert value == (us[0]+1/us[0])*prod(u+1/u for u in us[1:])
        rows.append(dict(rank=n,dimension=2**n,evaluation_points=2))
    return rows


def main():
    original = original_matrix_checks()
    report = dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                  original_matrix_checks=original,symbolic_critical_scheme=symbolic_critical_scheme(),
                  independent_weyl_characters=independent_weights(),
                  total_trace_pairings=sum(r['trace_pairings'] for r in original),
                  total_differential_pairings=sum(r['differential_pairings'] for r in original),
                  all_rank_proof_in_note=True,generic_morse_candidate_refuted=True,
                  generic_reducedness_refuted=False,nonregular_fibers_classified=False)
    sources = ['SPIN_FIBER_EXCESS.md','tools/verify_spin_fiber_excess.py',
               'tools/verify_spin_cayley_geometry.py','tools/verify_spin_polynomial_degree.py',
               'GENERAL_CAYLEY_FIBER_CRITICAL.md']
    report['source_sha256'] = {name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for name in sources}
    (ROOT/'spin-fiber-excess-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k.startswith('total_') or k=='result'},indent=2))


if __name__ == '__main__':
    main()
