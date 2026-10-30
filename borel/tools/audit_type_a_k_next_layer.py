"""Native audit of the fixed next-layer obstruction and relative lattices."""
from functools import lru_cache
from hashlib import sha256
from math import comb, prod
from pathlib import Path
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix
from binary_tensor_certificates_native import binary_tensor_certificate_native, cubic_offset_obstruction_native

ROOT = Path(__file__).resolve().parents[1]


@lru_cache(None)
def powers(weights, n):
    if not weights:
        return ((),) if n==0 else ()
    return tuple((a,)+tail for a in range(n//weights[0]+1)
                 for tail in powers(weights[1:],n-a*weights[0]))


def reduction(A,p,transform=False):
    rows=[[int(x)%p for x in row] for row in A.tolist()]
    n,m=A.shape
    U=[[int(i==j) for j in range(n)] for i in range(n)] if transform else None
    pivots=[]
    for j in range(m):
        r=len(pivots)
        pivot=next((i for i in range(r,n) if rows[i][j]),None)
        if pivot is None:
            continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        if transform:
            U[r],U[pivot]=U[pivot],U[r]
        inv=pow(rows[r][j],-1,p)
        rows[r]=[(x*inv)%p for x in rows[r]]
        if transform:
            U[r]=[(x*inv)%p for x in U[r]]
        for i in range(n):
            if i==r or not rows[i][j]:
                continue
            c=rows[i][j]
            rows[i]=[(x-c*y)%p for x,y in zip(rows[i],rows[r])]
            if transform:
                U[i]=[(x-c*y)%p for x,y in zip(U[i],U[r])]
        pivots.append(j)
        if len(pivots)==n:
            break
    return rows,pivots,U


def primary(A,p):
    """Exact first two Smith layers, with rational rank bounding all later ones."""
    rational_rank=len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1])
    modulus=p*p
    rows=[[int(x)%modulus for x in row] for row in A.tolist()]
    n,m=A.shape
    r=0
    while r<min(n,m):
        pivot=next(((i,j) for i in range(r,n) for j in range(r,m) if rows[i][j]%p),None)
        if pivot is None:
            break
        i,j=pivot
        rows[r],rows[i]=rows[i],rows[r]
        for row in rows:
            row[r],row[j]=row[j],row[r]
        inv=pow(rows[r][r],-1,modulus)
        rows[r]=[(x*inv)%modulus for x in rows[r]]
        for i in range(r+1,n):
            c=rows[i][r]
            for j in range(r,m):
                rows[i][j]=(rows[i][j]-c*rows[r][j])%modulus
        r+=1
    residual=s.Matrix([[rows[i][j]//p for j in range(r,m)] for i in range(r,n)])
    assert all(rows[i][j]%p==0 for i in range(r,n) for j in range(r,m))
    next_rank=len(reduction(residual,p)[1]) if residual.rows and residual.cols else 0
    assert r+next_rank==rational_rank
    return [p]*next_rank


def shared_torsion_witness(J,C,p):
    rows,pivots,_=reduction(C,p)
    _,jpivots,U=reduction(J,p,True)
    jrank=len(jpivots)
    columns=[[int(x) for x in C.col(j)] for j in range(C.cols)]
    for free in range(C.cols):
        if free in pivots:
            continue
        numerator=[columns[free][i]-sum(rows[r][free]*columns[j][i]
                                       for r,j in enumerate(pivots)) for i in range(C.rows)]
        assert all(x%p==0 for x in numerator)
        witness=[x//p for x in numerator]
        if any(sum(a*b for a,b in zip(row,witness))%p for row in U[jrank:]):
            # p*witness is in C by the displayed integer column combination,
            # but witness is outside J even modulo p.
            return True
    return False


def moment_sector(p,k):
    d=p+3
    variables=s.symbols('L '+' '.join('a%d'%i for i in range(1,d+1)))
    L=variables[0]
    aa=list(variables[1:])
    qs=[L**j+sum(comb(j,i)*aa[i-1]*L**(j-i) for i in range(1,j+1)) for j in range(1,d+1)]
    basis=[e for e in powers(tuple([1]+list(range(1,d+1))),d) if sum(e[1:])<=k]
    J,C=[],[]
    for e in powers(tuple(list(range(1,d+1))*2),d):
        if sum(e)>k:
            continue
        poly=s.Poly(prod(g**t for g,t in zip(aa+qs,e)),*variables)
        terms=dict(poly.terms())
        column=s.Matrix([terms.get(power,0) for power in basis])
        J.append(column)
        if any(e[:d]):
            C.append(column)
    J=s.Matrix.hstack(*J)
    C=s.Matrix.hstack(*C)
    a,b=primary(J,p),primary(C,p)
    assert a==b==([p] if k==3 else [])
    shared=shared_torsion_witness(J,C,p) if k==3 else None
    assert k!=3 or shared
    return dict(p=p,degree=d,coefficient_degree=k,ambient_rank=len(basis),
                source_primary_torsion=a,ideal_primary_torsion=b,relative_primary_torsion=[],
                shared_torsion_witness=shared,
                method='Exact modular Smith layers through p^2, rational rank, and an integer column-relation witness')


def main():
    tensors=[binary_tensor_certificate_native(k,3,proper) for k in [3,4] for proper in [False,True]]
    assert [t['gluing_obstruction_dimension'] for t in tensors]==[1,1,0,0]
    witness=cubic_offset_obstruction_native()
    moments=[]
    for p in [5,7]:
        for k in [3,4]:
            moments.append(moment_sector(p,k))
            print(json.dumps(moments[-1]),flush=True)
    names=['TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md','tools/audit_type_a_k_next_layer.py',
           'tools/binary_tensor_certificates_native.py']
    report=dict(result='passed',tensor_certificates=tensors,explicit_obstruction=witness,
                moment_sectors=moments,sage_imported=False,external_review=False,
                actual_elementary_degree_p_plus_three_formula=False,
                all_rank_claim_from_finite_checks=False,
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for name in names})
    (ROOT/'type-a-k-next-layer-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(result='passed',moment_sectors=len(moments))))


if __name__=='__main__':
    main()
