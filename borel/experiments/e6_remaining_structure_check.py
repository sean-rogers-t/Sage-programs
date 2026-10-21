"""Independent literal-Cartan check of the full-rank subsystem and P1 step.

Plain Python/SymPy; imports neither Sage nor the main image experiments.
"""
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
C=s.Matrix([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],
            [0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
alpha=[C[:,i] for i in range(6)]
zero=s.zeros(6,1)


def orbit(v):
    found={tuple(v)};queue=[v]
    for w in queue:
        for i in range(6):
            t=w-w[i]*alpha[i]
            if tuple(t) not in found:found.add(tuple(t));queue.append(t)
    return found


def local(A):return all(s.denom(x)%3 for x in A)


def v3(x):
    a,b=map(abs,map(int,s.fraction(x)));v=0
    while a%3==0:a//=3;v+=1
    while b%3==0:b//=3;v-=1
    return v


def main():
    roots=orbit(alpha[0]);weights=orbit(s.eye(6)[:,0]);assert len(roots)==72 and len(weights)==27
    reports=[]
    for kind in ('A5','A2A2A1'):
        if kind=='A5':
            chain=[alpha[i-1] for i in (1,3,4,5,6)];z=s.eye(6)[:,1];retained=(1,3,4,5,6)
        else:
            gamma=-alpha[0]-alpha[1]-2*alpha[2]-2*alpha[3]-2*alpha[4]-alpha[5]
            assert tuple(gamma) in roots and tuple(gamma)==(0,0,-1,1,-1,0)
            chain=[alpha[0],alpha[2],gamma,alpha[4],alpha[5]];z=alpha[1];retained=(1,2,3,5,6)
        gram=lambda u,v:(u.T*C.inv()*v)[0]
        assert tuple(z) in roots and gram(z,z)==2
        for i,a in enumerate(chain):
            assert gram(z,a)==0
            for j,b in enumerate(chain):assert gram(a,b)==(2 if i==j else -1 if abs(i-j)==1 else 0)
        # Solve differences and sum-zero, independently of the main formula.
        coeff=s.Matrix([[int(j==i)-int(j==i+1) for j in range(6)] for i in range(5)]+[[1]*6])
        xs=[s.zeros(6,1) for _ in range(6)]
        for coordinate in range(6):
            values=coeff.inv()*s.Matrix([a[coordinate] for a in chain]+[0])
            for j in range(6):xs[j][coordinate]=values[j]
        RB=s.Matrix.hstack(z,*[x-xs[-1] for x in xs[:-1]])
        WB=s.Matrix.hstack(z,*xs[:-1]);RQ=C.inv()*RB
        assert local(RQ) and local(WB) and v3(RQ.det())==v3(WB.det())==0
        proposed={tuple(xs[i]+xs[j]) for i,j in combinations(range(6),2)}
        proposed|={tuple(-x+sign*z/2) for x in xs for sign in (-1,1)}
        assert proposed==weights or {tuple(-a for a in w) for w in proposed}==weights
        actions=[]
        for node in retained:
            reflection=lambda w:w-w[node-1]*alpha[node-1]
            images=[reflection(x) for x in xs]
            perm=[xs.index(x)+1 for x in images]
            expected={1:(2,1,3,4,5,6),3:(1,3,2,4,5,6),
                      5:(1,2,3,5,4,6),6:(1,2,3,4,6,5)}
            if node==4:assert perm==[1,2,4,3,5,6]
            elif node==2:assert perm==[1,2,3,4,5,6]
            else:assert perm==list(expected[node])
            sign=-1 if kind=='A2A2A1' and node==2 else 1
            assert reflection(z)==sign*z
            actions.append(dict(node=node,x_permutation=perm,z_sign=sign))
        reports.append(dict(case=kind,root_basis_determinant=str(RQ.det()),
                            weight_basis_determinant=str(WB.det()),actions=actions))
    q=s.symbols('q')
    bracket=lambda d:sum(q**j for j in range(d))
    weyl=s.prod(bracket(d) for d in (2,5,6,8,9,12))
    betti={}
    for kind,levi,top,rank in (('A5',(2,3,4,5,6),21,72),
                            ('A2A2A1',(2,3,2,3,2),29,720),
                            ('A2A2',(2,3,2,3),30,1440)):
        polynomial=s.Poly(s.cancel(weyl/s.prod(bracket(d) for d in levi)),q)
        assert polynomial.degree()==top
        coefficients=[int(polynomial.nth(d)) for d in range(top+1)]
        assert sum(coefficients)==rank and min(coefficients)>0 and coefficients==coefficients[::-1]
        betti[kind]=coefficients
    base_path=ROOT/'experiments/e6-independent-rational-base.json'
    base=json.loads(base_path.read_text());assert base['result']=='passed'
    cs=s.symbols('c2 c3 c4 c5 c6')
    gb=s.groebner(list(map(s.sympify,base['relations'])),*cs,order='grlex',domain=s.QQ)
    leads=[f.terms(order='grlex')[0][0] for f in gb.polys]
    bounds=[min(e[i] for e in leads if e[i] and sum(e)==e[i]) for i in range(5)]
    standard=[e for e in product(*(range(b) for b in bounds))
              if not any(all(a>=b for a,b in zip(e,lm)) for lm in leads)]
    degrees=[sum(a*w for a,w in zip(e,(2,3,4,5,6))) for e in standard]
    assert len(standard)==36 and max(degrees)==20
    assert [degrees.count(d)+degrees.count(d-1) for d in range(22)]==betti['A5']
    result=dict(result='passed',cartan_determinant=int(C.det()),roots=72,minimal_weights=27,
                schubert_ranks=betti,
                rational_base_dimension=36,rational_base_top_degree=20,
                finite_basis_method='Pure-power leading terms bound every variable; enumerate every standard monomial inside those bounds.',
                cases=reports,p1_root_pairing=2,
                method='Literal Cartan matrix, solved coordinate equations, reflection orbits and exact lattice determinants; Python/SymPy only.',
                source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for p in (Path(__file__).resolve(),base_path)})
    out=ROOT/'experiments/e6-remaining-structure.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
