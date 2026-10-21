"""Plain-Python/SymPy reconstruction from the 27 literal branching weights.

No Sage or project imports. Save the five full-rank A5+A1 quotient
relations so the separate arithmetic audit can replay them in Sage's
older Python/SymPy environment without repeating the symbolic expansion.
"""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from time import perf_counter
import json
import sympy as sp


def main():
    started=perf_counter();x=sp.symbols('x1:7');z=sp.symbols('z')
    cs=sp.symbols('c2 c3 c4 c5 c6');relations=[]
    for k in (2,5,6,8,9,12):
        f=sp.expand(sum((a+b)**k for a,b in combinations(x,2))+
                    sum((-a+sign*z/2)**k for a in x for sign in (-1,1)))
        symmetric,remainder,elementary=sp.symmetrize(f,*x,formal=True)
        assert remainder==0
        f=sp.expand(symmetric.subs({e:(0 if i==0 else cs[i-1]) for i,(e,_) in enumerate(elementary)}))
        zp=sp.Poly(f,z);assert all(e[0]%2==0 for e,c in zp.terms())
        f=sp.expand(sum(c*(4*cs[0])**(e[0]//2) for e,c in zp.terms()))
        if f:relations.append(str(f))
        print('literal moment',k,round(perf_counter()-started,2),flush=True)
    result=dict(result='passed',sympy_version=sp.__version__,sage_imported=False,
        degrees=[5,6,8,9,12],relations=relations,
        recipe='sum_(i<j)(xi+xj)^k + sum_i[(-xi+z/2)^k+(-xi-z/2)^k]; elementary e1=0; z^2=4c2',
        seconds=round(perf_counter()-started,3),
        source_sha256={Path(__file__).name:sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    p=Path(__file__).parent/'e6-independent-rational-base.json'
    p.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
