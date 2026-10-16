from pathlib import Path
import sys,json,time
from itertools import combinations
from sage.all import QQ,prod

BASE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(BASE/'papers'/'table_audit'))
from recompute import Calculator

def elementary(values,i):
    return sum(prod(subset) for subset in combinations(values,i))

def complete(values,bound,R):
    h=[R(1)]+[R(0)]*bound
    for v in values:
        for i in range(1,bound+1):h[i]+=v*h[i-1]
    return h

def partition(calc,w,r):
    n=calc.rank;k=n-r;ambient=calc.W.domain();basis=ambient.basis()
    indices=[]
    for a in range(r):
        v=w.action(basis[a]);coords=list(v.to_vector())
        j=next(i for i,c in enumerate(coords) if c)
        indices.append(j+1 if coords[j]>0 else 2*n-j)
    indices.sort()
    dual=sorted(2*n+1-v for v in indices)
    lam=tuple(n+k+1-v+sum(dual[a]+v>2*n+1 for a in range(j))
              for j,v in enumerate(dual))
    if sum(lam)!=w.length():raise ArithmeticError(('Bad partition',n,k,w.reduced_word(),lam))
    return tuple(int(v) for v in lam if v)

def main():
    start=time.time();out=[];wrong_sign_examples=[]
    for n in range(2,8):
        calc=Calculator('C'+str(n));x=[calc.x[0]]+[calc.x[j]-calc.x[j-1] for j in range(1,n)]
        for r in range(1,n+1):
            k=n-r
            for i in range(1,r+1):
                f=calc.R(elementary(x[:r],i));computed=calc.expansion(f,r)
                expected_word=list(range(r-i+1,r+1))
                expected={calc.W.from_reduced_word(expected_word):QQ(1)}
                assert computed==expected,('GL',n,k,i,computed,expected)
                out.append({'n':n,'k':k,'i':i,'factor':'GL','passed':True,
                            'word':expected_word,'partition':partition(calc,next(iter(computed)),r)})
            if not k:continue
            h=complete(x[:r],2*k,calc.R)
            for i in range(1,k+1):
                f=calc.R(elementary([v*v for v in x[r:]],i))
                computed=calc.expansion(f,r)
                terms={partition(calc,w,r):str(c) for w,c in computed.items()}
                expected={}
                for b in range(i+1):
                    lam=tuple(v for v in (2*i-b,b) if v)
                    if len(lam)<=r and lam[0]<=n+k:expected[lam]=str((-1)**(i+b))
                assert terms==expected,('Partition',n,k,i,terms,expected)
                rhs=h[i]**2+2*sum((-1)**j*h[i+j]*h[i-j] for j in range(1,i+1))
                assert calc.expansion(rhs,r)==computed,('Chern',n,k,i)
                literal=h[i]**2+2*sum((-1)**i*h[i+j]*h[i-j] for j in range(1,i+1))
                if calc.expansion(literal,r)!=computed:
                    wrong_sign_examples.append({'n':n,'k':k,'i':i,'correct':{str(lam):c for lam,c in terms.items()},
                        'literal_chern':{str(partition(calc,w,r)):str(c) for w,c in calc.expansion(literal,r).items()}})
                out.append({'n':n,'k':k,'i':i,'factor':'C','passed':True,
                            'basis_classes_checked':len(calc.representatives(r,2*i)),
                            'computed':{str(lam):c for lam,c in terms.items()}})
        print('C'+str(n),'completed',flush=True)
    data={'checks':out,'literal_sign_counterexamples':wrong_sign_examples,'seconds':time.time()-start}
    (Path(__file__).parent/'theorem55_results.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print('PASS',len(out),'checks; literal-sign failures',len(wrong_sign_examples),'seconds',time.time()-start,flush=True)
    for row in out:
        if row['n']==4 and row['factor']=='C':print(row,flush=True)
    if wrong_sign_examples:print('First literal sign counterexample:',wrong_sign_examples[0],flush=True)

if __name__=='__main__':main()
