"""Independent exact BGG audit of all explicit Part III Schubert expansions.
Run with sage -python recompute.py. No historical calculation code is imported.
"""
from pathlib import Path
from itertools import combinations
import json,time
from sage.all import PolynomialRing, QQ, WeylGroup, prod

ROOT=Path(__file__).resolve().parent

class Calculator:
    def __init__(self,kind):
        self.kind=kind;self.rank=int(kind[1:])
        self.W=WeylGroup(kind,prefix='s');self.s=self.W.simple_reflections()
        roots=self.W.domain().simple_roots()
        self.R=PolynomialRing(QQ,['w'+str(i) for i in range(1,self.rank+1)])
        self.x=self.R.gens()
        self.alpha={i:sum(2*roots[i].dot_product(roots[j])/roots[j].dot_product(roots[j])*self.x[j-1]
                         for j in range(1,self.rank+1)) for i in range(1,self.rank+1)}
        self.levels={}

    def reflection(self,f,i):
        images=list(self.x);images[i-1]-=self.alpha[i]
        return f(*images)

    def divided(self,f,i):
        q,r=(f-self.reflection(f,i)).quo_rem(self.alpha[i])
        if r:raise ArithmeticError('Nonexact divided difference')
        return q

    def representatives(self,parabolic,degree):
        levels=self.levels.setdefault(parabolic,[{self.W.one()}])
        while len(levels)<=degree:
            d=len(levels)
            candidates={s*w for w in levels[-1] for s in self.s if (s*w).length()==d}
            levels.append({w for w in candidates if all((w*self.s[i]).length()>d
                    for i in range(1,self.rank+1) if i!=parabolic)})
        return levels[degree]

    def polynomial(self,row):
        spec=row['input'];ys={int(k):sum(c*x for c,x in zip(v,self.x))
                                for k,v in row.get('y_definitions',{}).items()}
        op=spec['operation']
        if op=='linear':return self.R(ys[spec['variable']])
        if op=='omega':return self.x[spec['index']-1]
        if op=='elementary':
            values=[ys[k]**(2 if spec.get('squared') else 1) for k in spec['variables']]
            return self.R(sum(prod(c) for c in combinations(values,spec['degree'])))
        if op=='psi':
            a=self.alpha
            xx=[5*a[1]+4*a[3]+3*a[4]+2*a[5]+a[6],
                -a[1]+4*a[3]+3*a[4]+2*a[5]+a[6],
                -a[1]-2*a[3]+3*a[4]+2*a[5]+a[6],
                -a[1]-2*a[3]-3*a[4]+2*a[5]+a[6],
                -a[1]-2*a[3]-3*a[4]-4*a[5]+a[6],
                -a[1]-2*a[3]-3*a[4]-4*a[5]-5*a[6]]
            z=-3*(a[1]+2*a[2]+2*a[3]+3*a[4]+2*a[5]+a[6])
            values=[v+z for v in xx]+[v-z for v in xx]+[-u-v for u,v in combinations(xx,2)]
            return self.R(sum(v**spec['degree'] for v in values))
        if op=='g2quadratic':return (self.x[1]-QQ(3)/2*self.x[0])**2
        if op=='g2product':return self.x[0]*self.x[1]-self.x[0]**2
        raise ValueError(spec)

    def expansion(self,f,parabolic):
        degree=f.total_degree();cache={():f};out={}
        for w in self.representatives(parabolic,degree):
            prefix=()
            for i in reversed(w.reduced_word()):
                nxt=prefix+(int(i),)
                if nxt not in cache:cache[nxt]=self.divided(cache[prefix],i)
                prefix=nxt
            value=cache[prefix]
            if value and value.total_degree()!=0:raise ArithmeticError('Nonscalar coefficient')
            if value:out[w]=QQ(value)
        return out

    def audit(self,row):
        f=self.polynomial(row);p=row['parabolic'];d=f.total_degree()
        noninvariant=[i for i in range(1,self.rank+1) if i!=p and f!=self.reflection(f,i)]
        printed={};word_errors=[]
        for term in row['published_terms']:
            word=term['word'];w=self.W.from_reduced_word(word)
            if w.length()!=d or len(word)!=d or w not in self.representatives(p,d):word_errors.append(word)
            printed[w]=printed.get(w,0)+QQ(term['coefficient'])
        computed=self.expansion(f,p)
        differences=[]
        for w in set(printed)|set(computed):
            if printed.get(w,0)!=computed.get(w,0):
                differences.append({'word':[int(i) for i in w.reduced_word()],
                                    'printed':str(printed.get(w,0)),'computed':str(computed.get(w,0))})
        return {'id':row['id'],'degree':int(d),'input_terms':len(f.monomials()),
                'basis_classes_checked':len(self.representatives(p,d)),
                'printed_terms':len(row['published_terms']),
                'noninvariant_nodes':noninvariant,'word_errors':word_errors,
                'differences':differences,'passed':not(noninvariant or word_errors or differences),
                'computed_terms':[{'word':[int(i) for i in w.reduced_word()],'coefficient':str(c)}
                                  for w,c in sorted(computed.items(),key=lambda pair:pair[0].reduced_word())]}

def g2_rows():
    rows=[]
    for j,p,spec,c,word in [(1,1,{'operation':'omega','index':1},'1',[1]),
                           (2,1,{'operation':'g2quadratic'},'-3/4',[2,1]),
                           (3,2,{'operation':'omega','index':2},'1',[2]),
                           (4,2,{'operation':'g2product'},'1',[1,2])]:
        rows.append({'id':f'G2/P{p}/row{j}','kind':'G2','parabolic':p,'input':spec,
                     'published_terms':[{'coefficient':c,'word':word}]})
    return rows

def main():
    rows=g2_rows()+json.loads((ROOT/'part_iii_rows.json').read_text())
    results=[];calculators={};start=time.time()
    for row in rows:
        t=time.time();calc=calculators.setdefault(row['kind'],None)
        if calc is None:calc=calculators[row['kind']]=Calculator(row['kind'])
        result=calc.audit(row);result['seconds']=round(time.time()-t,3);results.append(result)
        print(row['id'],'PASS' if result['passed'] else 'MISMATCH',result['differences'],
              'basis',result['basis_classes_checked'],'seconds',result['seconds'],flush=True)
        (ROOT/'part_iii_results.json').write_text(json.dumps(results,indent=2)+'\n')
    print('COMPLETE',len(results),'rows',sum(r['printed_terms'] for r in results),'printed terms',
          'mismatches',sum(not r['passed'] for r in results),'seconds',round(time.time()-start,2),flush=True)

if __name__=='__main__':main()
