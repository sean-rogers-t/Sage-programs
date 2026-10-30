"""Research prototype for normalized character indecomposables.

Implements the classical differentiated tensor-product recursion (Davis,
Section 2), then applies it to the K-Springer comparison data. It does
not claim to invent a representation-ring algorithm or compute new
ambient K-theory rings. Sage root numbering is recorded explicitly.
"""
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from sage.all import QQ, WeylCharacterRing, matrix, vector

ROOT = Path(__file__).resolve().parents[1]


class CotangentPilot:
    def __init__(self, kind):
        self.kind=kind;self.W=WeylCharacterRing(kind,style='coroots')
        self.space=self.W.space();self.r=self.space.root_system.cartan_type().rank()
        self.coroots=self.space.simple_coroots();self.weights=self.space.fundamental_weights()
        roots=self.space.roots();self.root_count=len(roots)
        self.scale=QQ(2)/max(a.inner_product(a) for a in roots)
        self.units=[tuple(int(i==j) for i in range(self.r)) for j in range(self.r)]
        self.rho=sum(self.weights.values())
        self.fund_indices=vector(QQ,[self.index(e) for e in self.units])

    def label(self, weight):
        return tuple(int(weight.inner_product(self.coroots[i+1])) for i in range(self.r))

    @lru_cache(None)
    def dimension(self, label):
        return int(self.W(*label).degree())

    @lru_cache(None)
    def index(self, label):
        w=sum((label[i]*self.weights[i+1] for i in range(self.r)),self.space.zero())
        if self.kind.startswith('A'):
            coordinates=w.to_vector();average=sum(coordinates)/(self.r+1)
            w=self.space.from_vector(vector(QQ,[x-average for x in coordinates]))
        return QQ(self.dimension(label))*self.scale*w.inner_product(w+2*self.rho)/(self.r+self.root_count)

    @lru_cache(None)
    def derivative(self, label):
        if not any(label):return vector(QQ,self.r)
        if label in self.units:return vector(QQ,label)
        i=next(i for i,c in enumerate(label) if c)
        mu=tuple(c-int(j==i) for j,c in enumerate(label))
        decomposition={self.label(w):int(c) for w,c in (self.W(*self.units[i])*self.W(*mu)).monomial_coefficients().items()}
        assert decomposition.pop(label)==1
        answer=self.dimension(mu)*vector(QQ,self.units[i])+self.dimension(self.units[i])*self.derivative(mu)
        for nu,c in decomposition.items():answer-=c*self.derivative(nu)
        assert all(a.denominator()==1 for a in answer)
        assert self.fund_indices.dot_product(answer)==self.index(label)
        answer.set_immutable()
        return answer

    def report(self,label):
        derivative=self.derivative(label);index=self.index(label)
        return dict(highest_weight=list(label),dimension=self.dimension(label),
                    trace_index=str(index),cotangent_coordinates=list(map(int,derivative)),
                    normalized_coordinates=list(map(str,derivative/index)))


def main():
    start=perf_counter();reports=[]
    cases=[('A2',[(2,0),(0,2),(1,1),(3,0)]),
           ('G2',[(2,0),(1,1),(3,0)]),
           ('F4',[(0,0,0,2)]),
           ('E6',[(2,0,0,0,0,0),(1,0,0,0,0,1)]),
           ('E7',[(0,0,0,0,0,0,2)]),
           ('E8',[(0,0,0,0,0,0,0,2)])]
    for kind,labels in cases:
        case_start=perf_counter();pilot=CotangentPilot(kind)
        rows=[]
        for label in labels:
            row=pilot.report(label);rows.append(row);print(kind,row,flush=True)
        # Differentiating a tensor square independently from the direct
        # recursion checks tau_(V tensor V)=tau_V for the selected fundamental.
        node={'A2':0,'G2':0,'F4':3,'E6':0,'E7':6,'E8':7}[kind]
        e=pilot.units[node];decomp=(pilot.W(*e)**2).monomial_coefficients()
        derivative=sum((int(c)*pilot.derivative(pilot.label(w)) for w,c in decomp.items()),vector(QQ,pilot.r))
        idx=sum(int(c)*pilot.index(pilot.label(w)) for w,c in decomp.items())
        assert derivative/idx==vector(QQ,e)/pilot.index(e)
        reports.append(dict(cartan_type=kind,fundamental_dimensions=[pilot.dimension(e) for e in pilot.units],
                            fundamental_indices=list(map(str,pilot.fund_indices)),irreducibles=rows,
                            tensor_square_normalized_cotangent_unchanged=True,
                            seconds=round(perf_counter()-case_start,3)))
    report=dict(result='passed',date='2026-10-22',cases=reports,
                seconds=round(perf_counter()-start,3),
                scope='Finite differentiated tensor-product checks of known indecomposable representation data used by the new image comparison; not an exceptional K-image classification.',
                source_sha256={'experiments/other_theories_cotangent_pilot.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    (ROOT/'experiments/other-theories-cotangent-pilot.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
