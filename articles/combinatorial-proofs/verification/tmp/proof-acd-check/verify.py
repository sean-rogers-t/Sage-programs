"""Exact audit of the proposed C proof and the D Chern-polynomial analogue.

Run with SageMath. This imports only the new root-data BGG Calculator, not any
historical scripts. Every cohomological equality is compared over QQ in the
complete coset basis; separate one-row word and formal-series checks are
recorded explicitly as such.
"""
from pathlib import Path
import sys, json, time
from itertools import combinations
from sage.all import QQ, prod

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
sys.path.insert(0, str(BASE / "papers" / "table_audit"))
sys.path.insert(0, str(BASE / "tmp" / "dissertation-proof-audit"))
from recompute import Calculator
from check_theorem55 import partition

def elementary(values, degree, R):
    return R(sum(prod(subset) for subset in combinations(values, degree)))

def complete(values, bound, R):
    h = [R(1)] + [R(0)] * bound
    for value in values:
        for degree in range(1, bound + 1):
            h[degree] += value * h[degree - 1]
    return h

def coordinates(calc):
    w = calc.x
    n = calc.rank
    if calc.kind[0] == "C":
        return [w[0]] + [w[i] - w[i-1] for i in range(1, n)]
    if calc.kind[0] == "D":
        return ([w[0]] + [w[i] - w[i-1] for i in range(1, n-2)]
                + [w[-2] + w[-1] - w[-3], w[-1] - w[-2]])
    raise ValueError(calc.kind)

def terms(expansion):
    return [{"word": [int(i) for i in w.reduced_word()], "coefficient": str(c)}
            for w, c in sorted(expansion.items(), key=lambda p: p[0].reduced_word())]

def main():
    started = time.time()
    checks = []
    one_row_checks = []
    gf_checks = []
    for n in range(2, 9):
        calc = Calculator("C" + str(n))
        x = coordinates(calc)
        bound = 6
        H = {m: complete(x[:m], bound, calc.R) for m in range(1, n+1)}
        F = {}
        for m in range(n):
            h = complete(x[:m], bound, calc.R)
            hs = complete([v*v for v in x[m:]], bound // 2, calc.R)
            F[m] = [sum(h[d-2*j]*hs[j] for j in range(d//2+1))
                    for d in range(bound+1)]
        # Coefficientwise verification of each claimed formal-series identity.
        for m in range(1, n+1):
            active = m
            target = H[m+1] if m < n else F[n-1]
            for d in range(1, bound+1):
                assert calc.divided(H[m][d], active) == target[d-1]
                for j in range(1, n+1):
                    if j != active:
                        assert calc.divided(H[m][d], j) == 0
                gf_checks.append({"kind":"C", "n":n, "series":"H", "m":m,
                                  "degree":d, "active_node":active, "passed":True})
        for m in range(1, n):
            for d in range(1, bound+1):
                assert calc.divided(F[m][d], m) == F[m-1][d-1]
                for j in range(1, n+1):
                    if j != m:
                        assert calc.divided(F[m][d], j) == 0
                gf_checks.append({"kind":"C", "n":n, "series":"F", "m":m,
                                  "degree":d, "active_node":m, "passed":True})
        for r in range(1, n+1):
            k = n-r
            h = complete(x[:r], max(2*k+1, n+k if n <= 5 else 0), calc.R)
            E = [elementary(x[:r], i, calc.R) for i in range(r+1)]
            for i in range(1, r+1):
                word = list(range(r-i+1, r+1))
                w = calc.W.from_reduced_word(word)
                out = calc.expansion(E[i], r)
                assert out == {w:QQ(1)}
                lam = partition(calc, w, r)
                assert lam == ((i,) if k == 0 else (1,)*i)
                checks.append({"kind":"C", "n":n, "r":r, "k":k, "i":i,
                               "factor":"GL", "passed":True,
                               "basis_classes_checked":len(calc.representatives(r,i)),
                               "word":word, "partition":list(lam)})
            path = list(range(r, n+1)) + list(range(n-1, 0, -1))
            for d in range(1, n+k+1):
                word = list(reversed(path[:d]))
                w = calc.W.from_reduced_word(word)
                assert w.length() == d
                assert w in calc.representatives(r,d)
                full_support_verified = n <= 5
                if full_support_verified:
                    out = calc.expansion(h[d], r)
                    assert out == {w:QQ(1)}
                assert partition(calc,w,r) == (d,)
                one_row_checks.append({"kind":"C", "n":n, "r":r, "k":k,
                                       "degree":d, "word":word, "passed":True,
                                       "full_support_verified":full_support_verified,
                                       "basis_classes_checked":len(calc.representatives(r,d)) if full_support_verified else 0})
            for i in range(1,k+1):
                f = elementary([v*v for v in x[r:]],i,calc.R)
                out = calc.expansion(f,r)
                chern = h[i]**2 + 2*sum((-1)**j*h[i+j]*h[i-j] for j in range(1,i+1))
                assert calc.expansion(chern,r) == out
                symmetric_identity = sum((-1)**q*(h[2*i-q]*h[q]-h[2*i-q+1]*(h[q-1] if q else 0)) for q in range(i+1))
                assert symmetric_identity == complete([v*v for v in x[:r]],i,calc.R)[i]
                got = {partition(calc,w,r):c for w,c in out.items()}
                expected = {tuple(v for v in (2*i-q,q) if v):QQ((-1)**(i+q))
                            for q in range(i+1) if r >= 2 or q == 0}
                assert got == expected
                checks.append({"kind":"C", "n":n, "r":r, "k":k, "i":i,
                               "factor":"squared_tail", "passed":True,
                               "basis_classes_checked":len(calc.representatives(r,2*i)),
                               "computed_terms":terms(out),
                               "partitions":{str(lam):str(c) for lam,c in got.items()}})
        print("C%d complete" % n, flush=True)
    for n in range(4,9):
        calc = Calculator("D"+str(n))
        x = coordinates(calc)
        for r in range(1,n-1):
            k = n-r
            h = complete(x[:r],2*k+1,calc.R)
            for i in range(1,k+1):
                f = elementary([v*v for v in x[r:]],i,calc.R)
                assert all(calc.reflection(f,j)==f for j in range(1,n+1) if j != r)
                out = calc.expansion(f,r)
                chern = h[i]**2 + 2*sum((-1)**j*h[i+j]*h[i-j] for j in range(1,i+1))
                assert calc.expansion(chern,r) == out
                rhs = sum((-1)**(i+q)*(h[2*i-q]*h[q]-h[2*i-q+1]*(h[q-1] if q else 0)) for q in range(i+1))
                assert rhs == chern
                # These are pullbacks of ordinary Schur polynomials, without
                # claiming that each is one type-D Schubert basis element.
                schur_out = []
                for q in range(i+1):
                    schur = h[2*i-q]*h[q]-h[2*i-q+1]*(h[q-1] if q else 0)
                    schur_out.append({"partition":[2*i-q,q],
                                      "terms":terms(calc.expansion(schur,r)) if schur else []})
                checks.append({"kind":"D", "n":n, "r":r, "k":k, "i":i,
                               "factor":"squared_tail", "passed":True,
                               "basis_classes_checked":len(calc.representatives(r,2*i)),
                               "computed_terms":terms(out),
                               "ordinary_schur_pullbacks":schur_out})
        print("D%d complete" % n, flush=True)
    data = {"method":"Exact root-data BGG, complete minimal-right-coset basis",
            "checks":checks, "one_row_checks":one_row_checks,
            "generating_function_checks":gf_checks,
            "counts":{"main":len(checks), "C_main":sum(r["kind"]=="C" for r in checks),
                      "D_main":sum(r["kind"]=="D" for r in checks),
                      "C_one_row":len(one_row_checks), "C_generating_coefficients":len(gf_checks),
                      "C_one_row_full_basis":sum(r["full_support_verified"] for r in one_row_checks),
                      "main_basis_coefficients":sum(r["basis_classes_checked"] for r in checks),
                      "one_row_basis_coefficients":sum(r["basis_classes_checked"] for r in one_row_checks)},
            "seconds":round(time.time()-started,3)}
    (HERE/"verification.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(data["counts"]), "seconds", data["seconds"], flush=True)

if __name__ == "__main__":
    main()
