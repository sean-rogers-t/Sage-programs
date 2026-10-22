"""Classical equivariant characteristic maps and integral projective K-images.

Type A uses minimal parabolic cosets, reduced-subword cohomology and
0-Hecke K restrictions. Actual character integrality is checked.
The symplectic projective cokernel has an even/odd rank-n presentation.
"""
from functools import lru_cache
from itertools import combinations
from math import comb, factorial
from operator import index
from sage.all import ZZ, QQ, PolynomialRing, LaurentPolynomialRing, matrix, vector, prod, block_diagonal_matrix


def _integer(value, name, minimum=0):
    if isinstance(value, bool):
        raise ValueError(name + " must be an integer, not a boolean")
    try:
        result = index(value)
    except TypeError as error:
        raise ValueError(name + " must be an integer") from error
    if result < minimum:
        raise ValueError(name + " must be at least " + str(minimum))
    return result


def _swap(w, i):
    v = list(w)
    v[i], v[i + 1] = v[i + 1], v[i]
    return tuple(v)


def _length(w):
    return sum(a > b for i, a in enumerate(w) for b in w[i + 1:])


def _word(w):
    current = list(range(1, len(w) + 1))
    result = []
    for position, value in enumerate(w):
        j = current.index(value)
        while j > position:
            current[j - 1], current[j] = current[j], current[j - 1]
            result.append(j)
            j -= 1
    return tuple(result)


def _cosets(values, sizes):
    if len(sizes) == 1:
        yield tuple(values)
    else:
        for head in combinations(values, sizes[0]):
            rest = tuple(v for v in values if v not in head)
            for tail in _cosets(rest, sizes[1:]):
                yield head + tail


class EquivariantExpansion:
    """Schubert coordinates over one parent's exact coefficient ring."""
    def __init__(self, parent, coefficients):
        self.parent = parent
        self.coefficients = tuple(parent.coefficient_ring(c) for c in coefficients)
        if len(self.coefficients) != parent.schubert_rank:
            raise ValueError("wrong number of Schubert coordinates")

    def items(self):
        return tuple((w, c) for w, c in zip(self.parent.basis, self.coefficients) if c)

    def restrictions(self):
        values = self.parent.restriction_matrix() * vector(self.parent.coefficient_ring, self.coefficients)
        return dict(zip(self.parent.basis, values))

    def specialize(self):
        """Forget equivariance, keeping the Schubert basis."""
        return {w: self.parent._augment(c) for w, c in self.items() if self.parent._augment(c)}

    def __add__(self, other):
        if not isinstance(other, EquivariantExpansion) or other.parent is not self.parent:
            return NotImplemented
        return EquivariantExpansion(self.parent, [a+b for a,b in zip(self.coefficients, other.coefficients)])

    def __mul__(self, other):
        if isinstance(other, EquivariantExpansion):
            if other.parent is not self.parent:
                raise ValueError("expansions must have the same context")
            left, right = self.restrictions(), other.restrictions()
            return self.parent.from_restrictions({w: left[w]*right[w] for w in self.parent.basis})
        return EquivariantExpansion(self.parent, [other*a for a in self.coefficients])

    __rmul__ = __mul__

    def __eq__(self, other):
        return (isinstance(other,EquivariantExpansion) and other.parent is self.parent
                and self.coefficients==other.coefficients)

    def __pow__(self, exponent):
        exponent = _integer(exponent,"exponent",0)
        result = self.parent.schubert(tuple(range(1,self.parent.m+1)))
        base = self
        while exponent:
            if exponent%2:
                result = result*base
            exponent //= 2
            if exponent:
                base = base*base
        return result

    def __repr__(self):
        return "EquivariantExpansion({!r})".format(dict(self.items()))


class TypeAEquivariantMap:
    """Defining SL_m on any standard parabolic, with ZZ or QQ coefficients.

    Cartan variables x1,...,x_(m-1) are epsilon weights; x_m=-sum(x_i).
    The literal linear lattice is the root lattice. Over ZZ, expansion
    checks actual Laurent character integrality. Theory H means ordinary
    equivariant cohomology; K means equivariant algebraic K-theory.
    Signs agree with ordinary beta(omega_i)=epsilon_si. At projective
    fixed points the positive line class restricts to u_i or t_i.
    """
    def __init__(self, m, omitted_nodes=None, theory="K", coefficients=ZZ, max_fixed_points=720):
        self.m = _integer(m, "m", 2)
        self.rank = self.m-1
        if coefficients not in (ZZ, QQ):
            raise ValueError("coefficients must be ZZ or QQ")
        self.coefficients = coefficients
        self.theory = str(theory).upper()
        if self.theory not in ("H", "K"):
            raise ValueError("theory must be H or K")
        nodes = tuple(range(1, self.m))
        if omitted_nodes is None:
            omitted = nodes
        elif isinstance(omitted_nodes, (int, type(ZZ(1)))):
            omitted = (_integer(omitted_nodes, "omitted node", 1),)
        else:
            omitted = tuple(_integer(v, "omitted node", 1) for v in omitted_nodes)
        if len(set(omitted)) != len(omitted) or any(v >= self.m for v in omitted):
            raise ValueError("omitted nodes must be distinct nodes from 1 to m-1")
        self.omitted_nodes = tuple(sorted(omitted))
        self.retained_nodes = tuple(v for v in nodes if v not in omitted)
        cuts = (0,)+self.omitted_nodes+(self.m,)
        self.block_sizes = tuple(b-a for a,b in zip(cuts,cuts[1:]))
        self.schubert_rank = int(factorial(self.m)//prod(factorial(a) for a in self.block_sizes))
        if self.schubert_rank > _integer(max_fixed_points, "max_fixed_points", 1):
            raise ValueError("{} fixed points; raise max_fixed_points explicitly".format(self.schubert_rank))
        self.basis = tuple(sorted(_cosets(tuple(range(1,self.m+1)),self.block_sizes),
                                  key=lambda w: (_length(w),w)))
        self.lengths = tuple(_length(w) for w in self.basis)
        self.dimension = max(self.lengths)
        self.cartan_ring = PolynomialRing(QQ,self.rank,names=["x"+str(i) for i in nodes])
        self.weights = tuple(self.cartan_ring.gens())+(-sum(self.cartan_ring.gens()),)
        self.linear_source_generators = tuple(x-self.weights[-1] for x in self.weights[:-1])
        self.linear_lattice_matrix = matrix(ZZ,self.rank,lambda i,j: 1+int(i==j))
        self._rational_torus = LaurentPolynomialRing(QQ,self.rank,names=["t"+str(i) for i in nodes])
        diagonal = tuple(self._rational_torus.gens())
        self._torus_values = diagonal+(prod(diagonal)**-1,)
        mean = sum(self._torus_values)/self.m
        self._theta = tuple(t-mean for t in self._torus_values)
        constructor = PolynomialRing if self.theory=="H" else LaurentPolynomialRing
        prefix = "u" if self.theory=="H" else "t"
        self.coefficient_ring = constructor(coefficients,self.rank,names=[prefix+str(i) for i in nodes])
        self.parameters = tuple(self.coefficient_ring.gens())
        last = -sum(self.parameters) if self.theory=="H" else prod(self.parameters)**-1
        self.diagonal_parameters = self.parameters+(last,)
        self._restrictions = None

    def _polynomial(self, f):
        try:
            return self.cartan_ring(f)
        except (TypeError, ValueError) as error:
            raise ValueError("input must belong to this rational Cartan polynomial ring") from error

    def _invariant(self, f):
        for i in self.retained_nodes:
            images = list(self.weights)
            images[i-1],images[i] = images[i],images[i-1]
            if f(*images[:-1]) != f:
                raise ValueError("not invariant under retained Levi node "+str(i))

    def pullback(self, f):
        """Actual rational Laurent character Theta^*(f), before K."""
        return self._rational_torus(self._polynomial(f)(*self._theta[:-1]))

    def is_integral_source(self, f):
        f = self._polynomial(f)
        self._invariant(f)
        return all(c in ZZ for c in self.pullback(f).coefficients())

    def source_generators(self):
        """Full literal invariant generators on G/B or SL_m/P1.

        General parabolics have exact input membership/evaluation, but
        their full integral invariant-ring generators are not enumerated.
        """
        if not self.retained_nodes:
            return self.linear_source_generators
        if self.omitted_nodes == (1,):
            roots = [self.weights[0]-x for x in self.weights[1:]]
            return tuple(sum(prod(c) for c in combinations(roots,j)) for j in range(1,self.m))
        raise NotImplementedError("full integral Levi generators are not enumerated for this parabolic")

    def restrictions(self, f):
        f = self._polynomial(f)
        self._invariant(f)
        if self.coefficients is ZZ and not self.is_integral_source(f):
            raise ValueError("not a literal integral Springer character; use QQ explicitly for normalized inputs")
        R = self.coefficient_ring
        value = f if self.theory=="H" else self.pullback(f)
        return {w: R(value(*[self.diagonal_parameters[w[i]-1] for i in range(self.rank)]))
                for w in self.basis}

    def restriction_matrix(self):
        """Rows are fixed points, columns are opposite Schubert classes."""
        if self._restrictions is not None:
            return self._restrictions
        R = self.coefficient_ring
        identity = tuple(range(1,self.m+1))
        rows = []
        for w in self.basis:
            word = _word(w)
            prefix = identity
            factors = []
            for node in word:
                a,b = prefix[node-1]-1,prefix[node]-1
                factor = (self.diagonal_parameters[b]-self.diagonal_parameters[a] if self.theory=="H"
                          else 1-self.diagonal_parameters[a]/self.diagonal_parameters[b])
                factors.append(factor)
                prefix = _swap(prefix,node-1)

            @lru_cache(None)
            def at(v,j):
                if j==0:
                    return R.one() if v==identity else R.zero()
                if _length(v)>j:
                    return R.zero()
                node = word[j-1]-1
                skip = at(v,j-1)
                if v[node]<v[node+1]:
                    return skip
                q = factors[j-1]
                take = at(_swap(v,node),j-1)
                return skip+q*take if self.theory=="H" else (1-q)*skip+q*take

            rows.append([at(v,len(word)) for v in self.basis])
            at.cache_clear()
        self._restrictions = matrix(R,rows)
        self._restrictions.set_immutable()
        return self._restrictions

    def from_restrictions(self, values):
        if set(values)!=set(self.basis):
            raise ValueError("supply exactly one value at every fixed point")
        R = self.coefficient_ring
        field = R.fraction_field()
        M = self.restriction_matrix()
        coefficients = []
        for i,w in enumerate(self.basis):
            numerator = R(values[w])-sum(M[i,j]*coefficients[j] for j in range(i))
            try:
                quotient = field(numerator)/field(M[i,i])
                if self.theory=="K":
                    coefficients.append(R(quotient.numerator())*R(quotient.denominator()).inverse_of_unit())
                else:
                    coefficients.append(R(quotient))
            except (ArithmeticError,TypeError,ValueError) as error:
                raise ValueError("restrictions are not an integral/polynomial class in this context") from error
        return EquivariantExpansion(self,coefficients)

    def expand(self, f):
        return self.from_restrictions(self.restrictions(f))

    def schubert(self, w):
        w = tuple(w)
        if w not in self.basis:
            raise ValueError("not a minimal right-coset representative")
        return EquivariantExpansion(self,[int(v==w) for v in self.basis])

    def _augment(self, c):
        return self.coefficients(c(*([0 if self.theory=="H" else 1]*self.rank)))


def type_a_equivariant(m, omitted_nodes=None, theory="K", coefficients=ZZ, max_fixed_points=720):
    return TypeAEquivariantMap(m,omitted_nodes,theory,coefficients,max_fixed_points)


class SymplecticProjectiveKImage:
    """Integral Sp_(2n)/P1 K-image and its rank-n cokernel presentation.

    Spectral parameters work over ZZ[s1,...,sn]; base change
    s_i=t_i+t_i^-1 gives the default torus-character coefficient ring.
    Fitting outputs are determinantal ideals, not multivariable Smith forms.
    """
    def __init__(self, n, spectral=False):
        self.rank = _integer(n,"n",1)
        if not isinstance(spectral,bool):
            raise ValueError("spectral must be a boolean")
        self.spectral = spectral
        names = [("s" if spectral else "t")+str(i) for i in range(1,self.rank+1)]
        constructor = PolynomialRing if spectral else LaurentPolynomialRing
        self.coefficient_ring = constructor(ZZ,self.rank,names=names)
        self.parameters = tuple(self.coefficient_ring.gens())
        self.s = self.parameters if spectral else tuple(t+t**-1 for t in self.parameters)
        self.ring = PolynomialRing(self.coefficient_ring,"b")
        self.b = self.ring.gen()
        auxiliary = PolynomialRing(self.coefficient_ring,"s")
        s = auxiliary.gen()
        P = prod(s-value for value in self.s)
        self.half_ring = PolynomialRing(self.coefficient_ring,"u")
        u = self.half_ring.gen()
        self.half_even = self.half_ring(sum(P[2*j]*u**j for j in range(self.rank//2+1)))
        self.half_odd = self.half_ring(sum(P[2*j+1]*u**j for j in range((self.rank+1)//2)))
        self.even = self.ring(sum(P[2*j]*(self.b**2+4)**j for j in range(self.rank//2+1)))
        self.odd = self.ring(sum(P[2*j+1]*(self.b**2+4)**j for j in range((self.rank+1)//2)))
        self.cokernel_polynomial = self.even+self.b*self.odd
        self.image_relation = prod(self.b**2-(value**2-4) for value in self.s)
        self._matrix = None
        self._half_matrix = None
        self._half_hankel_matrix = None
        self._half_minors = {0:(self.coefficient_ring.one(),)}

    def half_presentation_matrix(self):
        """The floor(n/2) block 2H whose cokernel occurs twice.

        H multiplies O modulo E for even n, and E modulo O for odd n.
        The full R-module cokernel adds one R/(2) summand when n is odd.
        This splitting does not assert compatibility with the b action.
        """
        if self._half_matrix is None:
            r = self.rank//2
            modulus,multiplier = ((self.half_even,self.half_odd) if self.rank%2==0
                                  else (self.half_odd,self.half_even))
            u = self.half_ring.gen()
            self._half_matrix = matrix(self.coefficient_ring,r,r,
                lambda i,j:((2*multiplier*u**j)%modulus)[i])
            self._half_matrix.set_immutable()
        return self._half_matrix

    def block_presentation_matrix(self):
        """Unimodularly equivalent R-presentation: optional (2), then 2H,2H."""
        blocks = [self.half_presentation_matrix()]*2
        if self.rank%2:
            blocks.insert(0,matrix(self.coefficient_ring,1,1,[2]))
        answer = block_diagonal_matrix(blocks)
        answer.set_immutable()
        return answer

    def half_hankel_presentation_matrix(self):
        """Integral Hankel presentation equivalent to the half block 2H.

        Moments extract the u^(r-1) coefficient of 2q*u^k modulo monic p.
        The coefficient pairing is unimodular even in characteristic two;
        this uses no field trace or division by a discriminant.
        """
        if self._half_hankel_matrix is None:
            r = self.rank//2
            p,q = ((self.half_even,self.half_odd) if self.rank%2==0
                   else (self.half_odd,self.half_even))
            u = self.half_ring.gen()
            moments = [((2*q*u**j)%p)[r-1] for j in range(max(2*r-1,0))]
            self._half_hankel_matrix = matrix(self.coefficient_ring,r,r,lambda i,j:moments[i+j])
            self._half_hankel_matrix.set_immutable()
        return self._half_hankel_matrix

    def half_minor_generators(self, j, max_minors=10000):
        """Generate I_j(2H) using only the first j columns.

        Integral Schur cyclicity proves equality with the full minor ideal.
        The guard counts binomial(floor(n/2),j) determinants, rather than
        its square. These generators can still be redundant.
        """
        j = _integer(j,"minor size",0)
        r,R = self.rank//2,self.coefficient_ring
        if j>r:
            return (R.zero(),)
        count = comb(r,j)
        if count>_integer(max_minors,"max_minors",1):
            raise ValueError("{} initial-column minors required; raise max_minors explicitly".format(count))
        if j not in self._half_minors:
            M = self.half_presentation_matrix()
            values = set()
            for rows in combinations(range(r),j):
                value = M.matrix_from_rows_and_columns(rows,list(range(j))).det(algorithm="df")
                if value and -value not in values:
                    values.add(value)
            self._half_minors[j] = tuple(sorted(values)) if values else (R.zero(),)
        return self._half_minors[j]

    def _closed_small_block_fitting(self, i):
        """Explicit quadratic/cubic half-block formulas in ranks five to seven."""
        n,R = self.rank,self.coefficient_ring
        e = [R.one()]+[sum(prod(c) for c in combinations(self.s,k)) for k in range(1,n+1)]
        if n==5:
            A,B = e[1]*e[4]-e[5],e[1]*e[2]-e[3]
            J = (A,B)
            if i==1:
                F = prod(self.s[a]+self.s[b] for a in range(n) for b in range(a+1,n))
                return tuple(R(16*F*x) for x in J)
            if i==2:
                return (R(8*A*A),R(8*A*B),R(8*B*B))
            return tuple(R(4*x) for x in J)
        P,Q,R0 = e[2],e[4],e[6]
        a,b,c = ((e[1],e[3],e[5]) if n==6 else
                 (e[1]*e[2]-e[3],e[1]*e[4]-e[5],e[1]*e[6]-e[7]))
        J = (a,b,c)
        K = (c*c-Q*a*c+R0*a*b,b*c-P*a*c+R0*a*a,b*b-P*a*b-a*c+Q*a*a)
        def square(generators):
            return tuple(generators[a]*generators[b] for a in range(len(generators))
                         for b in range(a,len(generators)))
        if i==1:
            F = prod(self.s[a]+self.s[b] for a in range(n) for b in range(a+1,n))
            return tuple(R(2**(n-1)*F*x) for x in K)
        if i==2:
            values = square(K)
        elif i==3:
            values = tuple(x*y for x in J for y in K)
        elif i==4:
            values = square(J)
        else:
            values = J
        return tuple(R(2**(n-i)*x) for x in values)

    def presentation_matrix(self):
        """Multiplication by 2*O(b^2+4) modulo monic G=E+bO."""
        if self._matrix is None:
            columns = []
            for j in range(self.rank):
                f = (2*self.odd*self.b**j).quo_rem(self.cokernel_polynomial)[1]
                columns.append([f[i] for i in range(self.rank)])
            self._matrix = matrix(self.coefficient_ring,columns).transpose()
            self._matrix.set_immutable()
        return self._matrix

    def determinant_generator(self):
        return 2**self.rank*prod((self.s[i]+self.s[j])**2
                                for i in range(self.rank) for j in range(i+1,self.rank))

    def fitting_generators(self, i, max_minors=10000):
        i = _integer(i,"Fitting index",0)
        n,R = self.rank,self.coefficient_ring
        if i>=n:
            return (R.one(),)
        if i==0:
            return (R(self.determinant_generator()),)
        if n==3 and i==1:
            F = prod(self.s[a]+self.s[b] for a in range(n) for b in range(a+1,n))
            return (R(4*F),)
        if n==4 and i in (1,2):
            a = sum(self.s)
            d = sum(prod(c) for c in combinations(self.s,3))
            if i==2:
                return (R(4*a*a),R(4*a*d),R(4*d*d))
            F = prod(self.s[a]+self.s[b] for a in range(n) for b in range(a+1,n))
            return (R(8*F*a),R(8*F*d))
        if i==n-1:
            if n%2:
                return (R(2),)
            return tuple(2*sum(prod(c) for c in combinations(self.s,j))
                         for j in range(1,n,2))
        if n in (5,6,7):
            return self._closed_small_block_fitting(i)
        # Integral Hankel balancing reduces the entire convolution to
        # one product of equal or adjacent half-block minor ideals.
        size = n-i
        r = n//2
        degree = size-n%2
        a,b = degree//2,(degree+1)//2
        pairs = [(a,b)]
        sizes = {j for pair in pairs for j in pair}
        count = sum(comb(r,j) for j in sizes if j)
        count += comb(comb(r,a)+1,2) if a==b else comb(r,a)*comb(r,b)
        if count>_integer(max_minors,"max_minors",1):
            raise ValueError("{} half-block minor/product operations required; raise max_minors explicitly".format(count))
        ideals = {j:self.half_minor_generators(j,max_minors) for j in sizes}
        # M=2H: its products contribute 2^degree. In odd rank the
        # additional (2) block contributes the one remaining factor.
        scale = 2 if n%2 else 1
        answer,seen = [],set()
        for a,b in pairs:
            for index,left in enumerate(ideals[a]):
                for right in ideals[b][index if a==b else 0:]:
                    value = scale*left*right
                    if value and value not in seen and -value not in seen:
                        answer.append(value)
                        seen.add(value)
        return tuple(answer) if answer else (R.zero(),)

    def fitting_ideal(self, i, max_minors=10000):
        return self.coefficient_ring.ideal(self.fitting_generators(i,max_minors))

    def cyclic_cokernel_factors(self):
        """Proved global cyclic decompositions in ranks 1, 2 and 3 only."""
        R = self.coefficient_ring
        F = R(prod(self.s[i]+self.s[j] for i in range(self.rank)
                   for j in range(i+1,self.rank)))
        if self.rank==1:
            return (R(2),)
        if self.rank==2:
            return (2*F,2*F)
        if self.rank==3:
            return (R(2),2*F,2*F)
        raise NotImplementedError("global cyclic decompositions are proved here only for n<=3")

    def ordinary_smith_factors(self):
        # Augment before constructing matrices; no symbolic Laurent matrix
        # or minors are needed for this integral specialization.
        n = self.rank
        ring = PolynomialRing(ZZ,"u")
        u = ring.gen()
        even = ring(sum(comb(n,2*j)*(-2)**(n-2*j)*u**j for j in range(n//2+1)))
        odd = ring(sum(comb(n,2*j+1)*(-2)**(n-2*j-1)*u**j for j in range((n+1)//2)))
        modulus,multiplier = (even,odd) if n%2==0 else (odd,even)
        r = n//2
        M = matrix(ZZ,r,r,lambda i,j:((2*multiplier*u**j)%modulus)[i])
        paired = tuple(abs(int(s)) for s in M.elementary_divisors() for _ in range(2))
        return (1,)*n+((2,) if n%2 else ())+paired


def symplectic_projective_k_image(n, spectral=False):
    return SymplecticProjectiveKImage(n,spectral)
