# An integral equivariant image whose specialization has a kernel

Status, 20 October 2026: explicit repository counterexample with a written
integral proof and independently reconstructed Sage/SymPy integer lattices.
It resolves the first Track E1 pilot negatively for **automatic faithful
specialization**. It does not classify every specialization defect or
assert novelty. See the [proof/source audit](LOCAL_IMAGE_AUDIT.md).

## 1. The precise Springer–Borel source

Take G=SL_4(C), rho its defining representation, P=P_1, so G/P=P^3.
Let X=X(T) be the SL_4 weight lattice and Q its root lattice.
The [ordinary-family proof](LOCAL_IMAGE_BOUND.md#6-a-complete-projective-family-for-the-defining-sl_n-representation)
identifies the full literal source as

\[
\mathcal C_\rho=\operatorname{Sym}_{\mathbb Z}(Q),\qquad
A_{\rho,L}=\mathbb Z[E_1,E_2,E_3].
\]

Here E_j=e_j(x_1-x_2,x_1-x_3,x_1-x_4), and the Levi Weyl group is S_3.
These root differences form an integral basis of Q, and S_3 permutes them, so
the displayed invariant ring is the full ring, not a selected subalgebra.
This uses the defining representation of the simply connected group;
there is no change to an adjoint representation.

## 2. Ambient equivariant ring and coefficient conventions

Write

\[
S=\mathbb Z[t_1,t_2,t_3],\qquad t_4=-t_1-t_2-t_3,\qquad
a=e_2(t),\quad b=e_3(t),\quad c=e_4(t).
\]

These are integral weight coordinates for the torus: the first three
diagonal characters form a basis of X, and the fourth is their negative
sum. Choose the matching characteristic-class signs with beta(x_1)=h
and the projective bundle presentation

\[
H=H_T^{2*}(\mathbb P^3,\mathbb Z)
 =S[h]/\bigl(\prod_i(h-t_i)\bigr)
 =S[h]/(h^4+a h^2-bh+c).
\tag{1}
\]

Changing simultaneously the line-bundle and equivariant-parameter signs
does not change the arithmetic conclusions. Every t_i and h has
polynomial degree one; a,b,c have degrees two, three and four.
The parameter ideal is I=(t_1,t_2,t_3). The ambient ring H is free
over S on 1,h,h^2,h^3, and H/IH=H^(2*)(P^3,Z).

Formula (1) is the standard equivariant projective bundle relation.
The [equivariant oriented characteristic-map framework](https://arxiv.org/abs/1409.7111)
provides general background; it is not being claimed as new.

The shift identity from the ordinary-family proof now retains the
global equivariant symmetric functions. Consequently

\[
f_1=\beta_T(E_1)=4h,\quad
f_2=\beta_T(E_2)=6h^2+a,\quad
f_3=\beta_T(E_3)=4h^3+2ah-b.
\tag{2}
\]

At the fixed point indexed by i, h specializes to t_i and (2) equals
e_j(t_i-t_k : k ne i). Thus the formulas also follow independently by
their fixed-point restrictions.

Define the equivariant image, **including the full coefficient ring**, by

\[
J_T=S[f_1,f_2,f_3]\subset H.
\]

The ordinary image is

\[
J=\mathbb Z\oplus4\mathbb Z h\oplus2\mathbb Z h^2
                     \oplus4\mathbb Z h^3
\subset\mathbb Z[h]/(h^4).
\]

Specialization induces a surjection J_T/IJ_T -> J. The question is
whether this map is injective; the ambient freeness of H does not answer it.

## 3. Convenient exact generators and low-degree products

Put

\[
u=4h,\qquad v=2h^2,\qquad w=4h^3+2ah.
\]

Then J_T=S[u,v,w]. Indeed
v=3f_2-f_1^2-3a, f_2=3v+a and f_3=w-b.
Let M be the S-submodule with basis 1,u,v,w. These vectors are linearly
independent over S, by their triangular h-coefficients.

Reducing with (1) gives

\[
\begin{aligned}
u^2&=8v,& uv&=2w-au,\\
uw&=-4av+4bu-16c,&v^2&=-2av+bu-4c.
\end{aligned}
\tag{3}
\]

Every source monomial of degree at most four reduces into M by (3).
Conversely M is contained in J_T. Therefore (J_T)_d=M_d for d<=4.
The same calculation shows these specializations agree with J through
degree four: free rank one in degrees zero through three and zero in
degree four, with no extra torsion.

## 4. A nonzero specialization class of order two

Define the degree-five class

\[
\nu=2a^2h.
\]

The next mixed product gives the two exact identities

\[
\boxed{\nu=vw+aw-4bv+2cu\in J_T,\qquad
                    2\nu=a^2u\in IJ_T.}
\tag{4}
\]

On the other hand nu does not belong to M. Its h^3 coefficient forces
the coefficient of w in a putative M-expression to be zero; its h^2
coefficient then forces the coefficient of v to be zero. The h coefficient
would require 2a^2 to be divisible by 4 in S. This is false:

\[
a=-(t_1^2+t_2^2+t_3^2+t_1t_2+t_1t_3+t_2t_3)
\]

and a^2 has coefficient one at t_1^4.

Every degree-five element of IJ_T is a sum of parameters times elements
of J_T of degrees at most four. Hence (IJ_T)_5 is contained in M_5.
Conversely every coefficient of an M_5-expression has positive parameter
degree, since the four basis degrees are 0,1,2,3. Thus

\[
(IJ_T)_5=M_5,\qquad \nu\notin(IJ_T)_5.
\]

All degree-five source monomials except vw reduce into M by (3), and
(4) reduces vw to nu modulo M. Therefore

\[
(J_T)_5=M_5+\mathbb Z\nu,\qquad
\boxed{(J_T/IJ_T)_5\cong\mathbb Z/2.}
\tag{5}
\]

But J_5=0, since dim(P^3)=3. The first failure of injectivity is
therefore exactly a Z/2 in polynomial degree five, cohomological degree
ten. This statement describes the first defect; it does not assert that
there are no further defects in higher equivariant degrees.

## 5. Meaning and remaining questions

This is torsion in the **specialized image module** J_T/IJ_T.
The ambient cohomology of P^3 is torsion-free, and the ordinary image J
is torsion-free as a subgroup. The parameter ideal times J_T must not
be confused with J_T intersect IH; their difference is exactly the
specialization kernel. Ordinary specialization of the ambient ring
still works normally.

In homological terms, for C=H/J_T the kernel is Tor_1^S(C,S/I),
because H is free over S. The example witnesses a failure of flatness
of the ambient quotient C with respect to this specialization.
It is not evidence that the original Springer map or Borel map is incorrect.

The first Track E1 milestone is thus a precise negative integral result,
alongside the [general good-prime comparison](LOCAL_IMAGE_BOUND.md#6a-equivariant-good-primes-and-the-support-of-specialization-defects).
For every full literal source in that theorem's scope, each fixed-degree
kernel is finite with primes dividing t_G e_rho; it vanishes after
localization away from those primes. Here t_G=1 and e_rho=4.
Inverting two makes h=u/4 available, so the whole equivariant ambient
ring is already recovered. The next question is a sharp criterion or
description at the allowed bad primes.

The exact [Sage calculation](experiments/local_image_equivariant_pilot.py)
checks (2) at twelve fixed-point/generator pairs, all displayed product
identities, and integer lattices in degrees zero through five.
A separate [SymPy calculation](experiments/verify_local_image_independent.py)
imports neither Sage nor sage_borel, rebuilds both degree-five lattices
in the rank-52 ambient group, and obtains 51 unit factors and a single
factor 2 for (IJ_T)_5 inside (J_T)_5.
See the [Sage report](experiments/local-image-equivariant-pilot.json)
and [independent report](experiments/local-image-independent.json).

The ambient characteristic-map/projective-bundle results and the general
base-change/Tor mechanism are established mathematics. Targeted literature
searches did not identify a prior statement of this exact literal-source
example; this limited search does not establish novelty or priority.
