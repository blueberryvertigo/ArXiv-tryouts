# An integer polynomial Hurwitz counterexample in dimension six

**Run the certificate** (standard library only, no dependencies):

```bash
python3 weak_markus_yamabe_six_independent.py
```

It verifies `det(lambda I - DX) = (lambda + 1)^6` exactly in Z[x,y,z,lambda] and
the three equilibria with exact fractions, then writes a JSON report beside the
script.

**Date:** 6 September 2026.

**Scope.** The result below is an explicit algebraic counterexample to the
six-dimensional polynomial weak Markus–Yamabe statement. It is not a result
about the Riemann Hypothesis. The identities have been checked with two
independent exact-arithmetic implementations. Publication priority and
independent expert review have not been established.

## 1. Statement and construction

A real differentiable map is called Hurwitz if every eigenvalue of its
Jacobian has strictly negative real part at every point. The weak
Markus–Yamabe statement asserts that such a map is injective.

**Theorem.** There is an integer-coefficient polynomial vector field
\(X:\mathbb R^6\to\mathbb R^6\), of total degree ten, such that

\[
\det(\lambda I_6-DX(Z))=(\lambda+1)^6
\quad\text{for every }Z\in\mathbb R^6,
\]

and \(X\) has exactly three distinct real zeros. Two of the zeros have
integer coordinates. Consequently the polynomial weak Markus–Yamabe
statement fails in dimension six and in every higher dimension.

For coordinates \(Z=(x,y,z,a,b,c)\), put

\[
q=y-6xz
\]

and define

\[
\begin{aligned}
h_1&=x^2(3q+2xz),\\
h_2&=24xq(q+xz),\\
h_3&=2q(2q+3xz),\\
h_4&=12x^2q^2(3q+2xz),\\
h_5&=2xq^2(7q+6xz),\\
h_6&=4x^2q^3(3q+2xz).
\end{aligned}
\]

The vector field is

\[
X(Z)=
\begin{pmatrix}
-x+h_1\\
1-y+a-h_2\\
-z+b-h_3\\
-a-h_4\\
-b+c-h_5\\
-c-h_6
\end{pmatrix}.
\tag{1}
\]

Its three zeros are

\[
\begin{aligned}
Z_0&=(0,1,-4,0,0,0),\\
Z_1&=(-1,-11,2,12,-6,4),\\
Z_2&=(1,7,5/4,-3,-3/2,1/2).
\end{aligned}
\tag{2}
\]

## 2. The associated three-variable map has constant determinant one

Define

\[
G(x,y,z)=
\begin{pmatrix}
x-h_1\\
y+h_2+h_4\\
z+h_3+h_5+h_6
\end{pmatrix}.
\tag{3}
\]

The starting ingredient is the three-variable map attributed to Levent
Alpöge in 2026. With \(u=1+sr\), write it as

\[
\mathcal K(s,r,w)=
\begin{pmatrix}
u^3w+r^2u(4+3sr)\\
r+3su^2w+3sr^2(4+3sr)\\
2s-3s^2r-s^3w
\end{pmatrix}.
\tag{4}
\]

Its determinant is

\[
\det D\mathcal K=-2.
\tag{5}
\]

For an elementary coefficient-wise verification, write
\(\mathcal K=A_0(s,r)+wB_0(s,r)\). Direct differentiation gives

\[
\begin{aligned}
\det(A_{0,s},A_{0,r},B_0)&=-2,\\
\det(B_{0,s},A_{0,r},B_0)&=-6s^2(1+sr),\\
\det(A_{0,s},B_{0,r},B_0)&=6s^2(1+sr),\\
\det(B_{0,s},B_{0,r},B_0)&=0.
\end{aligned}
\]

Expanding the determinant in powers of \(w\) proves (5).

Expansion of (3) gives the polynomial identity

\[
G(x,y,z)=
\begin{pmatrix}
\mathcal K_3(2x,y-6xz,z)/4\\
\mathcal K_2(2x,y-6xz,z)\\
\mathcal K_1(2x,y-6xz,z)
\end{pmatrix}.
\tag{6}
\]

The input map \((x,y,z)\mapsto(2x,y-6xz,z)\) has determinant two.
The output linear map in (6) has determinant \(-1/4\). Thus

\[
\det DG=(-1/4)(-2)(2)=1.
\tag{7}
\]

## 3. The six-dimensional Jacobian has determinant one

The last three equations \(X_4=X_5=X_6=0\) eliminate the auxiliary
variables uniquely:

\[
a=-h_4,\qquad c=-h_6,\qquad b=-h_5-h_6.
\tag{8}
\]

Substituting into the first three equations gives

\[
(0,1,0)^{\mathsf T}-G(x,y,z).
\tag{9}
\]

For the determinant, partition \(DX\) according to the base variables
\((x,y,z)\) and auxiliary variables \((a,b,c)\). The auxiliary diagonal
block is

\[
D=\begin{pmatrix}-1&0&0\\0&-1&1\\0&0&-1\end{pmatrix},
\qquad \det D=-1.
\]

The Schur complement is precisely \(-DG\), by the elimination (8).
Therefore

\[
\det DX=\det D\,\det(-DG)=(-1)(-1)=1.
\tag{10}
\]

## 4. Weighted similarity forces the entire spectrum

Let

\[
M=I_6+DX.
\]

All entries of \(M\) depend only on \((x,y,z)\), not on \((a,b,c)\).
Give \(x,y,z\) weights \(0,1,1\). Then \(q\) has weight one, and
\(h_1,\ldots,h_6\) have respective weights

\[
1,2,2,3,3,4.
\]

Set

\[
\nu=(0,1,1,2,2,3),
\qquad D_t=\operatorname{diag}(1,t,t,t^2,t^2,t^3).
\]

Each nonzero entry satisfies

\[
M_{ij}(x,ty,tz)=t^{1+\nu_i-\nu_j}M_{ij}(x,y,z).
\]

Indeed, differentiation with respect to \(x\) does not decrease weighted
degree, whereas differentiation with respect to \(y\) or \(z\) decreases
it by one. The three constant entries connecting the auxiliary chains
also satisfy the displayed identity. Consequently, for every \(t\ne0\),

\[
M(x,ty,tz)=tD_tM(x,y,z)D_t^{-1}.
\tag{11}
\]

By (10),

\[
\det(I_6-M)=\det(-DX)=1.
\]

Evaluate this at \((x,ty,tz)\), and use (11). Similarity invariance yields

\[
\det(I_6-tM(x,y,z))=1
\qquad(t\ne0).
\tag{12}
\]

Both sides are polynomials in \(t\), so (12) is an identity for every
\(t\). It follows that

\[
\det(\lambda I_6-M)=\lambda^6.
\]

Since \(DX=M-I_6\),

\[
\boxed{\det(\lambda I_6-DX)=(\lambda+1)^6.}
\tag{13}
\]

This argument never divides by \(x\). In particular, the zero weight
assigned to \(x\) creates no exception at \(x=0\).

## 5. The collisions and the complete zero set

At \(Z_0\), the auxiliary value is \(q=1\), and

\[
(h_1,h_2,h_3,h_4,h_5,h_6)=(0,0,4,0,0,0).
\]

At \(Z_1\), again \(q=1\), and

\[
(h_1,h_2,h_3,h_4,h_5,h_6)=(-1,24,-8,-12,10,-4).
\]

Substitution into (1) proves \(X(Z_0)=X(Z_1)=0\), already enough to
establish non-injectivity.

For completeness, (8) and (9) identify the entire zero set with
\(G^{-1}(0,1,0)\). By (6), this is equivalent to

\[
\mathcal K(s,r,w)=(0,1,0).
\]

When \(s=0\), the unique solution is \((s,r,w)=(0,1,-4)\).
When \(s\ne0\), set \(R=sr\), \(S=s^2w\). The equation
\(\mathcal K_3=0\) gives \(S=2-3R\). On this branch, direct substitution
gives

\[
s\mathcal K_2=6+4R,
\qquad
s^2\mathcal K_1=(R+1)(R+2).
\]

Thus \(R\) is either \(-1\) or \(-2\), while \(s=6+4R\). The two
remaining solutions are

\[
(s,r,w)=(2,-1/2,5/4),\qquad(-2,1,2).
\]

Converting through (6) and then (8) gives exactly \(Z_2\) and \(Z_1\)
in (2). This proves the theorem.

## 6. Consequences, verification, and provenance

For \(n>6\), append \(n-6\) components \(-u_j\). The spectrum remains
entirely \(-1\), and all three zeros survive with \(u=0\). If the
formulation requires a zero at the origin, replace \(X(Z)\) by
\(X(Z+Z_0)\); this preserves integer coefficients and the Jacobian
identity.

The 5 August 2026 paper by Castañeda, Honorato and Valenzuela-Henríquez
constructs a counterexample in dimension fourteen and explicitly leaves
dimensions three through thirteen unsettled. It supplies the ordinary
chain-realization idea and a graded-similarity lemma. The construction
here instead uses the weighted grading with a zero-weight base variable,
after a triangular change of coordinates. Only three auxiliary variables
are required.

This comparison is based on the sources located on 6 September 2026. It
is not a certification of priority or independent expert acceptance.

Two executable certificates accompany this note:

* `weak_markus_yamabe_six_verify.py`: exact SymPy verification of the full
  characteristic polynomial, determinant, weighted identity, Schur
  complement, starting-map identity, and equilibria.
* `weak_markus_yamabe_six_independent.py`: a separate implementation using
  only Python's standard library. Sparse integer polynomial arithmetic and
  recursive determinant expansion verify (13) in
  \(\mathbb Z[x,y,z,\lambda]\); exact fractions verify (2).

Both executed successfully. The independent characteristic-polynomial
residual had zero nonzero monomials. These are exact algebraic checks,
not numerical interpolation or testing at finitely many sample points.
They are not proof-assistant formalizations. The symbolic argument in
Sections 2–5 supplies the proof independently of the programs.

### References

1. Álvaro Castañeda, Gerardo Honorato, Francisco Valenzuela-Henríquez,
   *The weak Markus–Yamabe conjecture fails in dimension 14*,
   arXiv:2608.05392, 5 August 2026.
2. *Formal Verification of an Explicit Counterexample to the Jacobian
   Conjecture*, Archive of Formal Proofs, entry
   `Jacobian_Counterexample`, July 2026. This verifies the three-variable
   starting map attributed to Levent Alpöge; it does not verify the new
   six-variable construction here.
3. Guowu Meng and Liang Yang, *A five-variable counterexample to the
   Hessian conjecture, and the low-dimensional status of the Jacobian
   and Hessian conjectures*, arXiv:2607.22198. Formula (1) records the
   same three-variable starting map.
