#!/usr/bin/env python3
"""Independent exact certificate using only the Python standard library.

A minimal sparse polynomial implementation verifies
    det(lambda I_6 - DX) = (lambda+1)^6
in Z[x,y,z,lambda]. Separate Fraction evaluations verify the three zeros.
No SymPy, numerical interpolation, random sampling, or floating point is used.
"""
from __future__ import annotations
from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path
from typing import Dict, Tuple

Monomial = Tuple[int, int, int, int]
ZERO_MONOMIAL: Monomial = (0, 0, 0, 0)


class Poly:
    """Sparse polynomials in x,y,z,lambda with integer coefficients."""
    def __init__(self, value=0):
        if isinstance(value, Poly):
            self.terms = dict(value.terms)
        elif isinstance(value, int):
            self.terms = {} if value == 0 else {ZERO_MONOMIAL: value}
        elif isinstance(value, dict):
            self.terms: Dict[Monomial, int] = {k: v for k, v in value.items() if v}
        else:
            raise TypeError(type(value))

    @staticmethod
    def variable(index: int) -> 'Poly':
        exponent = [0, 0, 0, 0]
        exponent[index] = 1
        return Poly({tuple(exponent): 1})

    def __add__(self, other):
        other = Poly(other)
        result = dict(self.terms)
        for key, value in other.terms.items():
            result[key] = result.get(key, 0) + value
            if result[key] == 0:
                del result[key]
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + (-Poly(other))

    def __rsub__(self, other):
        return Poly(other) - self

    def __mul__(self, other):
        other = Poly(other)
        result = {}
        for ka, va in self.terms.items():
            for kb, vb in other.terms.items():
                key = tuple(a+b for a, b in zip(ka, kb))
                result[key] = result.get(key, 0) + va*vb
        return Poly(result)

    __rmul__ = __mul__

    def __pow__(self, exponent: int):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("A nonnegative integer exponent is required")
        result, base = Poly(1), self
        while exponent:
            if exponent & 1:
                result = result*base
            base = base*base
            exponent >>= 1
        return result

    def derivative(self, index: int):
        result = {}
        for key, value in self.terms.items():
            degree = key[index]
            if degree:
                new_key = list(key)
                new_key[index] -= 1
                result[tuple(new_key)] = degree*value
        return Poly(result)

    def __bool__(self):
        return bool(self.terms)


def main():
    x, y, z, lam = [Poly.variable(i) for i in range(4)]
    theta = y-6*x*z

    # Signed nonlinear terms in the six equations of X.
    nonlinear = [
        x**2*(3*theta+2*x*z),
        -24*x*theta*(theta+x*z),
        -2*theta*(2*theta+3*x*z),
        -12*x**2*theta**2*(3*theta+2*x*z),
        -2*x*theta**2*(7*theta+6*x*z),
        -4*x**2*theta**3*(3*theta+2*x*z),
    ]

    # J = DX. The field's constant target does not affect its derivative.
    J = [[Poly(-1 if i == j else 0) for j in range(6)] for i in range(6)]
    for i in range(6):
        for j in range(3):
            J[i][j] = J[i][j] + nonlinear[i].derivative(j)
    J[1][3] += 1  # +a in the second equation
    J[2][4] += 1  # +b in the third equation
    J[4][5] += 1  # +c in the fifth equation
    A = [[(lam if i == j else Poly(0))-J[i][j]
          for j in range(6)] for i in range(6)]

    @lru_cache(None)
    def determinant(columns: Tuple[int, ...]) -> Poly:
        row = 6-len(columns)
        if not columns:
            return Poly(1)
        result = Poly(0)
        for position, column in enumerate(columns):
            if not A[row][column]:
                continue
            minor_columns = columns[:position]+columns[position+1:]
            term = A[row][column]*determinant(minor_columns)
            result = result + (term if position % 2 == 0 else -term)
        return result

    computed = determinant(tuple(range(6)))
    expected = (lam+1)**6
    if computed.terms != expected.terms:
        raise AssertionError(f"Characteristic identity failed: {(computed-expected).terms}")

    def field(point):
        x, y, z, a, b, c = map(Fraction, point)
        q = y-6*x*z
        h1 = x*x*(3*q+2*x*z)
        h2 = 24*x*q*(q+x*z)
        h3 = 2*q*(2*q+3*x*z)
        h4 = 12*x*x*q*q*(3*q+2*x*z)
        h5 = 2*x*q*q*(7*q+6*x*z)
        h6 = 4*x*x*q**3*(3*q+2*x*z)
        return (-x+h1, 1-y+a-h2, -z+b-h3,
                -a-h4, -b+c-h5, -c-h6)

    points = [
        (0,1,-4,0,0,0),
        (-1,-11,2,12,-6,4),
        (1,7,Fraction(5,4),-3,Fraction(-3,2),Fraction(1,2)),
    ]
    values = [field(p) for p in points]
    if not all(all(v == 0 for v in value) for value in values):
        raise AssertionError("An equilibrium identity failed")

    report = {
        "status": "independent standard-library verification passed",
        "polynomial_ring": "Z[x,y,z,lambda]",
        "method": "sparse integer polynomial arithmetic and exact Laplace expansion",
        "verified_identity": "det(lambda I - DX) = (lambda + 1)^6",
        "characteristic_coefficients_ascending": [1,6,15,20,15,6,1],
        "number_of_nonzero_residual_terms": len((computed-expected).terms),
        "three_distinct_equilibria_verified": True,
        "external_libraries": [],
        "floating_point_used": False,
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
