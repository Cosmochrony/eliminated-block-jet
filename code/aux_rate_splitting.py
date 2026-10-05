"""Exact splitting of the auxiliary (frozen-embedding) rate on the witness of Theorem (realisation).

The auxiliary first-order residue is E1_aux = -G0 Qdot(0) G0^*, G0 = Pi_0 D, and its rate is
<J_3, {E0, E1_aux}> / <J_3, J_3>.  Splitting Qdot(0) = Q1 into its chirally diagonal part (LL and RR blocks)
and its chirally off-diagonal part (LR and RL blocks), the rate is the sum of two exact contributions.
Normalisation T = <J_3, G0 Qdot(0) G0^*>:  rate = sqrt(2) T / <J_3, J_3>  (E0 = -diag(1, c, c), c = 2^-1/2).

Run:  python -W error code/aux_rate_splitting.py     (SymPy >= 1.12)
"""

import sympy as sp
from sympy import I, Matrix, Rational, eye, sqrt, zeros

S = Matrix([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
K = Matrix([[2, 0, 0], [0, 0, 2 ** Rational(3, 4)], [0, 2 ** Rational(3, 4), 0]])
z = Rational(1, 2) + I / 3
H = Matrix([[1, 0, 0], [0, Rational(1, 3), z], [0, sp.conjugate(z), Rational(1, 3)]])
Dp = S * (-I / 2 * K + H)
D = Matrix(sp.BlockMatrix([[zeros(3), Dp.H], [Dp, zeros(3)]]))
Pi0 = Matrix(sp.BlockMatrix([[eye(3), S]])) / sqrt(2)
a = Matrix([[I / 2, 0, 0], [0, I / 3, Rational(2, 3) + I], [0, -Rational(2, 3) + I, -I / 5]])
b = Matrix([[Rational(1, 2) - I, 0, 0],
            [0, Rational(1, 3) + I / 2, 1 - I / 4],
            [0, 1 - I / 4, -Rational(1, 7) + 2 * I / 3]])
A = Matrix(sp.BlockMatrix([[a, b], [-b.H, -a.conjugate()]]))

Pi1 = -Pi0 * A                      # first-order coefficient of Pi_S(s) = Pi_0 exp(-sA)
Q1 = -(Pi1.H * Pi0 + Pi0.H * Pi1)   # Qdot(0), first-order coefficient of 1 - P(s)
G0 = Pi0 * D
J3 = sp.diag(0, 1, -1)
E0 = -sp.diag(1, 1 / sqrt(2), 1 / sqrt(2))


def hs(X, Y):
    return sp.simplify((X.H * Y).trace())


def chiral_part(M, diagonal):
    out = zeros(6)
    for (r, c) in ((0, 0), (1, 1)) if diagonal else ((0, 1), (1, 0)):
        out[3 * r:3 * r + 3, 3 * c:3 * c + 3] = M[3 * r:3 * r + 3, 3 * c:3 * c + 3]
    return out


Qd, Qo = chiral_part(Q1, True), chiral_part(Q1, False)
assert sp.simplify(Qd + Qo - Q1) == zeros(6)

n = hs(J3, J3)
assert n == 2
T = {k: sp.nsimplify(sp.simplify(hs(J3, G0 * X * G0.H)))
     for k, X in (("full", Q1), ("diag", Qd), ("offdiag", Qo))}
rate = {k: sp.simplify(sqrt(2) * v / n) for k, v in T.items()}

E1aux = -G0 * Q1 * G0.H
rate_aux = sp.simplify(hs(J3, E0 * E1aux + E1aux * E0) / n)

two34, two14 = 2 ** Rational(3, 4), 2 ** Rational(1, 4)
checks = [
    ("T diagonal part = 11 2^(3/4)/18", sp.simplify(T["diag"] - 11 * two34 / 18) == 0),
    ("T off-diagonal part = 4 2^(3/4)/15", sp.simplify(T["offdiag"] - 4 * two34 / 15) == 0),
    ("T total = 79 2^(3/4)/90", sp.simplify(T["full"] - 79 * two34 / 90) == 0),
    ("rate diagonal part = 11 2^(1/4)/18", sp.simplify(rate["diag"] - 11 * two14 / 18) == 0),
    ("rate off-diagonal part = 4 2^(1/4)/15", sp.simplify(rate["offdiag"] - 4 * two14 / 15) == 0),
    ("auxiliary rate = 79 2^(1/4)/90", sp.simplify(rate_aux - 79 * two14 / 90) == 0),
    ("auxiliary rate = sum of the two parts", sp.simplify(rate_aux - rate["diag"] - rate["offdiag"]) == 0),
    ("off-diagonal part is non-zero", sp.simplify(rate["offdiag"]) != 0),
]
for name, ok in checks:
    print(("PASS  " if ok else "FAIL  ") + name)
print("T:", T)
print("rates:", rate, "aux:", rate_aux)
assert all(ok for _, ok in checks)
print("%d/%d checks passed" % (len(checks), len(checks)))
