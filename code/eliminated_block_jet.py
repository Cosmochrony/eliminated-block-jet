"""Exact verification for the eliminated-block-jet note (sections 2-5 and 7).

Everything below is exact (SymPy rationals and algebraic numbers, no floating-point sampling).
The checks follow the statements of the note:

  Part A   Theorem (constrained jet): projector-jet identities for a general second-order transport
           U(s) = 1 + sA + s^2 (A^2 + B)/2, and the family exp(s^2 K) Q0 exp(-s^2 K) that no constant
           generator reproduces.
  Part A2  Theorem (ii)-(iii): parity decomposition of the tangent under an antiunitary parity, the odd part
           of Q2 for a mixed-parity constant generator, and alternating parities for an equivariant family.
  Part B   Lemma (parity propagation with a moving embedding): product rule for E1 against the exact family,
           parity of E0 and E1 on a generic covariant family, the C^2 instance in which the frozen embedding is
           wrong, and a covariant G = Pi_S D with D not J-even whose block is not equivariant.
  Part C   Proposition (split source): general form of a Hermitian antiunitary-odd operator on C^3_gen,
           evenness of i R_mix, vanishing of both R_mix pairings, a non-zero tangent annihilated by a frozen
           compression, the normal-form criterion u'(0) = -2c (E1)_{++}, and the Remark on constant blocks:
           a covariant moving embedding with 1-P(s) constant and u'(0) = 19/8, the formula
           u'(0) = tr([E0^2, J_3] kappa) / <J_3, J_3>, and its vanishing in the normal form.
  Part D   Theorem (realisation): the exact covariant moving family with E0^2 = diag(1,1/2,1/2),
           u'(0) = 11 2^{1/4} / 9, no first-order mixing, no first-order singlet leakage, frozen-embedding rate
           79 2^{1/4} / 90, independence of E0 from a general admissible H, the J-even s^2 coefficient of
           E_Pi(s)^2 with a non-zero (+,-) entry (second-order mixing), and the vanishing of the split
           functional at H = 0.

Convention on C^3_gen: basis (e_0, e_+, e_-), J_3 = diag(0, 1, -1), R_mix = [[0,0,0],[0,0,1],[0,-1,0]],
antiunitary parity J^(2) = S o conj with S : e_0 -> -e_0, e_+ <-> e_-.

Run:  python code/eliminated_block_jet.py     (SymPy >= 1.12; about a minute)
"""

import sympy as sp
from sympy import I, Matrix, Rational, conjugate, eye, simplify, sqrt, zeros

s = sp.symbols("s", real=True)


def dag(M):
    return M.H


def cj(M):
    return M.applyfunc(conjugate)


def comm(A, B):
    return A * B - B * A


def anticomm(A, B):
    return A * B + B * A


def hs(X, Y):
    return (dag(X) * Y).trace()


def is_zero(M):
    return M.applyfunc(simplify) == zeros(*M.shape)


def coeff(M, k):
    return M.applyfunc(lambda x: sp.expand(x).coeff(s, k)).applyfunc(simplify)


def rank1_proj(v):
    v = Matrix(v)
    return (v * dag(v)) / (dag(v) * v)[0]


def hstack(*Ms):
    return Matrix.hstack(*Ms)


def vstack(*Ms):
    return Matrix.vstack(*Ms)


# Generation-triplet structure (frozen convention).
S3 = Matrix([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
J3 = sp.diag(0, 1, -1)
RMIX = Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])


def parity_gen(X):
    """Antiunitary parity on operators of C^3_gen: J X J^{-1} = S conj(X) S."""
    return S3 * cj(X) * S3


CHECKS = {}


def record(name, ok):
    CHECKS[name] = bool(ok)


# ---------------------------------------------------------------------------------------------------
# Part A: projector-jet identities for a general second-order transport.
# ---------------------------------------------------------------------------------------------------
def part_A():
    n = 4
    In = eye(n)
    instances = [
        (Matrix([[1, 2 + I, 0, 1], [0, 1, I, 2]]).T, Matrix([[0, 1, 2, I], [-1, 0, 1, 0], [-2, -1, 0, 3], [I, 0, -3, 0]]),
         Matrix([[I, 2, 0, 1 - I], [-2, 0, 1, 0], [0, -1, 2 * I, 1], [-1 - I, 0, -1, 0]])),
        (Matrix([[2, 0, 1, I], [1, 1, 0, 0]]).T, Matrix([[0, I, 1, 0], [I, 0, 0, 2], [-1, 0, 0, 1], [0, -2, -1, 0]]),
         Matrix([[0, 1, 1, 1], [-1, 0, 2, 0], [-1, -2, 0, 3], [-1, 0, -3, I]])),
    ]
    for idx, (V, A, B) in enumerate(instances):
        # Q0: orthogonal projector onto the span of the columns of V.
        Q0 = V * (dag(V) * V).inv() * dag(V)
        A = A - dag(A)  # anti-Hermitian
        B = B - dag(B)
        U = In + s * A + s**2 * (A * A + B) / 2
        Q = U * Q0 * dag(U)
        Q1, Q2 = coeff(Q, 1), coeff(Q, 2)
        tag = f"A.inst{idx}"
        record(f"{tag}.Q0_projector", is_zero(Q0 * Q0 - Q0) and is_zero(Q0 - dag(Q0)))
        record(f"{tag}.U_unitary_to_order2", is_zero(coeff(U * dag(U), 1)) and is_zero(coeff(U * dag(U), 2)))
        record(f"{tag}.Q1_is_comm_A_Q0", is_zero(Q1 - comm(A, Q0)))
        record(f"{tag}.Q1_offdiag_upper", is_zero(Q0 * Q1 * Q0))
        record(f"{tag}.Q1_offdiag_lower", is_zero((In - Q0) * Q1 * (In - Q0)))
        record(f"{tag}.Q1_nonzero", not is_zero(Q1))
        record(f"{tag}.order1_constraint", is_zero(anticomm(Q0, Q1) - Q1))
        record(f"{tag}.order2_constraint", is_zero(Q0 * Q2 + Q2 * Q0 + Q1 * Q1 - Q2))
        record(f"{tag}.Q2_block_upper", is_zero(Q0 * Q2 * Q0 + Q0 * Q1 * Q1 * Q0))
        record(f"{tag}.Q2_block_lower", is_zero((In - Q0) * Q2 * (In - Q0) - (In - Q0) * Q1 * Q1 * (In - Q0)))
        record(f"{tag}.Q2_formula", is_zero(Q2 - comm(B, Q0) / 2 - comm(A, comm(A, Q0)) / 2))
        record(f"{tag}.B_matters", not is_zero(comm(B, Q0)))
    # A smooth projector family with Q1 = 0 and Q2 != 0: no constant generator.
    Q0 = sp.diag(1, 0)
    K = Matrix([[0, -1], [1, 0]])
    U = eye(2) + s**2 * K  # exp(s^2 K) to second order
    Q = U * Q0 * dag(U)
    record("A.s2family.Q1_zero", is_zero(coeff(Q, 1)))
    record("A.s2family.Q2_is_comm_K_Q0", is_zero(coeff(Q, 2) - comm(K, Q0)) and not is_zero(comm(K, Q0)))
    # Any constant generator with [A, Q0] = 0 gives Q2 = [A,[A,Q0]]/2 = 0: contradiction with Q2 != 0.
    a, b = sp.symbols("a b", real=True)
    Agen = Matrix([[I * a, 0], [0, I * b]])  # general anti-Hermitian commuting with diag(1,0)
    record("A.s2family.no_constant_generator", is_zero(comm(Agen, comm(Agen, Q0))))


# ---------------------------------------------------------------------------------------------------
# Part A2: parity decomposition under an antiunitary parity; equivariance is needed for alternation.
# ---------------------------------------------------------------------------------------------------
def part_A2():
    # Antiunitary parity on C^4: J = Sig o conj with Sig a real symmetric involution.
    Sig = Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])

    def par(X):
        return Sig * cj(X) * Sig

    def odd(X):
        return (X - par(X)) / 2

    def even(X):
        return (X + par(X)) / 2

    # Even projector: Q0 = diag(q, conj(q)) with q a 2x2 projector.
    q = rank1_proj([3, 4 + I])
    Q0 = vstack(hstack(q, zeros(2)), hstack(zeros(2), cj(q)))
    record("A2.Q0_even_projector", is_zero(Q0 * Q0 - Q0) and is_zero(par(Q0) - Q0))
    # Generic anti-Hermitian generator with both parities.
    Araw = Matrix([[I, 2 + I, 1, 3], [0, 2 * I, I, 1], [1, 0, 0, 2 - I], [I, 1, 1, 0]])
    A = Araw - dag(Araw)
    Ap, Am = even(A), odd(A)
    record("A2.parity_split", is_zero(A - Ap - Am) and is_zero(par(Ap) - Ap) and is_zero(par(Am) + Am))
    record("A2.comm_Aplus_even", is_zero(par(comm(Ap, Q0)) - comm(Ap, Q0)))
    record("A2.comm_Aminus_odd", is_zero(par(comm(Am, Q0)) + comm(Am, Q0)))
    record("A2.odd_part_Q1_is_comm_Aminus", is_zero(odd(comm(A, Q0)) - comm(Am, Q0)))
    record("A2.even_only_generator_sources_nothing",
           (not is_zero(comm(Ap, Q0))) and is_zero(odd(comm(Ap, Q0))))
    # Mixed-parity constant generator: Q2 = [A,[A,Q0]]/2 has an odd part (no alternating parity).
    Q2 = comm(A, comm(A, Q0)) / 2
    record("A2.mixed_generator_Q2_has_odd_part", not is_zero(odd(Q2)))
    # The 3x3 real instance of the note (linear J = diag(1,1,-1) and antiunitary J o conj coincide on real data).
    Q0r = sp.diag(1, 0, 0)
    Jr = sp.diag(1, 1, -1)
    Ar = Matrix([[0, -1, -1], [1, 0, 0], [1, 0, 0]])
    Q1r = comm(Ar, Q0r)
    Q2r = comm(Ar, Q1r) / 2
    record("A2.note_instance_Q1_mixed", not is_zero(Jr * Q1r * Jr + Q1r) and not is_zero(Jr * Q1r * Jr - Q1r))
    record("A2.note_instance_Q2_odd_23_entries", Q2r == Matrix([[-2, 0, 0], [0, 1, 1], [0, 1, 1]])
           and not is_zero(Jr * Q2r * Jr - Q2r))
    # Equivariant family (odd generator): Q1 odd, Q2 even, to second order.
    U = eye(4) + s * Am + s**2 * Am * Am / 2
    Q = U * Q0 * dag(U)
    record("A2.equivariant_Q1_odd", is_zero(par(coeff(Q, 1)) + coeff(Q, 1)))
    record("A2.equivariant_Q2_even", is_zero(par(coeff(Q, 2)) - coeff(Q, 2)))
    record("A2.equivariant_family_to_order2",
           is_zero(coeff(par(Q), 1) + coeff(Q, 1)) and is_zero(coeff(par(Q), 2) - coeff(Q, 2)))


# ---------------------------------------------------------------------------------------------------
# Part B: moving embedding. Product rule versus frozen embedding.
# ---------------------------------------------------------------------------------------------------
def part_B():
    t = sp.symbols("t", real=True)
    Pi = Matrix([[sp.cos(t), sp.sin(t)]])
    Q = eye(2) - Pi.T * Pi
    true_residue = sp.simplify((Pi * Q * Pi.T)[0])
    frozen = sp.simplify((Pi.subs(t, 0) * Q * Pi.subs(t, 0).T)[0])
    record("B.C2_true_residue_zero", true_residue == 0)
    record("B.C2_frozen_residue_is_sin2", sp.simplify(frozen - sp.sin(t) ** 2) == 0)
    # Product rule against the exact s-coefficient of a generic covariant family (structure of Part D, random data).
    Pi0 = hstack(eye(3), S3) / sqrt(2)
    Draw = Matrix([[1 + I, 2, I], [2, 3 - I, Rational(1, 2)], [I, Rational(1, 2), -1]])
    Dp = (Draw + Draw.T) / 2  # complex symmetric
    D = vstack(hstack(zeros(3), dag(Dp)), hstack(Dp, zeros(3)))
    a = Matrix([[I, 1 + I, 2], [-1 + I, 2 * I, I], [-2, I, 0]])
    b = Matrix([[1, I, 2], [I, 3, 1 - I], [2, 1 - I, I]])
    A = vstack(hstack(a, b), hstack(-dag(b), -cj(a)))
    U = eye(6) + s * A + s**2 * A * A / 2
    Pis = Pi0 * dag(U)
    Qs = eye(6) - dag(Pis) * Pis
    Es = -Pis * D * Qs * D * dag(Pis)
    Pi1, Q0, Q1, E1 = coeff(Pis, 1), coeff(Qs, 0), coeff(Qs, 1), coeff(Es, 1)
    formula = -(Pi1 * D * Q0 * D * dag(Pi0) + Pi0 * D * Q1 * D * dag(Pi0) + Pi0 * D * Q0 * D * dag(Pi1))
    record("B.product_rule_matches_family", is_zero(E1 - formula))
    record("B.Q1_from_embedding", is_zero(Q1 + dag(Pi1) * Pi0 + dag(Pi0) * Pi1))
    record("B.frozen_differs_from_moving", not is_zero(E1 + Pi0 * D * Q1 * D * dag(Pi0)))
    E0 = coeff(Es, 0)
    record("B.generic_E0_even_E1_odd", is_zero(parity_gen(E0) - E0) and is_zero(parity_gen(E1) + E1))
    # Covariance of G alone does not suffice: D = diag(1,1,1,1,1,-1) is Hermitian and unitary but not J-even.
    Sig = vstack(hstack(zeros(3), eye(3)), hstack(eye(3), zeros(3)))
    Dbad = sp.diag(1, 1, 1, 1, 1, -1)
    G = Pi0 * dag(U)
    record("B.G_covariant_to_order2",
           is_zero(coeff(S3 * cj(G) * Sig, 1) + coeff(G, 1)) and is_zero(coeff(S3 * cj(G) * Sig, 2) - coeff(G, 2)))
    Pi_bad = Pi0 * Dbad
    Qbad0 = eye(6) - dag(Pi_bad) * Pi_bad
    record("B.D_not_even", not is_zero(Sig * cj(Dbad) * Sig - Dbad))
    record("B.block_not_equivariant_when_D_not_even", not is_zero(Sig * cj(Qbad0) * Sig - Qbad0))


# ---------------------------------------------------------------------------------------------------
# Part C: antiunitary-odd Hermitian operators on C^3_gen; criterion for the split rate.
# ---------------------------------------------------------------------------------------------------
def part_C():
    re = sp.symbols("r0:9", real=True)
    im = sp.symbols("m0:9", real=True)
    X = Matrix(3, 3, lambda i, j: re[3 * i + j] + I * im[3 * i + j])
    eqs = list((X - dag(X)).applyfunc(sp.expand)) + list((parity_gen(X) + X).applyfunc(sp.expand))
    sol = sp.solve([sp.re(e) for e in eqs] + [sp.im(e) for e in eqs], list(re) + list(im), dict=True)[0]
    Xodd = X.subs(sol).applyfunc(simplify)
    free = sorted(Xodd.free_symbols, key=str)
    record("C.odd_hermitian_family_three_real_parameters", len(free) == 3)
    record("C.odd_hermitian_00_zero", Xodd[0, 0] == 0)
    record("C.odd_hermitian_plusminus_zero", Xodd[1, 2] == 0 and Xodd[2, 1] == 0)
    record("C.odd_hermitian_minusminus_is_minus_plusplus", simplify(Xodd[2, 2] + Xodd[1, 1]) == 0)
    record("C.odd_hermitian_singlet_row_conjugate", simplify(Xodd[0, 2] - conjugate(Xodd[0, 1])) == 0)
    record("C.Rmix_pairing_zero", simplify(hs(RMIX, Xodd)) == 0)
    record("C.iRmix_pairing_zero", simplify(hs(I * RMIX, Xodd)) == 0)
    record("C.Rmix_odd", is_zero(parity_gen(RMIX) + RMIX))
    record("C.iRmix_even_not_odd", is_zero(parity_gen(I * RMIX) - I * RMIX))
    # Criterion in the normal form E0 = -diag(1, c, c): u'(0) = -2c (E1)_{++}.
    c = sp.symbols("c", positive=True)
    E0 = -sp.diag(1, c, c)
    anti = anticomm(E0, Xodd)
    u1 = simplify(hs(J3, anti) / hs(J3, J3))
    record("C.criterion_u1_is_minus_2c_E1pp", simplify(u1 + 2 * c * Xodd[1, 1]) == 0)
    record("C.singlet_row_of_anticommutator", simplify(anti[0, 1] + (1 + c) * Xodd[0, 1]) == 0)
    # Non-zero odd tangent, zero split rate (frozen even compression G = diag(1,0), 2x2 model).
    Q0 = sp.diag(1, 0)
    K = Matrix([[0, -1], [1, 0]])
    G = sp.diag(1, 0)
    Q1 = comm(K, Q0)
    record("C.nonzero_tangent_annihilated_by_frozen_compression", (not is_zero(Q1)) and is_zero(G * Q1 * dag(G)))
    # Remark on constant blocks: a covariant moving embedding with constant block and non-zero split rate.
    Sig = vstack(hstack(zeros(3), eye(3)), hstack(eye(3), zeros(3)))
    Pi0 = hstack(eye(3), S3) / sqrt(2)
    Dp = Matrix([[1 + I, 1, 0], [1, 2, I], [0, I, -1]])
    D = vstack(hstack(zeros(3), dag(Dp)), hstack(Dp, zeros(3)))
    record("C.constant_block.D_even", is_zero(Sig * cj(D) * Sig - D) and is_zero(D - dag(D)))
    record("C.constant_block.Rmix_odd_antihermitian", is_zero(parity_gen(RMIX) + RMIX) and is_zero(RMIX + dag(RMIX)))
    Pis = (eye(3) + s * RMIX + s**2 * RMIX * RMIX / 2) * Pi0
    Qs = eye(6) - dag(Pis) * Pis
    record("C.constant_block.Q1_Q2_zero", is_zero(coeff(Qs, 1)) and is_zero(coeff(Qs, 2)))
    Es = -Pis * D * Qs * D * dag(Pis)
    E0, E1 = coeff(Es, 0), coeff(Es, 1)
    record("C.constant_block.E1_is_comm", is_zero(E1 - comm(RMIX, E0)))
    u1 = simplify(hs(J3, anticomm(E0, E1)) / hs(J3, J3))
    record("C.constant_block.u1_is_19_over_8", u1 == Rational(19, 8))
    record("C.constant_block.E0sq_plusminus", simplify((E0 * E0)[1, 2] - (Rational(19, 16) + 21 * I)) == 0)
    # Odd orders of the mixing datum vanish (m(-s) = m(s)); checked to third order on this covariant family.
    U3 = eye(3) + s * RMIX + s**2 * RMIX * RMIX / 2 + s**3 * RMIX**3 / 6
    Pis3 = U3 * Pi0
    Es3 = -Pis3 * D * (eye(6) - dag(Pis3) * Pis3) * D * dag(Pis3)
    m_s = sp.expand((Es3 * Es3)[1, 2])
    record("C.constant_block.mixing_orders_1_3_zero",
           simplify(m_s.coeff(s, 1)) == 0 and simplify(m_s.coeff(s, 3)) == 0)
    record("C.constant_block.mixing_orders_0_2_nonzero",
           simplify(m_s.coeff(s, 0)) != 0 and simplify(m_s.coeff(s, 2)) != 0)
    record("C.constant_block.formula", simplify(u1 - (comm(E0 * E0, J3) * RMIX).trace() / hs(J3, J3)) == 0)
    # Normal form: [E0^2, J_3] = 0, so every constant-block motion has zero first-order rate.
    mr = sp.symbols("mr0:9", real=True)
    mi = sp.symbols("mi0:9", real=True)
    mm = Matrix(3, 3, lambda i, j: mr[3 * i + j] + I * mi[3 * i + j])
    mm = mm - dag(mm)
    En = -sp.diag(1, 1 / sqrt(2), 1 / sqrt(2))
    E1n = comm(mm, En)
    record("C.normal_form_constant_block_rate_zero", simplify(hs(J3, anticomm(En, E1n))) == 0)


# ---------------------------------------------------------------------------------------------------
# Part D: exact covariant moving family realising the even normal form with a non-zero split rate.
# ---------------------------------------------------------------------------------------------------
def part_D():
    Sig = vstack(hstack(zeros(3), eye(3)), hstack(eye(3), zeros(3)))

    def par_H(X):  # antiunitary parity on S_L (+) S_R
        return Sig * cj(X) * Sig

    def par_Pi(X):  # intertwined action on maps H -> C^3_gen
        return S3 * cj(X) * Sig

    Pi0 = hstack(eye(3), S3) / sqrt(2)
    record("D.Pi0_coisometry", is_zero(Pi0 * dag(Pi0) - eye(3)))
    record("D.Pi0_equivariant", is_zero(par_Pi(Pi0) - Pi0))
    Q0 = eye(6) - dag(Pi0) * Pi0
    r = 2 ** Rational(-1, 4)
    K = Matrix([[2, 0, 0], [0, 0, 2 * r], [0, 2 * r, 0]])
    z = Rational(1, 2) + I * Rational(1, 3)
    H = Matrix([[1, 0, 0], [0, Rational(1, 3), z], [0, conjugate(z), Rational(1, 3)]])
    record("D.H_admissible", is_zero(H - dag(H)) and is_zero(S3 * H * S3 - cj(H)))

    def dirac(Hh):
        Dp = S3 * (-I * K / 2 + Hh)
        return Dp, vstack(hstack(zeros(3), dag(Dp)), hstack(Dp, zeros(3)))

    Dp, D = dirac(H)
    record("D.Dplus_symmetric", is_zero(Dp - Dp.T))
    record("D.D_hermitian_even", is_zero(D - dag(D)) and is_zero(par_H(D) - D))
    a = Matrix([[I * Rational(1, 2), 0, 0],
                [0, I * Rational(1, 3), Rational(2, 3) + I],
                [0, -(Rational(2, 3) - I), -I * Rational(1, 5)]])
    b = Matrix([[Rational(1, 2) - I, 0, 0],
                [0, Rational(1, 3) + I * Rational(1, 2), 1 - I * Rational(1, 4)],
                [0, 1 - I * Rational(1, 4), -Rational(1, 7) + I * Rational(2, 3)]])
    A = vstack(hstack(a, b), hstack(-dag(b), -cj(a)))
    record("D.A_antihermitian_odd", is_zero(A + dag(A)) and is_zero(par_H(A) + A))

    def family(D_, A_):
        U = eye(6) + s * A_ + s**2 * A_ * A_ / 2
        Pis = Pi0 * dag(U)
        Qs = eye(6) - dag(Pis) * Pis
        Es = -Pis * D_ * Qs * D_ * dag(Pis)
        return Pis, Qs, Es

    Pis, Qs, Es = family(D, A)
    E0, E1 = coeff(Es, 0), coeff(Es, 1)
    Pi1, Q1 = coeff(Pis, 1), coeff(Qs, 1)
    c = 1 / sqrt(2)
    record("D.E0_normal_form", E0 == -sp.diag(1, c, c))
    record("D.E0_squared", (E0 * E0).applyfunc(simplify) == sp.diag(1, Rational(1, 2), Rational(1, 2)))
    record("D.covariance_Pi1_odd_Q1_odd", is_zero(par_Pi(Pi1) + Pi1) and is_zero(par_H(Q1) + Q1))
    record("D.Q1_is_comm_A_Q0", is_zero(Q1 - comm(A, Q0)))
    record("D.E1_hermitian_odd", is_zero(E1 - dag(E1)) and is_zero(parity_gen(E1) + E1))
    formula = -(Pi1 * D * Q0 * D * dag(Pi0) + Pi0 * D * Q1 * D * dag(Pi0) + Pi0 * D * Q0 * D * dag(Pi1))
    record("D.E1_product_rule", is_zero(E1 - formula))
    record("D.E1_value", E1 == -Rational(11, 18) * 2 ** Rational(3, 4) * sp.diag(0, 1, -1))
    anti = anticomm(E0, E1).applyfunc(simplify)
    u1 = simplify(hs(J3, anti) / hs(J3, J3))
    v1 = simplify(hs(RMIX, anti) / hs(RMIX, RMIX))
    w1 = simplify(hs(I * RMIX, anti))
    record("D.u1_value", u1 == Rational(11, 9) * 2 ** Rational(1, 4))
    record("D.u1_criterion", simplify(u1 + 2 * c * E1[1, 1]) == 0)
    record("D.v1_zero_and_iRmix_zero", v1 == 0 and w1 == 0)
    record("D.no_singlet_leakage", anti[0, 1] == 0 and anti[0, 2] == 0)
    E1aux = (-Pi0 * D * Q1 * D * dag(Pi0)).applyfunc(simplify)
    u1aux = simplify(hs(J3, anticomm(E0, E1aux)) / hs(J3, J3))
    record("D.frozen_rate_value", u1aux == Rational(79, 90) * 2 ** Rational(1, 4))
    record("D.frozen_rate_differs", simplify(u1 - u1aux) != 0)
    # H = 0: the split functional vanishes on every J-odd generator (21 real dimensions, symbolic).
    _, D0 = dirac(zeros(3))
    record("D.H0_normal_form_too", coeff(family(D0, A)[2], 0) == -sp.diag(1, c, c))
    ar = sp.symbols("x0:9", real=True)
    ai = sp.symbols("y0:9", real=True)
    br = sp.symbols("p0:9", real=True)
    bi = sp.symbols("q0:9", real=True)
    ag = Matrix(3, 3, lambda i, j: ar[3 * i + j] + I * ai[3 * i + j])
    bg = Matrix(3, 3, lambda i, j: br[3 * i + j] + I * bi[3 * i + j])
    ag = ag - dag(ag)
    bg = bg + bg.T
    Ag = vstack(hstack(ag, bg), hstack(-dag(bg), -cj(ag)))
    record("D.generic_generator_odd", is_zero(par_H(Ag) + Ag) and is_zero(Ag + dag(Ag)))
    Q1g = comm(Ag, Q0)
    Pi1g = -Pi0 * Ag
    E1g = -(Pi1g * D0 * Q0 * D0 * dag(Pi0) + Pi0 * D0 * Q1g * D0 * dag(Pi0) + Pi0 * D0 * Q0 * D0 * dag(Pi1g))
    u1g = simplify(hs(J3, anticomm(E0, E1g)) / hs(J3, J3))
    record("D.H0_split_functional_vanishes", u1g == 0)
    # With H != 0 the same generic functional is non-zero.
    E1h = -(Pi1g * D * Q0 * D * dag(Pi0) + Pi0 * D * Q1g * D * dag(Pi0) + Pi0 * D * Q0 * D * dag(Pi1g))
    u1h = simplify(hs(J3, anticomm(E0, E1h)) / hs(J3, J3))
    record("D.H_split_functional_nonzero", u1h != 0)
    # E0 is independent of a general admissible H (Hermitian, S H S = conj(H)).
    hr = sp.symbols("hr0:9", real=True)
    hi = sp.symbols("hi0:9", real=True)
    Hg = Matrix(3, 3, lambda i, j: hr[3 * i + j] + I * hi[3 * i + j])
    eqs = list((Hg - dag(Hg)).applyfunc(sp.expand)) + list((S3 * Hg * S3 - cj(Hg)).applyfunc(sp.expand))
    sol = sp.solve([sp.re(e) for e in eqs] + [sp.im(e) for e in eqs], list(hr) + list(hi), dict=True)[0]
    Hg = Hg.subs(sol).applyfunc(simplify)
    _, Dg = dirac(Hg)
    E0g = (-Pi0 * Dg * Q0 * Dg * dag(Pi0)).applyfunc(simplify)
    record("D.E0_independent_of_general_H", len(Hg.free_symbols) == 6 and E0g == -sp.diag(1, c, c))
    # Second order: the s^2 coefficient of E_Pi(s)^2 is J-even and carries a non-zero (+,-) entry.
    E2 = coeff(Es, 2)
    F2 = (E1 * E1 + E0 * E2 + E2 * E0).applyfunc(simplify)
    record("D.second_order_even", is_zero(parity_gen(F2) - F2))
    record("D.second_order_mixing_nonzero", simplify(F2[1, 2]) != 0 and simplify(hs(RMIX, F2)) != 0)
    extras2 = {"m2": sp.nsimplify(simplify(F2[1, 2]))}
    return {"u1": u1, "u1aux": u1aux, "E1": E1, "m2": extras2["m2"]}


def main():
    part_A()
    part_A2()
    part_B()
    part_C()
    extras = part_D()
    failed = [k for k, v in CHECKS.items() if not v]
    for name, ok in CHECKS.items():
        print(f"  [{'OK ' if ok else 'FAIL'}] {name}")
    print()
    print("Realisation:  E1 =", list(extras["E1"]))
    print("              u'(0) =", extras["u1"], " frozen-embedding rate =", extras["u1aux"])
    print("              s^2 coefficient of the (+,-) entry of E^2 =", extras["m2"])
    print()
    print(f"{len(CHECKS)} checks, {len(failed)} failed")
    assert not failed, failed
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
