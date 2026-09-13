"""Exact verification of the chiral-diagonal defect-rate identity (section 6 of the note).

Exact complex-rational / symbolic arithmetic, no floating-point sampling.

Frame: S_Pi = S_L (+) S_R, each C^2; antiunitary Born--Infeld parity J = Sigma_tau o conj with
Sigma_tau = [[0, tau], [tau^{-1}, 0]]; tau = I is the chirality swap.
Transported defect: Delta_chi(P) = pi_LL - tau conj(pi_RR) tau^{-1}.

Checked:
  (1) Lemma (tau = I): [Pi_{J-odd}(1-P)]_LL = (1/2) Delta_chi(P), and its RR block is -(1/2) conj(Delta_chi(P)).
  (2) Proposition (transported frame): the same LL identity for one real orthogonal tau, with Sigma_tau an involution.
  (3) At an equivariant point (pi_LL = conj(pi_RR), tau = I) the defect and the odd diagonal blocks vanish.
  (4) Chiral block products for Pi_S = 1: E_LL = -D_- pi_RR D_+, E_RR = -D_+ pi_LL D_- with E = -D (1-P) D.
  (5) Rate identity: for a J-even orthogonal projector Q0 and a J-odd anti-Hermitian generator A, the Grassmann
      tangent Q1 = [A, Q0] is odd and off-diagonal with respect to Q0, and its LL block equals
      (1/2) d/ds Delta_chi at s = 0 for the family exp(sA) Q0 exp(-sA).

Not checked, and not claimed: any identification of a unique generator, any localisation of the odd tangent
onto its LL block, the passage from the tangent to the split coefficient u, its sign, or its normalisation.

Run:  python code/front_d0_generator.py     (SymPy >= 1.12; under a second)
"""

import sympy as sp

s = sp.symbols("s", real=True)
I2 = sp.eye(2)
Z2 = sp.zeros(2, 2)


def block(a, b, c, d):
    return sp.Matrix(sp.BlockMatrix([[a, b], [c, d]]))


def is_zero(M):
    return sp.Matrix(M).applyfunc(sp.simplify) == sp.zeros(*M.shape)


def comm(A, B):
    return A * B - B * A


def sigma_tau(tau):
    return block(Z2, tau, tau.inv(), Z2)


def parity(X, tau=I2):
    St = sigma_tau(tau)
    return St * X.conjugate() * St


def odd(X, tau=I2):
    return (X - parity(X, tau)) / 2


def hermitian(r1, r2, x, y):
    return sp.Matrix([[r1, x + sp.I * y], [x - sp.I * y, r2]])


def main():
    checks = {}

    pi_LL = hermitian(3, 7, 1, 2)
    pi_RR = hermitian(5, 11, 4, 3)
    pi_LR = sp.Matrix([[1 + sp.I, 2], [sp.I, 3 - sp.I]])
    elim = block(pi_LL, pi_LR, pi_LR.H, pi_RR)

    # (1) tau = I frame.
    delta = pi_LL - pi_RR.conjugate()
    o = odd(elim)
    checks["bridge_tau_id.LL_is_half_defect"] = is_zero(o[0:2, 0:2] - delta / 2)
    checks["bridge_tau_id.RR_is_minus_half_conj_defect"] = is_zero(o[2:4, 2:4] + delta.conjugate() / 2)

    # (2) General tau: real orthogonal, so tau conj(tau)^{-1} = 1 and Sigma_tau is an involution.
    tau = sp.Matrix([[sp.Rational(3, 5), -sp.Rational(4, 5)], [sp.Rational(4, 5), sp.Rational(3, 5)]])
    checks["bridge_tau.Sigma_tau_involution"] = is_zero(sigma_tau(tau) * sigma_tau(tau) - sp.eye(4))
    checks["bridge_tau.LL_is_half_transported_defect"] = is_zero(
        odd(elim, tau)[0:2, 0:2] - (pi_LL - tau * pi_RR.conjugate() * tau.inv()) / 2
    )

    # (3) Equivariant point.
    elim_eq = block(pi_RR.conjugate(), pi_LR, pi_LR.H, pi_RR)
    checks["equivariant_point.odd_LL_zero"] = is_zero(odd(elim_eq)[0:2, 0:2])
    checks["equivariant_point.odd_RR_zero"] = is_zero(odd(elim_eq)[2:4, 2:4])

    # (4) Chiral block products with Pi_S = 1.
    Dm = sp.Matrix([[2, sp.I], [0, 1]])
    Dp = Dm.H
    D = block(Z2, Dm, Dp, Z2)
    E = -D * elim * D
    checks["block_products.E_LL"] = is_zero(E[0:2, 0:2] + Dm * pi_RR * Dp)
    checks["block_products.E_RR"] = is_zero(E[2:4, 2:4] + Dp * pi_LL * Dm)

    # (5) Rate identity on a genuine Grassmann tangent.
    v = sp.Matrix([3, 4 + sp.I])
    q = (v * v.H) / (v.H * v)[0]
    Q0 = block(q.conjugate(), Z2, Z2, q)  # J-even orthogonal projector for tau = I
    checks["tangent.Q0_even_projector"] = is_zero(Q0 * Q0 - Q0) and is_zero(parity(Q0) - Q0)
    a = sp.Matrix([[sp.I, 2 + sp.I], [-2 + sp.I, -3 * sp.I]])  # anti-Hermitian
    b = sp.Matrix([[1 - sp.I, 2], [2, sp.I / 2]])  # symmetric
    A = block(a, b, -b.H, -a.conjugate())
    checks["tangent.A_antihermitian_odd"] = is_zero(A + A.H) and is_zero(parity(A) + A)
    U = sp.eye(4) + s * A + s**2 * A * A / 2
    Qs = U * Q0 * U.H
    Q1 = Qs.applyfunc(lambda x: sp.expand(x).coeff(s, 1))
    checks["tangent.Q1_is_comm"] = is_zero(Q1 - comm(A, Q0))
    checks["tangent.Q1_odd_offdiagonal"] = (
        is_zero(parity(Q1) + Q1) and is_zero(Q0 * Q1 * Q0) and is_zero((sp.eye(4) - Q0) * Q1 * (sp.eye(4) - Q0))
    )
    checks["tangent.Q1_nonzero"] = not is_zero(Q1)
    ddelta = (Qs[0:2, 0:2] - Qs[2:4, 2:4].conjugate()).applyfunc(lambda x: sp.expand(x).coeff(s, 1))
    checks["tangent.LL_block_is_half_defect_rate"] = is_zero(odd(Q1)[0:2, 0:2] - ddelta / 2)
    checks["tangent.defect_rate_nonzero"] = not is_zero(ddelta)

    for name, ok in checks.items():
        print(f"  [{'OK ' if ok else 'FAIL'}] {name}")
    failed = [k for k, ok in checks.items() if not ok]
    print()
    print("d Delta_chi/ds at 0 =", list(ddelta.applyfunc(sp.simplify)))
    print(f"{len(checks)} checks, {len(failed)} failed")
    assert not failed, failed
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
