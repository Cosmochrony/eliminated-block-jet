# Eliminated-Block-Jet — The Constrained Jet of the Lorentzian Eliminated Block and the Chiral Generator

*From the Schur Reduction to the Transported Equivariance Defect*

J. Beau, Independent Researcher, France

## Status

Working note, v1.1. Not yet deposited. (v1.1 adds the chiral-generator identification, Section 6.)

## Abstract

This note isolates a single operator-level lemma in the fermionic-matter sub-programme, logically
between the **Projective Residue Schur reduction** (PRS) and the physical identification of the
chiral generator. The Schur reduction writes the projective endomorphism as
$E_\Pi = -\Pi_S D (1-P) D \Pi_S^*$, so the eliminated block $1-P$ controls the three-generation
split coefficient $u$; the companion genus analysis fixes the electric sign $\mu_\chi^2 < 0$, while
$|u|$, the Yukawa sector, and the mass spectrum stay downstream.

The note freezes the structural form the Lorentzian eliminated block must take. Writing the
self-adjoint projector family as a jet $1-P(s) = Q_0 + s Q_1 + s^2 Q_2 + O(s^3)$ along the
$J_\Pi$-odd modulus $s$, the projector constraints force $Q_1$ to be a purely off-diagonal Grassmann
tangent $Q_1 = [A, Q_0]$ and pin the diagonal blocks of $Q_2$ to $Q_1^2$ with opposite signs.

## Results

1. **Constrained-jet theorem.** $Q_0 Q_1 Q_0 = 0 = (1-Q_0) Q_1 (1-Q_0)$, $Q_1 = [A, Q_0]$; and
   $Q_0 Q_2 Q_0 = - Q_0 Q_1^2 Q_0$, $(1-Q_0) Q_2 (1-Q_0) = + (1-Q_0) Q_1^2 (1-Q_0)$.

2. **Chiral source.** With $A = A_+ + A_-$ by $J_\Pi$-parity, the split-sourcing part is
   $\Pi_{J_\Pi\text{-odd}} Q_1 = [A_-, Q_0]$; the split is dead iff $[A_-, Q_0] = 0$, not merely iff
   $Q_1 = 0$.

3. **Split source.** Parity propagation grades the residue jet $E_0$ even, $E_1$ odd, $E_2$ even, so
   $u(s) = s\,\langle J_3, \{E_0, E_1\}\rangle / \langle J_3, J_3\rangle + O(s^3)$; with phase-free
   data the mixing $v$ vanishes identically while $u$ is generically non-zero. Hence $u \neq 0$ is
   admissible with no $R_{\mathrm{mix}}$ leakage and with $|u|$ left free.

4. **Chiral generator (Section 6).** Because the antiunitary Born–Infeld parity exchanges chiralities,
   $[\Pi_{J_\Pi\text{-odd}}(1-P)]_{LL} = \tfrac12 \Delta_\chi(P)$ and
   $[\Pi_{J_\Pi\text{-odd}}\dot Q(0)]_{LL} = \tfrac12 \partial_s \Delta_\chi(P)|_0$. The generator $A_-$
   is therefore not a new selector but exactly the infinitesimal transported chiral equivariance defect
   rate. Its absolute value is not normalised here, leaving a single downstream target: the admissible
   normalisation of $\partial_s \Delta_\chi(P)|_0$.

## Reproducibility

All identities are verified by exact symbolic / exact-rational computation, no sampling, in
`simulation/fermionic-matter/eliminated_block_jet.py` (jet, Sections 3–5) and
`simulation/fermionic-matter/front_d0_generator.py` (chiral generator, Section 6).

## Build

```bash
bash compile.sh   # pdflatex -> bibtex -> pdflatex x2, output in out/
```
