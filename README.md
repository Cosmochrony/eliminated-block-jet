# Eliminated-Block-Jet — The Constrained Jet of the Lorentzian Eliminated Block and the Chiral Generator

*From the Schur Reduction to the Transported Equivariance Defect*

J. Beau, Independent Researcher, France

## Status

Working note. Concept DOI: [10.5281/zenodo.20763532](https://doi.org/10.5281/zenodo.20763532).

## Abstract

This note isolates an operator-level lemma in the fermionic-matter sub-programme, logically between the
**Projective Residue Schur reduction** (PRS) and the physical identification of the chiral generator. The
Schur reduction writes the projective endomorphism as $E_\Pi = -\Pi_S D (1-P) D \Pi_S^*$, so the eliminated
block $1-P$ controls the three-generation split coefficient $u$; the companion genus analysis fixes the
electric sign $\mu_\chi^2 < 0$, while $|u|$, the Yukawa sector, and the mass spectrum stay downstream.

Writing the self-adjoint projector family as a jet $1-P(s) = Q_0 + s Q_1 + s^2 Q_2 + O(s^3)$ along the
modulus $s$, the projector constraints alone force $Q_1 = [A, Q_0]$ to be a purely off-diagonal Grassmann
tangent and pin the diagonal blocks of $Q_2$ to $Q_1^2$ with opposite signs, for every smooth family and
every smooth unitary transport. The parity grading of the jet holds when the family is equivariant under the
antiunitary Born–Infeld parity $J_\Pi$. The contract is realisable: an explicit covariant family
with a moving embedding has $E_0^2 = \mathrm{diag}(1, \tfrac12, \tfrac12)$ and a non-zero split rate.

## Results

1. **Constrained-jet theorem** (every smooth family, general transport $U(s) = 1 + sA + \tfrac{s^2}{2}(A^2+B)$).
   $Q_0 Q_1 Q_0 = 0 = (1-Q_0) Q_1 (1-Q_0)$, $Q_1 = [A, Q_0]$;
   $Q_0 Q_2 Q_0 = - Q_0 Q_1^2 Q_0$, $(1-Q_0) Q_2 (1-Q_0) = + (1-Q_0) Q_1^2 (1-Q_0)$,
   $Q_2 = \tfrac12 [B, Q_0] + \tfrac12 [A, [A, Q_0]]$. A $J_\Pi$-equivariant family,
   $J_\Pi Q(s) J_\Pi^{-1} = Q(-s)$, has alternating jet parity $(+,-,+)$, and for a real-analytic family
   the converse holds; then the odd tangent is the chiral commutator $[A_-, Q_0]$. With a moving embedding,
   a vanishing odd tangent does not by itself kill the first-order split.

2. **Parity propagation with a moving embedding.** Since $P(s) = \Pi_S(s)^* \Pi_S(s)$, the residue jet
   carries the derivatives of the embedding:
   $E_1 = -(G_1 Q_0 G_0^* + G_0 Q_1 G_0^* + G_0 Q_0 G_1^*)$ with $G = \Pi_S D$. For a $J_\Pi$-even $D$
   and a covariant embedding, $J_\Pi \Pi_S(s) J_\Pi^{-1} = \Pi_S(-s)$, the jet is graded $E_0$ even,
   $E_1$ odd. Freezing the embedding is an auxiliary model with a different rate.

3. **Split source and vanishing of mixing.** The first-order $J_\Pi$-odd part of $E_\Pi^2$ is
   $\{E_0, E_1\}$. A Hermitian operator odd under the antiunitary parity has zero $(+,-)$ entry on
   $\mathbb{C}^3_{\mathrm{gen}}$, so generation mixing vanishes at first order in every covariant family,
   with or without complex phases. Mixing at second order lies in the even sector and is not excluded.
   In the normal form $E_0 = -\mathrm{diag}(1, c, c)$, $c = 2^{-1/2}$, the split rate is
   $u'(0) = -2c\,(E_1)_{++}$.

4. **Realisation.** An explicit covariant family on $\mathbb{C}^3 \oplus \mathbb{C}^3$ with a moving
   embedding has $E_0^2 = \mathrm{diag}(1, \tfrac12, \tfrac12)$, $u'(0) = \tfrac{11}{9}\,2^{1/4} \neq 0$,
   neither mixing nor singlet leakage at first order, non-zero mixing at second order, and frozen-embedding
   rate $\tfrac{79}{90}\,2^{1/4}$.
   This is an existence statement under explicit hypotheses; the programme is not shown to select this family.

5. **Chiral-diagonal defect-rate identity.** Because the antiunitary Born–Infeld parity exchanges
   chiralities, $[\Pi_{J_\Pi\text{-odd}}(1-P)]_{LL} = \tfrac12 \Delta_\chi(P)$ with
   $\Delta_\chi(P) = \pi_{LL} - \tau\,\overline{\pi_{RR}}\,\tau^{-1}$, and
   $[\Pi_{J_\Pi\text{-odd}}\dot Q(0)]_{LL} = \tfrac12 \partial_s \Delta_\chi(P)|_0$. This determines one block
   of the odd tangent, not a unique generator: $A$ is defined only modulo the centraliser of $Q_0$.

**Open**: the selection of a covariant family by the programme; the uniqueness of the split carrier and the
injective localisation of the odd tangent onto its $LL$ block; the normalisation of
$\partial_s \Delta_\chi(P)|_0$, hence $|u|$; generation mixing at second and higher order.

## Reproducibility

All identities, counterexamples, and the realisation are verified by exact symbolic and exact algebraic
computation (SymPy, no sampling). From a fresh clone:

```bash
pip install -r code/requirements.txt
python code/eliminated_block_jet.py    # jet, propagation, split source, realisation (93 exact checks)
python code/front_d0_generator.py      # defect bridge and rate identity (15 exact checks)
```

## Build

```bash
bash compile.sh   # pdflatex -> bibtex -> pdflatex x2, output in out/
```
