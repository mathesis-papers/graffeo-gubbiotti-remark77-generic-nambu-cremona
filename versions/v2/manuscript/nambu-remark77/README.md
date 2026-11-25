# Graffeo–Gubbiotti, Remark 7.7 (arXiv:2301.09129v2, pp. 35–37)

## Nambu conjugacy and the twelve-plane configuration

Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang (National University of Singapore; ChronoAI Pte Ltd, Singapore).

License: Creative Commons Attribution 4.0 International (CC BY 4.0).

The note proves the conjecture in the exact posing version, Remark 7.7 of [arXiv:2301.09129v2](https://arxiv.org/abs/2301.09129v2), for a positive-definite real symmetric first matrix, a real symmetric second matrix with simple generalized spectrum, and every nonzero real step. An explicit complex projective change sends the homogeneous Nambu Kahan map to the cubic Cremona map followed by `2 I_4 - 11^T`.

The same change describes the eight base points, the six reducible net members and all twelve planes from (7.16), (7.20)–(7.23). The eight points are projective sign vectors; the planes are `v_i ± v_j = 0`. Each point belongs to six planes, each plane contains four points, and the conjugator sends the point set to the coordinate points together with their images under `2 I_4 - 11^T`.

The written proofs also extend the conjugacy to complex symmetric pencils whose binary determinant has three distinct projective roots. For real positive-definite first matrix, every eigenvalue multiplicity is covered: the iterate degrees are `2n^2+1` for simple spectrum, `2n+1` for exactly two distinct eigenvalues, and `1` for scalar pencils. The first degree law is transported from Graffeo–Gubbiotti, Proposition 5.3. The integrals defining the net are credited to Celledoni et al., Proposition 5. Singular pencils, nonsemisimple collisions and the separate theta map in (8.5)–(8.7) are outside this coverage.

## Build

With Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is inside `note.tex`; BibTeX and external figures are unnecessary.

## Exact identities and finite enumeration

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python verify.py --output verification.json
```

The checker uses SymPy 1.14.0 and rational symbolic arithmetic. Its 23 identity groups have 107 zero polynomial residuals. They include the universal Euler cubic, the conjugator determinant, cross-product covariance, reversal and composition constants, six reducible net members and a rational nonorthogonal example with two noncommuting positive-definite forms. The sign-point enumeration checks all 48 incidences and all eight images in the normal form. Four explicit collision iterates have degrees 3, 5, 7 and 9. The spectral arguments, general collision degree law and exhaustive classification of reducible quadrics have written proofs in the note.

## Lean declaration coverage

The [pinned Lean file](https://github.com/the-omega-institute/trureturing/blob/86172888e706177b6c1d135cb12aec3b6d2cdce8/D5/S3/FluidDynamics/Solitons/NambuKahanCremonaDecomposition.lean) contains `D5.S3.FluidDynamics.Solitons.NambuKahanCremonaDecomposition.result`. Its proposition is the existence of a two-projectivity factorization for two real positive-definite matrices, the rank condition, and positive step. The auxiliary `euler_top_identity` proves the complex Euler cubic with specified square roots; `covariance` proves the coordinate-change identity for a given simultaneous congruence. The written note proves the broader conjugacy, incidence and collision statements. The matrix conventions agree since the Jacobian of `4(Ax × Bx)` is four times the source's `nambuLinearPart`.

## Contents

- `note.tex`, `note.pdf`: exact statements and complete focused proofs.
- `metadata.json`: stable scientific identity, creators, license and mathematical coverage.
- `verify.py`, `requirements.txt`, `verification.json`: symbolic identities and finite enumeration.
- `coverage.json`: posing-source locators, statuses and proof scope.
- `MANIFEST.sha256`: member hashes.
