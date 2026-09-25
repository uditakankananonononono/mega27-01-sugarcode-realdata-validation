# Numbered, used formula audit - September 25, 2026

Scope: item 1 validation paper in `paper/paper.tex` plus its included `paper/derivations2.tex`. This gate counts **displayed equations that receive a PDF equation number and are tied to an implemented, committed benchmark or result**. It does not count inline formulas, explanatory restatements, software tests, repository/build commands, citations, or a source's mere availability as a research/data tool. The count is separate from the user's external-scientific-tool and direct-dataset gates.

| PDF equation | Formula and use in the study | Source in TeX |
|---:|---|---|
| 1 | Position-weight log-odds and score for splice-site window checks | `paper/paper.tex`, F1 |
| 2 | Rank-sum/AUC identity for held-out discrimination | `paper/paper.tex`, F5 |
| 3 | Class-weighted logistic objective used in the splice CNN baseline | `paper/paper.tex`, F7 |
| 4 | Binding free-energy difference from paired dissociation constants, used for the committed SKEMPI parity check | `paper/paper.tex`, F9 |
| 5 | Stoichiometric flux primal/dual pair, used for the BiGG *E. coli* core feasible optimum/parity check | `paper/paper.tex`, F10 |
| 6 | ACMG point-to-odds posterior conversion, checked against the committed six-point/published 0.900 calculation | `paper/paper.tex`, F12 |
| 7 | Maximum-entropy splice-site model and normalization in the second-round comparison | `paper/derivations2.tex`, F15 |
| 8 | Fisher transform for the reported approximate correlation assessment | `paper/derivations2.tex`, F17 |
| 9 | Shrake-Rupley solvent-accessible surface area used in the structural cross-check | `paper/derivations2.tex`, F21 |
| 10 | Nearest-neighbor melting temperature and salt correction used for oligo comparison | `paper/derivations2.tex`, F22 |

Equation numbers 4-6 were previously inline F9/F10/F12 derivations and are now displayed and numbered, with explicit use and caveats. This is not three new findings or extra scientific services. Rendered `paper/paper.pdf` is 77 pages; equations 4-5 were visually inspected on page 3 and equation 6 on page 4. The literal source count is ten `equation` environments, all rendered with numbers 1-10 in order. A future edit should re-render and recheck before citing the count.
