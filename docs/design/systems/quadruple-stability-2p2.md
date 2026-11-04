# Stability of a 2+2 quadruple star system (GEN.1031)

Part of GEN.1006 (overlap-free multi-star systems, [overlap-free-systems.md](overlap-free-systems.md)). This is the 2+2 counterpart of [quadruple-stability-3p1.md](quadruple-stability-3p1.md) (GEN.1029); how the network is ported, why a plain numpy forward pass, and the reading of "stable" are the same and are not repeated. The module is `physics/quadruple_stability.py` (functions with the suffix `_2p2`); GEN.1016 (Core 3) wires it into the multi-star limits. Stand-in rows QT1 to QT4 in section 6.81 of [gap-filling-plan.md](../research/gap-filling-plan.md).

## 1. The problem and the form chosen

A 2+2 quadruple is two binaries, `(m1 >= m2)` and `(m3 >= m4)` with `m1 + m2 >= m3 + m4`, that orbit each other. `physics/stability.py` can test it as two nested triples (each binary with the other as one point mass) and that is right for 83 percent of the paper's systems ([R4575: Table 3], Mardling-Aarseth), 88 percent with their triple network. The paper's 2+2 network is right for 94 percent. It is the same kind of network as the 3+1 one: 11 inputs, four hidden layers of 50 logistic units, one output (the probability of "unstable"), 8,301 numbers. `scripts/build_quadruple_stability.py 2p2 <model>` makes `data/quadruple_stability_2p2_weights.json` from the authors' `mlp_model_2p2_ghost_v1.2.2.pkl` without running the pickle (the restricted reader of GEN.1029); the numpy forward pass gives the probabilities of the same network built as a scikit-learn model from that state to 1e-6 on five stored cases (the test file carries them).

Inputs (dimensionless): `q_in1 = m2/m1`, `q_in2 = m4/m3`, `q_out = (m3+m4)/(m1+m2)`, `alpha_in1_out = a_in1/a_out`, `alpha_in2_out = a_in2/a_out`, the three eccentricities, the three mutual inclinations (radians, over pi inside). `ratios_2p2` makes the five ratios from masses and axes, `outside_range_2p2` names inputs outside what the paper documents (mass ratios under 0.1, `q_out` above 1, eccentricity above 0.95, an axis ratio outside 0 to 1), `hierarchy_holds_2p2` is the apoapsis-inside-periapsis test the training sample kept for both inner binaries.

## 2. Accuracy

**The authors' figures** (their Table 3, 20 percent test set of 1e5 systems): overall 0.94; precision 0.94 stable and 0.95 unstable; recall 0.94 and 0.95. Nested Mardling-Aarseth: 0.83 overall, precision 0.77 stable.

**Independent N-body check (this lane):** REBOUND 5.2.2 (IAS15), the paper's labelling rule (100 outer orbits, a twin with the inner axis of the lighter binary larger by a relative 1e-8, unstable if a binary unbinds, the outer orbit unbinds, or the twin's inner axis differs by more than 1 percent), masses log-uniform over a decade, isotropic angles, eccentricities uniform to 0.4 (to 0.95 for a fifth of the draws), the heavier binary first. Scripts and results: `/mnt/project-files/data/gen1031-nbody/`. Two samples, none dropped:

| Sample | Systems | N-body stable | Network right | Nested MA01 right | All "unstable" right | Network precision stable | Network recall stable |
|---|---|---|---|---|---|---|---|
| paper-like: axis ratios uniform 0.03 to 0.6, 30 percent of MA01-unstable draws kept (the paper keeps 30) | 185 | 26 | 91.4 percent (90 percent range 87.3 to 94.2) | 81.6 percent | 85.9 percent | 0.62 (16 false stable of 42) | 1.00 |
| near the boundary: axis ratios 0.06 to 0.35, MA01 margin 0.8 to 3, 50 percent of MA01-unstable kept | 787 | 124 | 68.6 percent | 39.5 percent | 84.2 percent | 0.33 (246 false stable of 369) | 0.99 |

What it shows and does not show:

1. **On the paper's kind of sample the result is close to the paper's.** 91 percent for the network against the paper's 94, and 82 percent for the nested test against the paper's 83 (the second number matching is the best sign that my sample and labels resemble theirs). The network is better than the nested test, by 10 points here and 11 in the paper.
2. **The network is lenient.** It never missed a truly stable system (one false unstable in 972), but it calls a share of the unstable ones stable: precision for "stable" is 0.62 on the paper-like sample (paper 0.94) and 0.33 near the boundary. This is the same direction as for the 3+1 network (precision 0.57 near the boundary) but stronger.
3. **The false stables are marginal.** Of 30 near-boundary systems that MA01 calls stable (margin above 1.3) and the N-body run broke, which I re-ran to see how, 29 were still bound with the twin's inner axis diverging by mostly 1 to 5 percent, not disrupted; in one the twin unbound. Re-labelling the paper-like sample with a divergence limit of 3 percent (instead of 1) makes the network right in 95 percent of the systems, with 10 percent it is 98 percent. The 1 percent limit is the paper's, so the leniency is partly the weak chaos of systems on the edge of the label, not a failure of the network. Other differences I cannot close: their code (MSTAR), their twin offset (1e-6 au on an unstated inner axis, against my relative 1e-8), their exact sampling.
4. **Few stable systems in my samples** (14 and 16 percent against about half in the paper), so the stable-class scores have wide error and always saying "unstable" scores 84 to 86 percent. The comparison of the network with the nested test on the same systems is the sturdier statement.
5. **Conclusion for the generator:** the network is a better 2+2 test than two nested triples, which are already good (83 percent). Use it after the nested screen, treat "stable" as likely and keep a margin on the axis ratio, as for the 3+1 case. Where the generator needs a hard guarantee of no overlap, the nested test with a margin is the cheaper first screen.

## 3. What it does not cover

1. **Stable means 100 outer orbits and not chaotic**, a short-span chaos test; slow secular evolution is not in the label.
2. **Training range.** Masses within a factor 10 (inner mass ratios at least 0.1), `q_out` at most 1, eccentricities to 0.95, isotropic angles, systems near the Mardling-Aarseth boundary of the two nested triples (the paper excludes 70 percent of the MA01-unstable draws). Far from the boundary the network extrapolates.
3. **No arguments of periapsis or mean anomalies** are inputs.
4. **Only the two quadruple shapes of the paper.** Higher multiples are the nested-triple test at each node.
5. **Licence.** As for the 3+1 file: the authors' repository has no licence file; Boss handles DATA-LICENSE.md. Cite [R4575] wherever the rule is used.

## 4. For the wiring (GEN.1016)

Call `quadruple_stability.is_stable_2p2(*ratios_2p2(...), e_in1, e_in2, e_out, i_in1_in2, i_in1_out, i_in2_out)` on every 2+2 node of a tree after the nested-triple screen, treat a non-empty `outside_range_2p2(...)` as lower confidence, and keep a margin on the axis ratios of systems it calls stable.

## 5. The axis margin (GEN.1155)

`axis_margin` multiplies `alpha_in1_out` and `alpha_in2_out` before the network is asked, so a system is accepted only when it is stable with both inner axes that much wider. The GEN.1016 value of 1.3 was the 3+1 margin carried over. GEN.1155 ran 302 more systems through the same REBOUND label (IAS15, 100 outer orbits, a twin with a relative 1e-8 inner-axis offset, unstable above 1 percent divergence or an unbound orbit; `nbody2p2_margin.py`): masses log-uniform over a decade (inner mass ratios at least 0.1), both inner axis ratios log-uniform from 0.02 to 0.5, eccentricities to 0.95 for a fifth of the draws and to 0.4 for the rest, isotropic angles. Only systems the network calls stable with no margin and that sit inside its documented range were run (the others cannot become stable with a margin); those still stable at a margin of 2.5 were run one in seven and weighted by seven.

| margin | new sample: precision, recall | paper-like (185): precision, recall | near boundary (787): precision, recall |
|---|---|---|---|
| 1.0 | 0.83, 1.00 | 0.62, 1.00 | 0.33, 0.99 |
| 1.3 | 0.95, 0.93 | 0.81, 0.85 | 0.62, 0.77 |
| 1.5 | 0.99, 0.83 | 0.94, 0.58 | 0.80, 0.51 |
| 1.7 | 1.00, 0.73 | 1.00, 0.38 | 0.90, 0.31 |
| 2.0 | 1.00, 0.63 | 1.00, 0.27 | 1.00, 0.12 |

The paper reports 0.94 precision for stable calls. 1.5 is the smallest margin that reaches it on the paper-like sample and is above it on the wider one; the near-boundary set is built from the hardest systems and stays below at every margin under 1.7. The value is tuned to a published statistic (tier T6, tweak level L2).
