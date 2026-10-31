# Stability of a 3+1 quadruple star system (GEN.1029)

Part of GEN.1006 (overlap-free multi-star systems, [overlap-free-systems.md](overlap-free-systems.md) section 5 and 13). This note says which criterion the generator uses for a 3+1 quadruple, how it was ported, how accurate it is, and what it does not cover. The module is `physics/quadruple_stability.py`; GEN.1016 (Core 3) wires it into the multi-star limits. Stand-in rows QS1 to QS4 in section 6.80 of [gap-filling-plan.md](../research/gap-filling-plan.md).

## 1. The problem

A 3+1 quadruple is an inner binary `(m1 >= m2)`, a third star `m3` on the intermediate orbit around the binary and a fourth star `m4` on the outer orbit. `physics/stability.py` tests a triple at each node (Mardling and Aarseth 2001, Eggleton and Kiseleva 1995, Vynatheya and others 2022 Equation 4). Applied to a quadruple as two nested triples that is right for 66 percent of 3+1 systems (88 percent of 2+2) and calls too many of them stable (precision 0.59): the outer star excites the eccentricity of the intermediate orbit, which destabilises the inner triple, and a triple test cannot see that (Vynatheya, Mardling and Hamers 2023, arXiv 2301.09930, Tables 3 and 4).

## 2. The form chosen: the authors' network as plain numbers

The paper's 3+1 classifier is a multi-layer perceptron: 11 inputs, four hidden layers of 50 logistic units, one logistic output (the probability of "unstable"), trained on 5e5 N-body integrations. The authors publish it as a scikit-learn model file (`mlp_model_3p1_ghost_v1.2.2.pkl`, 209 KB, github.com/pavanvyn/quadruple-stability). Three ways to use it were compared:

| Way | Cost | Verdict |
|---|---|---|
| scikit-learn at run time (the authors' script) | a heavy dependency, a pickle loaded at run time (a pickle can run code), version warnings | no |
| a fitted formula or a table of the network's answers | loses accuracy; a grid over 11 dimensions is impossible | no |
| the weights as numbers and a numpy forward pass | 8,301 numbers (111 KB of JSON), five matrix products, no new dependency | **chosen** |

The network is a plain feed-forward net, so the port is exact, not a fit. `scripts/build_quadruple_stability.py` reads the pickle with an unpickler that refuses every global except numpy's array constructors and keeps only the state of every scikit-learn object, so the pickle's code never runs; it refuses any network that is not 11-50-50-50-50-1 logistic. The forward pass equals scikit-learn's `predict_proba` to 1e-8 on five stored cases and 2000 random inputs (rounded to nine significant digits), and the test file carries the stored cases. No machine-learning library is needed to run it; scikit-learn is not a dependency of the project.

Inputs (all dimensionless): `q_in = m2/m1`, `q_mid = m3/(m1+m2)`, `q_out = m4/(m1+m2+m3)`, `alpha_in_mid = a_in/a_mid`, `alpha_mid_out = a_mid/a_out`, the three eccentricities and the three mutual inclinations (radians, divided by pi inside). The verdict is "stable" when the probability of "unstable" is below 0.5, as in the authors' script. `ratios()` makes the five ratios from masses and axes; `outside_range()` names inputs outside the documented range; `hierarchy_holds()` is the apoapsis-inside-periapsis test the training sample kept.

## 3. Accuracy

**The authors' figures** (their table 4, 20 percent test set of 1e5 systems, labelled by N-body runs with their ghost criterion): overall 0.93; precision 0.91 stable and 0.95 unstable; recall 0.91 stable and 0.95 unstable. Nested triples (MA01): 0.56 overall, precision 0.54 stable. These are the authors' numbers; their data set is shared "upon reasonable request" and was not obtained, so I could not re-score their test set.

**Reproduced here:** the numpy network gives the scikit-learn probabilities (test: five stored cases to 1e-6, and the builder checks the shape of the net). That is the whole claim about the port; it says nothing about the network's accuracy.

**Independent N-body check (my sample, this lane):** REBOUND 5.2.2 (IAS15), the paper's own labelling rule (100 outer orbits, a twin with the inner axis larger, unstable if the pair diverges by more than 1 percent in the inner axis or any orbit unbinds), masses log-uniform over a decade, isotropic angles, systems near the Mardling-Aarseth boundary. The results are below.

### 3.1 Result of my sample

Two samples, 331 systems that finished (scripts and results: `/mnt/project-files/data/gen1029-nbody/`, `nbody.py`, `analyse.py`). Each system ran twice (original and twin, relative offset 1e-8 of the inner axis, outputs four times per outer orbit), 100 outer orbits; systems that did not finish in 90 s were dropped (none in the second sample).

| Sample | Systems | N-body stable | Network right | Nested MA01 right | Network false stable | Network false unstable |
|---|---|---|---|---|---|---|
| wide: alpha 0.18 to 0.45, margin 0.6 to 5 times the MA01 boundary, 20 percent of MA01-unstable kept | 228 | 5 | 97.8 percent | 61.0 percent | 5 | 0 |
| near the boundary: alpha 0.10 to 0.30, margin 0.8 to 3, 50 percent of MA01-unstable kept | 103 | 14 | 89.3 percent (90 percent range 84.5 to 94.2) | 23.3 percent | 10 | 1 |
| both | 331 | 19 | 95.2 percent | 49.2 percent | 15 | 1 |

What it shows and does not show:

1. **The network is far better than the nested test on the same systems** (95 against 49 percent over both samples; the nested test calls 78 of the 103 near-boundary systems stable that the N-body run breaks). This agrees with the paper's finding that nested triples overestimate the stable systems.
2. **Its accuracy matches the paper's 93 percent only roughly, and its precision for "stable" is lower.** On the near-boundary sample 13 of 14 truly stable systems were found (recall 0.93) but 10 of 23 predicted stable broke (precision 0.57, the paper has 0.91). The sample has few stable systems (19 of 331) against about half in the authors' set, so the stable-class numbers have wide error; all-unstable guessing scores 86 percent on the near-boundary sample, so the score alone proves little and the false-stable count matters more.
3. **My labels are not the authors' labels.** Their code (MSTAR, a regularised integrator), the size of the twin offset (1e-6 au on an unknown inner axis), the output cadence, the 5-hour cutoff and their exact sampling of the boundary are not all stated; the divergence rule is sensitive to the offset (phase drift of a regular system alone gave a divergence of 0.6 percent at a relative offset of 1e-6, which is why I used 1e-8). Treat 89 to 95 percent as agreement within the uncertainty of the method, not a re-measurement of 93.
4. **Conclusion for the generator:** the network is the right tool for a 3+1 node; use it as the verdict after the nested-triple screen, and treat a "stable" verdict as likely rather than certain (a quarter to a half of near-boundary "stable" calls may break within 100 outer orbits by this check), so keep a margin on the axis ratio of systems it calls stable. The margin is a stand-in for GEN.1016 to register.

## 4. What it does not cover

1. **Stable means 100 outer orbits and not chaotic.** It is a short-span chaos test, not a guarantee for a stellar lifetime. Slow secular evolution (Lidov-Kozai cycles longer than 100 outer orbits) is not in the label.
2. **Training range.** Masses within a factor 10, eccentricities to 0.95, isotropic angles, systems kept near the Mardling-Aarseth boundary of the two nested triples (20 percent of the MA01-unstable draws kept). Far from the boundary (very wide or very tight) the network extrapolates; `outside_range` flags mass ratios above 3 (`q_mid`) and 2 (`q_out`) as the authors' script advises. The axis ratios are not given a range by the paper.
3. **Only 3+1.** The 2+2 network is a separate module function set (`*_2p2`, GEN.1031, [quadruple-stability-2p2.md](quadruple-stability-2p2.md)).
4. **No arguments of periapsis or mean anomalies** are inputs, so two systems that differ only in them get the same verdict (the paper says they can matter).
5. **Licence.** The authors' repository has no licence file (the raw `LICENSE` returned 404 on 10-30). The weights are numbers from a published, citable model; Boss handles DATA-LICENSE.md. Cite Vynatheya, Mardling and Hamers 2023 wherever the rule is used.

## 5. For the wiring (GEN.1016)

Call `quadruple_stability.is_stable(*ratios(...), e_in, e_mid, e_out, i_in_mid, i_in_out, i_mid_out)` on every 3+1 node of a tree after the nested-triple screen; treat `outside_range(...)` non-empty as "lower confidence" (the verifier of overlap-free-systems.md section 4.7 flags it) and keep the nested-triple test with a margin as the first screen. Use a different criterion in the verifier from the generator: EK95 for the triples, this network for the 3+1 node in the generator only.

## 6. As built: GEN.1016 part 1 and GEN.1028 (Core 3, 10-30)

**What is in the code now, lock-free** (nothing the generator calls changes; `physics/hierarchy_stability.py` is read by tests only until the wiring merges on its train).

1. **The Mardling and Aarseth range (GEN.1028).** `stability.MA01_Q_OUT_MAX` is 5, the paper's own limit (registry row GM12). `stability.criterion` picks Equation 4 inside its fit, Mardling and Aarseth where `q_out <= 5` holds, and Eggleton and Kiseleva beyond that; a `required_axis_ratio` asked for `"ma01"` by name past 5 raises. Generation almost never reached the old case (150 systems: 1 of 294 triple tests used Mardling and Aarseth at all, none past 5); the case it did change is a satellite pair about a star in the hosting gate (`hosting_gate.stable_orbit_exists`), where Equation 4 stops at `q_out` 100 and the old fallback asked the Earth and the Moon for an axis ratio of 430 against the real 390.
2. **Triple levels with a margin.** `hierarchy_stability.verdicts(nodes, members)` judges every node that has a sub-node for a child as a triple at `level_margin` times the criterion's axis ratio (1.0, none, row GM8, because the nested test is right for 88 percent of 2+2 systems).
3. **The 3+1 quadruple.** Both of its triple levels must pass, then the network must call the system stable with `alpha_in_mid` and `alpha_mid_out` both multiplied by `axis_margin` 1.3 (row GM13). The 1.3 is read off the 331 REBOUND systems of section 3 (`/mnt/project-files/data/gen1029-nbody/`): the plain network was right on 18 of its 33 stable calls (0.55), at 1.1 on 14 of 19 (0.74), at 1.3 on 8 of 9 (0.89, the paper reports 0.91) keeping 0.42 of the truly stable systems, at 1.5 on 5 of 5 keeping 0.26. Outside what the authors documented (`outside_range`, or apsides that overlap) both triples must clear `fallback_margin` 1.5 (row GM7) and the verdict says `confident=False`.
4. **The 2+2 quadruple** (GEN.1031 ported its network): both nested triples (each pair with the other pair as the third body) must pass at the level margin, then the 2+2 network must call it stable with `alpha_in1_out` and `alpha_in2_out` both multiplied by its `axis_margin` 1.3 (row GM14). That 1.3 is the 3+1 margin carried over: there is no 2+2 N-body sample of our own behind it (tweak L4), and an own sample replaces it. Outside the network's range or with overlapping apsides the fallback is the same as for 3+1 (both triples at 1.5, `confident=False`).

**What the generator draws today** (`/mnt/project-files/reports/gen1016-quadruple-need-2026-10-30.md`, 300 seeds of `STARS = 4`, no planets): 166 trees are 3+1 and 134 are 2+2. Every triple level passes, as the generator promises. Of the 166 3+1 trees, none is outside the network's range, 68 (41 percent) are called unstable by the plain network and 91 (55 percent) with the margin. The middle orbit is eccentric (median `e_mid` 0.33) and the outer one more so (median `e_out` 0.65), which is the case the nested test cannot see. So the wiring has to enforce the 3+1 verdict while the system is drawn, or most 3+1 systems of the galaxy would be tight enough to be disrupted in the network's 100 outer orbits. With the 2+2 network the 134 2+2 trees give 5 unstable without the margin (4 percent) and 16 with it (12 percent): the nested test said stable for all of them.
