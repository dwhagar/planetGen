**Headline:** An air with a haze that takes more of the sunlight aloft is never read as warmer at the ground, so a hazy world no longer has several fixed points that depend on where its solve began.

### Changed

- GEN.1780: `greenhouse_table.warming` reads an air with sunlight taken aloft (`aloft_share` over 0) as the running minimum of the network's readings at the shares of the same haze from none up to the share asked (`ALOFT_MONOTONE_NODES` 10 nodes per unit share, interpolated between the two nodes around the share), so the ground is never warmer at a larger share. The table gave 171, 199 and 229 K where the full model gave 166, 155 and 139 K for shares 0, 0.3 and 0.7 of a nitrogen and methane air under a faint star; the rise made a second and third fixed point that the solve found or missed by its start. Stand-in row IR20 (T5, L0).
- The table file itself is unchanged (a refit that is better on every held-out set but does not yet meet the stated-accuracy bounds goes with GEN.1781). `scripts/fit_greenhouse_table.py` gains the `hazy`, `cold` (with `--deep`) and `designed` families, held-out `validate` figures with hydrogen-rich airs stated apart (`hydrogen_rich`) and `--rng-offset`, which the refit needs.
- The haze replay test (`tests/test_haze_in_the_air.py`) no longer skips the world whose replays disagreed (seed 106, 118 K): its replays now converge to one answer.

generation-output: changed
generation-affects: fill
