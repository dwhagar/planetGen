# The magma-ocean redox model (GEN.1005, part 1)

Status: module and registry rows, written 2026-10-30 (as the lane reported) by the Aux 3 lane. `src/planetgen/physics/magma_ocean_redox.py`; tests in `src/tests/test_magma_ocean_redox.py`; registry rows MR1 to MR7 (section 6.57 of the gap-filling plan). Part 1 changed nothing the generator draws. Part 2 (the digest flip, 2026-10-30, with the dry-carbon and comet-species rows) has rows a1 and a2 call `magma_ocean_diw`, `magma_ocean_diw_sigma` and `later_mantle_diw`; SSR1 and SSR2 are retired and the two rows are `partly verified`.

## 1. What replaces what

Rows a1 and a2 drew the redox of a rocky body's outgassing mantle from lines in the log of the mass fitted to the Moon, Mars and the Earth (`magma-ocean-redox-against-mass.md`, GEN.482). GEN.1005 asks for a model that carries the physics those three points sample, so that no Solar System planet is an input.

The model has one input, the mass of the body, and four steps:

1. **Base pressure.** `P = 55 GPa x M^0.55`, held at 140 GPa (MR1, MR2). 55 GPa is [R1075]'s base pressure of a one-Earth-mass ocean; the power is inferred from their three bases (5, 14 and 55 GPa) and gives 4.9 and 16.1 GPa at the first two. The 140 GPa is where their first-principles calculations stop, at 5.4 Earth masses.
2. **Ferric share.** Fe3+/Fe2+ at the ocean's surface is the reference value times (P / 55 GPa)^0.75 (MR3, MR4). Disproportionation of FeO into Fe and Fe2O3 favours the smaller molar volume, which is Fe2O3 up to 140 GPa [R1075]. The reference Fe3+/FeT is 0.068, the middle of [R1951]'s 0.034 to 0.103 for the Earth. The power is the middle of 0.5 (the slope of [R1951]'s ocean-base conditions) and 1.0, which holds Deng's 0.88 (their text: the ferric share halves from 55 to 25 GPa).
3. **From ratio to fugacity.** `dIW = 4 log10(Fe3+/Fe2+) + 5.62` (MR5). The factor 4 is the stoichiometry of 4 FeO + O2 = 2 Fe2O3. The 5.62 is read off two pairs of [R1951] for the Earth's ocean surface (Fe3+/FeT 0.034 gives IW-0.2, 0.103 gives IW+1.8; the relation gives -0.19 and +1.86). The temperature term of Hirschmann's equation is not in it.
4. **Width.** The standard deviation of the two uniform ranges (the reference share, the power), in quadrature, times four (MR6): 0.6 at the reference pressure, 0.85 at 4.9 GPa.

The later upper mantle (row a2) is the ocean surface of the body (at most the one-Earth-mass value) plus an increment of 1.2 at 0.1 Earth masses or less up to 2.4 at one (MR7), with width 1.0.

## 2. What it reproduces (tests, never inputs)

| Body (test point) | Base pressure, model / published | Surface dIW, model | Published surface range | Later mantle, model | Published later mantle |
|---|---|---|---|---|---|
| Moon (0.0123 Earth masses) | 4.9 / 5 GPa | -2.08 | -2.2 to -1.5 [R1951]; IW-2 [R1075] | -0.88 | IW to IW-2 |
| Mars (0.107) | 16.1 / 14 GPa | -0.53 | -1.4 to +0.2 (Hirschmann); IW-0.3 (Deng) | +0.71 | IW to IW+1 |
| Earth (1) | 55 / 55 GPa | +1.07 | -1.2 to +1.8 (Hirschmann); IW+2 (Deng) | +3.47 | IW+3.5 |

The model's surface values sit within 0.25 of the middle of [R1951]'s range for the Moon and Mars and in the upper half for the Earth (the middle is +0.3); the old line was within 0.2 of every middle. The difference is the price of a single physical form for all three. Above one Earth mass the surface keeps rising (+1.9 at three Earth masses, +2.3 at 5.4) and then stops: no source reaches further. The old line rose 1.0 per decade of mass to the same point and reached +1.2.

## 3. What was not read

`par.nsf.gov` refuses automated reads now (robots), so Eq. 21 and Eq. 34 of [R1951] and the equation bodies of [R1075] (the Nature Communications text on the publisher's page has the equation numbers but not the equations, and the supplement was not reachable) were not available this session. The model uses the numbers GEN.482 recorded when the paper was read (the ranges of section 3.1 of `magma-ocean-redox-against-mass.md`, the printed pairs) and the sentences of the Deng text that are readable (factor of two from 55 to 25 GPa; bases of 55, 14 and 5 GPa; calculations to 140 GPa). When the equations are read, MR3 to MR5 become derivations and their tiers and tweaks move.

## 4. What it does not do

- It is not an interior model: the base pressure is a power of the mass. A body with a shallower or deeper ocean (hotter, wetter, younger) gets the wrong pressure.
- It carries no temperature (the surface temperature is one value) and no composition (Cr and Ti effects, [R1951] sections 6.4 and 6.5).
- The graphite-flotation alternative ([R2329], down to IW-5) is another mechanism and is not drawn.
- The later mantle is a line between two numbers, not a model of the oxidation (hydrogen loss, lower-mantle disproportionation, recycling).
