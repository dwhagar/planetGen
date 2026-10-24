# The design-flag registry (GEN.525)

Status: built 2026-10-23 by lane M (builder B17). `planetgen/physics/design_flags.py`, the data file `planetgen/physics/data/design_flags.json`, tests in `src/tests/test_design_flags.py`.

## 1. What a design flag is

A row of the atmosphere recipes whose quantity no ensemble can define: no exoplanet has a measured value, or the quantity is whether life makes a gas. Tuning such a row to the only anchor there is (the Earth's 400 ppbv of nitrous oxide, Venus's contested 20 ppb of phosphine, the modern Earth's lightning) would fix the Solar System into the generator, which NO FIXED CONTENT forbids. Decision D4 of `tunable-gaps.md` keeps a labelled design bracket, off unless a named flag allows it, never tuned.

## 2. The registered flags

| Row | Stands for | Bracket (constant) | On only if |
|---|---|---|---|
| m3 | nitrous oxide (N2O) | 0.05 to 0.5 log (`c_m3`) | the life stage is chemical or later, or stellar energetic particles drive weakly reducing air |
| m4 | phosphine (PH3) | 0 to 0.05 uniform (`c_m4`) | the life stage is chemical or later |
| m11 | nitrogen dioxide (NO2) | 0 to 0.3 uniform (`c_m11`) | the life stage is chemical or later, or lightning acts |

Rows j2, j3 and j4 (the Mars-, Titan- and Pluto-like atmospheres) were listed with these by the gap plan, but Boss's decision of 2026-10-18 (fork 7) already made them derived outputs of the retention physics; they are not design flags. The body-level items the gap plan also names (the within-system comet scatter of k7, the belt taxon mix E3, the moon origin T17 and T18) are not rows of the atmosphere recipes and belong to the builders of those tables (GEN.520, GEN.521, GEN.522); each registers its flag in the same file.

## 3. What changed

- Rows m3, m4 and m11 were stand-ins naming this item. They are now `partly verified` with the limit "design flag", the honest label for a part that is a design bracket.
- Row m4 was unconditional (up to 5 percent of every atmosphere had phosphine). It is now zero unless the life stage is chemical or later, as D4 says.
- `design_flags.check()` holds the recipes against the registry: each row exists, is not `verified`, holds the registered constant, is not tunable, has no target, reads every enabling flag in its weight and is zero when none is set. A later edit that ungates a row or tunes it fails the test.

## 4. What it does not do

The recipes are not yet read by the generator, so no generated value changes and no epoch is taken. Whether a world's life stage is chemical is decided by the life tables of other builders; this registry only fixes that the gases wait for it.
