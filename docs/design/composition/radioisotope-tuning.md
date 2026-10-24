# Radioisotopes: the retune of rows o2, o4 and o6 (GEN.514)

Status: built 2026-10-23 by lane M (builder B6). Rows o2 (thorium scale), o4 (Th/U) and o6 (aluminium-26 present) of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`); two targets in `physics/data/tuning_targets.json`; tests in `src/tests/test_tuned_draws.py`.

## 1. o2: thorium of non-thin-disc hosts

The row drew the thorium-to-magnesium ratio, over its solar value, log-uniform from 59 to 251 percent: the range of 14 solar twins. The Hypatia catalogue has 32 thick-disc stars with both elements: mean -0.063 dex, standard deviation 0.193. Thin disc (215 stars) gives -0.044 and 0.153, so there is no thick-thin offset in the ratio planets care about. The row now draws `log10 f ~ N(-0.063, 0.19)` (range -0.64 to 0.52). The width is the widest of the intrinsic spread of 53 solar twins (0.056, Botelho 2019) and the catalogue's 0.15 to 0.19, as the widest-range rule asks; it is an upper bound on the intrinsic width, because the catalogue mixes papers. Both targets pass (the mean -0.065, z -0.05; the spread 0.187, z -0.25).

Halo, young (under 1 Gyr) and M-dwarf hosts have no thorium measurement; they take the same draw as an extrapolation the validity range flags.

**Finding for the thin-disc row o1.** Row o1 centres the thin disc at +0.08 dex (the solar-twin scale of Botelho 2019, where the Sun is thorium-poor), the catalogue puts the thin disc at -0.044. The two scales differ by 0.12 dex, less than the catalogue's own spread, and the rows are not on one scale. Choosing the scale is open (filed to the TODO thread).

## 2. o4: thorium to uranium by formula

The row drew Th/U uniformly from 3.8 to 4.5, the Solar System value at 4.5 Gyr, whatever the planet's age. Th/U grows with age because uranium-238 decays faster (x1.17 over 1.5 Gyr). The row now uses the published formula (Frebel 2007 Table 1 for the production ratio 1.667 to 2.0, the half-lives 14.05 and 4.468 Gyr): `Th/U (mass) = 0.97483 PR F(T) exp[(l238 - l232) age]`, with `F(T)` the integral of constant r-process input over `T` before the host forms. The effective duration `T` is calibrated on one Solar System datum (CI chondrites need 6.5 Gyr), so the row draws between `T = 3 Gyr` and `T = 13 Gyr` (the widest bracket) and between the two production ratios: 3.04 to 5.39 at 4.5 Gyr. The row stays partly verified with that limit; the birth-epoch `T` needs a chemical-evolution model.

## 3. o3 and o6

- **o3 (K/U)** stays a labelled stand-in: planetary potassium is volatile-depleted and only Solar System bodies give the depletion. The sodium proxy of Harrison 2021 needs Boss's acceptance; it is listed as an open choice.
- **o6 (aluminium-26 present)** drew a share log-uniform from 0.1 to 1. The evidence is a mixture and not one distribution: at least 4 percent of systems enriched (Curry 2022), 9 to 28 percent of discs in cluster models (Lichtenberg 2016), tens of percent for white-dwarf debris, 66 percent differentiated (Bonsor 2020). The bracket is now 0.04 to 0.66 and no target exists (nothing states a count of a defined sample); the bracket is pinned by a test.

## 4. What it does not do

- The recipes are not yet read by the generator, so this changes no generated value and takes no generator epoch.
- The 5.2e-5 level of aluminium-26 itself is a Solar System value and is not tuned to.
