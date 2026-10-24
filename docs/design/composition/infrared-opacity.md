# The grey infrared optical depth of an air from its mixture (GEN.703)

Module: `src/planetgen/physics/infrared_opacity.py`. Tests: `src/tests/test_infrared_opacity.py`, `src/tests/test_atmosphere_state.py`,
`src/tests/test_greenhouse.py`. Lane E, 2026-10-24. Research: part 03 of `atmosphere-composition-parts/` (sections 2, 3, 5.4, 6.2).

## What it does

`atmosphere_state` (GEN.423) gave every gas one grey mass opacity, so forty parts per million of carbon dioxide on a temperate world
contributed a thousandth of an optical depth and nothing depended on temperature. Here `tau0` is the equivalent grey thickness of
Koll (2019, Eq. 13), route 1 of part 03 section 6.2:

    tau0 = -ln( sum_k w_k(T) exp(-tau_k) )

on a grid of 10 cm-1 bins from 0 to 6000 cm-1 (and a last bin to infinity). `w_k(T)` is the share of a blackbody at the surface
temperature in bin `k` (exact Planck integrals), `tau_k` the column optical depth in it. A window holds `tau0` down however thick the
rest is; a band under the Planck peak of the surface matters far more than one far from it; overlapping gases share a bin.

| Part | Form | Source |
|---|---|---|
| Bands | each band of the gas table (centre or range) is a plateau `kappa0 u` with wings `exp(-d / 10 cm-1)`; `u` the gas's column, kg/m2; `kappa0` scales with the pressure | band positions: gas table; the 10 cm-1 e-fold from the logarithm of the CO2 forcing law (`F = 5.32 ln(C/C0)`, part 03 section 2.2) |
| Far wings | `kappa_w u (P/P_ref)` within `600 (P/P_ref)` cm-1 of a band | Koll 2019 text: pure H2O thick above ~1e-2 bar, pure CO2 of order unity at 1 bar, CO thin |
| Water continuum | `(sigma_self x_w + sigma_for (1-x_w)) x_w n0^2 H / (2 n_ref)` per bin, 200 to 3600 cm-1 | MT_CKD 4.3 table of part 03 section 3.2 (computed from the AER file) |
| Collision-induced pairs | `ln(1 + sum kP(T) x1 x2 n0^2 H / 2)` over each pair's range | HITRAN table of part 03 section 3.1 |

`atmosphere_state.greenhouse_warming` takes the opacity at the surface temperature, which the closure gives, so the two are solved
together (four to six steps of the fixed point; 1 to 3 ms for the opacity, 3 ms for the closure each step). The closure's pressure
exponent `n` is 2 for every air; the lapse scaling `alpha` runs from 0.6 (condensable vapour supplies all the opacity) to 0.8 (none
does), and 0.85 for an H2/He air (part 03 section 5.3: Earth 0.6, Venus 0.8, Titan 0.77, Jupiter 0.85).

The closure itself (`greenhouse.py`) had a bug this found: the continuity integral stepped 400 times over the whole of a thick
atmosphere, so above a `tau0` of about 500 it stepped over the part that counts and returned temperatures several times too high.
It now integrates over the 60 e-folds that carry the integral, and reproduces the table of part 03 section 5.4 to its two digits
from `tau0` 1 to 2e5 (`test_a_very_thick_atmosphere_follows_the_table_of_part_03`); it is also a few times faster (numpy).

## Untuned results (checks, not inputs)

| Air | `tau0` | `T_s` K | `G` | Measured or published |
|---|---|---|---|---|
| Earth, 1 kPa well-mixed vapour | 1.81 | 289 | 0.65 | `tau0` 1.82, `G` 0.625 |
| Earth, 400 Pa (the measured 25 kg/m2 of vapour, well-mixed) | 1.28 | 282 | 0.49 | `G` 0.625 |
| Earth, no vapour | 0.27 | 263 | 0.13 | |
| Venus, 92 bar CO2 | 2.8e4 | 682 | 76 | `G` 99.6 (735 K) |
| Mars, 600 Pa | 0.17 | 214 | 0.075 | `G` 0.10 to 0.21 |
| Titan, 1.47 bar | 3.4 | 102 | 1.4 | `G` 0.69 (a haze world; the paper's `tau0` is 5.6) |
| H2/He, 1 bar | 2.0 | 144 | 1.0 | `tau0` 6.3 for Jupiter at 1 bar |

The collision-induced pairs reproduce the cases of part 03 section 3.1 near the Planck peak to 25 to 50 percent (3.3 bar N2 + 10 percent
H2: 1.26; 1.1 bar: 0.25; the Ramirez 2014 CO2 + H2 case: 2.4; the hycean H2 case: 10) and are within a factor of 3 for the cold ones.
Doubling carbon dioxide on a state with its water held fixed: the band model changes the outgoing flux by 3.7 W/m2 (the measured
forcing) for columns from 140 to 1120 ppm, but the grey depth that results warms the surface by 0.3 to 0.6 K, a third to a half of the
1.1 K of a full radiative-convective model (the equivalent depth counts a saturated band as black where a real one radiates from the
cold air above; a well-mixed water column also masks the 15 um band further).

## Limits (all stated in the module)

- The repo holds band centres and no line strengths. `kappa0`, the far-wing coefficients, the hot-band (1e-6) and Fermi dyad (1e-5) factors,
  the 600 cm-1 wing reach and the fingerprint plateau of gases without a listed band are STAND-INS, labelled in the gap plan. `kappa0`
  enters only through `ln(kappa0 u)`, so a factor of a few in it moves nothing much. The 15 um band of CO2 alone gives the 3.7 W/m2 a
  doubling is measured to; a hot band or dyad counted fully would add two to three times that.
- No condensation: the stored vapour is taken as it is, well-mixed. Earth's measured vapour column is that of about 400 Pa well-mixed
  and the surface vapour pressure about 1 kPa; both land within 25 percent of Earth's `G`. A vapour profile that falls with height
  (the Clausius-Clapeyron scale height) would unmask the 15 um band and raise the sensitivity to CO2; it is a follow-up.
- The cold haze world (Titan) is two to three times too warm: the pair `ln(1 + x)` overstates N2-N2 at 94 K (3.3 against 1.24), and the
  single-channel closure has no haze (the paper itself needs two channels and `n` varying with pressure).
- One `T` for the Planck weights (the surface); no temperature dependence of the band strengths.
- Real k-tables (correlated-k, offline with ExoMol or HITRAN data) would replace the band and wing numbers one for one; nothing else
  here changes. That is a data item for Boss, since it adds a download at setup.

## Wiring

`atmosphere_state.greenhouse_warming` is the one entry point. It is not on the generation path yet: lane I's GEN.411 rewrite of
`physics/planets.py` draws the mixture first and calls it, which takes the generator epoch. No generated value moves with this change.
