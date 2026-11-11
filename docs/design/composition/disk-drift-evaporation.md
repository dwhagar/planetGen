# Drift and evaporation of water and the other ices in the disc (GEN.1472)

Item GEN.1472, Lane 14. A design note, built by GEN.1473 (`physics/disk_drift.py`; section 7 says what the build changed and what it measured). Boss's rule for this note (7:44 PM EST 11-08): build the model, secondary and tertiary citations are accepted by default, and speculative theory that fits the physics is allowed, as long as every part says which it is. Every stand-in is a row in the registry (section 6.242 of `docs/design/research/gap-filling-plan.md`, rows DE1 to DE10).

## 0. What this note decides

1. The disc's gas does not keep the composition of the static budget (`disk_condensation.gas_and_solid_c_over_o`). Pebbles grow, drift inward faster than the gas, and evaporate where they cross each ice line. The gas inside a line gets that line's ice at a rate set by the pebble flux over the gas flux, not at the ratio the budget holds.
2. One number does the work: the **enrichment ratio** `E_L(t)`, the flux of a line's ice carried in by pebbles over the flux the static budget would carry with the gas. `E = 1` is the static model of GEN.1363 and must reproduce it exactly (a test the build carries). `E > 1` while pebbles are plentiful, `E = 0` before the front arrives and after the pebbles are gone.
3. Everything else in the model is a clock (when the pebbles arrive, when they stop, how long the gas takes to carry the vapor to a place) or a leak (how much vapor goes back out and freezes again, the cold finger).
4. The model has no new random draw. Its inputs are already drawn and stored: the star's mass and luminosity, its drawn snow line, the disc's outer edge (18 times the snow line, `tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER`), its metal ratio, and a giant's core mass.

## 1. What is missing

`disk_condensation` holds a disc that stays as it was born. Water is a vapor inside 182 K and ice outside, the CO2 and CO lines do the same further out, and the C/O of the gas at a place is a function of the temperature at that place only. That is right for a disc in which solids and gas move together. They do not. [R5208] (abstract read) says pebbles carry volatiles inward and enrich the gas at snow lines; [R5209] (abstract read) says particles drifting through an evaporation front raise the local vapor by more than an order of magnitude under nominal conditions, that the enhancement spreads inward, and that solids are enhanced outside the front, for any condensible; [R5210] (abstract read) says an atmosphere dominated by vapor from evaporating pebbles lets a giant hold a large heavy-element content, that planets formed farther out tend to have higher C/O, and that low C/O goes with high heavy-element content. [R3322] (read) gave the static picture we already use. [R5212] (title only) is the same question for the effect of drift and viscous accretion on the lines.

Two things in the generator want this: the C/O of the gas a giant gathers (GEN.1363 part C, GEN.1471), and the ice-to-rock ratio of the solids just outside the water line (the jump `SNOW_LINE_SOLID_JUMP`, GEN.1396).

## 2. The model

### 2.1 The clocks

**The pebble front (DE1, DE2).** Dust grows where the growth time is shorter than the age, so a front sweeps outward, `r_g(t) = 122 AU x (t / 1 Myr)^(2/3) x (M* / Msun)^(1/3)` (the form and exponents of [R5211], recalled; the number follows from a growth efficiency of 0.05 and a dust-to-gas ratio of 0.01, both recalled). The pebbles the front leaves behind drift inward. The line of an ice at distance `r_L` (`disk_condensation.line_au`, for the Sun 2.3 AU for water, 35 AU for CO2, 194 AU for CO) sees its pebbles from `t_L = 1 Myr x (r_L / (122 AU x M*^(1/3)))^(3/2)`.

**The pebble flux (DE1).** `Mdot_F(t) = 9.5e-5 Earth masses per year x (Z / 0.01)^(5/3) x (M* / Msun)^(1/3) x (t / 1 Myr)^(-1/3)` ([R5211], recalled; `Z` is the solids over gas of the star, `disk_condensation.masses_percent(15 K)` summed, 1.43 percent for the Sun, scaled by the star's metal ratio). The flux ends when the front reaches the disc's edge, `t_end = t_L(edge)`: the dust beyond is gone. The edge is the drawn 18 snow lines, so `t_end` is 0.25 Myr for the Sun (edge 49 AU), 0.017 Myr for a 0.3 solar-mass star and 0.63 Myr for 1.5 solar masses. This is speculative (DE3): a real disc loses its pebbles gradually, over a few Myr, not at a step. The step is the simplest form that conserves the dust. The build keeps the step and puts a decay time of one growth time at the edge on the stand-in list if the measure asks for it.

**The gas flux (DE4).** The disc accretes viscously, `Mdot_g(t) = Mdot_0 x (1 + t / t_nu)^(-3/2)`, with `Mdot_0 = Md_0 / (2 t_nu)`, `Md_0 = 0.05 Msun` and `t_nu = 0.5 Myr` for the Sun (typical Class II accretion rates are 1e-8 Msun per year at 1 to 2 Myr, which this gives). The disc mass scales as `M*^1.8` (`tuning.DISK_MASS_STELLAR_MASS_EXPONENT`, [R0114]), and **cancels**: the flux of pebbles of [R5211] is proportional to the same disc mass, so `E` does not depend on it.

**The lag (DE5).** The gas that crosses a line takes time to reach a place inside it: `tau(r_L to r) = (r_L - r) / (3 nu_1)` with `nu = nu_1 x (r / 1 AU)` (a steady disc with `alpha = 1e-3`, 0.05 Myr per AU). That is 1.2 Myr from the CO2 line (35 AU) to 10 AU, and 0.1 Myr from 5 AU to 3 AU. The composition at `r` at time `t` is the mean of `E_L` over the window `[t - tau, t]`. Inside a few AU the lag is under 0.1 Myr and the quasi-steady form is right; beyond 10 AU it is not, and the CO2 and CO signal arrives a Myr after the pebbles that carried it.

### 2.2 The enrichment ratio

For the ice of species `L` (water, CO2, CH4, CO; ammonia, H2S and the noble-gas clathrates follow the same form and ride with their parent line in the first build):

`E_L(t) = gate_L(t) x [Mdot_F(t) / (Mdot_g(t) x Z)]`

`gate_L = 1` for `t_L <= t <= t_end` (and the line inside the disc's edge), else 0. A pebble that left the disc beyond the line carries the ice of every species whose line is outside where it started, so the gate is on the front, not on the pebble. The ratio is the pebbles' mass flux over the gas's, in units of the solids-to-gas ratio of the disc: 1 when the solids and the gas move together. For the Sun it comes out near 2 to 4 for the whole pebble epoch (2.3 at 0.05 Myr, 2.1 at 0.25 Myr).

### 2.3 The gas at a place and a time

The static gas at temperature `T` (`_gas_and_solid_atoms`) is the **supply** (the gas at the disc's edge, where the species with lines warmer than the edge are all frozen) plus the ice atoms of every species whose line lies between the edge and `T`. The drift model keeps the supply and weights each ice by its `E_L` (lagged as in 2.1):

`gas(T) = supply(T_edge) + sum over L in (T_edge, T] of  E_L x ice_atoms_L`

That has the property the build tests: `E_L = 1` gives `_gas_and_solid_atoms(T)` exactly. The C/O of the gas is the carbon over the oxygen of that sum; the rock and the refractory carbon are untouched (500 K, DSC3).

### 2.4 The cold finger (DE6)

Vapor that crosses a line outward freezes again on the grains just outside it, and those grains drift back in: a loop. If a fraction `leak` of the vapor goes outward, the steady vapor abundance at the line is `1 / (1 - leak)` times what the incoming pebbles alone would give, the ice just outside is enriched by the same factor, and the **net** flux inward is unchanged (it equals the pebble flux). So the cold finger does not change the gas composition inward of a line; it changes two other things: the ice-to-rock ratio of the solids in a band just outside the line, and the vapor at the line. [R5209] gives an order of magnitude (`leak` 0.9); the build takes `leak = 0.9` (DE6) and a band as wide as the gas scale height at the line. The band's rise of the water ice over rock is `1 / (1 - leak) = 10` at most, to be cut by the share of the pebble flux that is ice. It feeds the jump `SNOW_LINE_SOLID_JUMP` only if the measure of GEN.1473 shows the C/O of the solids needs it; the first build leaves the jump at its GEN.1396 value and reports the band.

### 2.5 What a giant gathers

A giant takes the gas of its place and the solids it swallows (`envelope_c_over_o`, GEN.1363 part C). The change here is only where the gas comes from: the gas of section 2.3 at the giant's formation temperature and **at the age it gathers its gas** (DE9). That age is not drawn today. The build derives it without a draw: the time the cumulative pebble mass `integral of Mdot_F dt` from the front's arrival at the giant's formation distance reaches the core mass in the recipe (the critical core, 10 to 20 Earth masses), capped at the disc's lifetime (3 Myr, `scatter.DISK_LIFETIME_YR`). A core that the pebble epoch cannot build in time gathers its gas at `t_end` with a post-pebble core (planetesimal accretion) and sees the pebble-depleted gas.

## 3. First numbers (an offline script, not repo code; the static values it departs from are 0.895 outside the water line and 0.381 inside)

Solar-mass star, edge 49 AU, no lag (so the CO2 numbers are an upper bound on how early it arrives):

| age | `E` (water) | gas C/O inside water line | gas C/O 4 snow lines out |
|---|---|---|---|
| 0.05 Myr | 2.3 | 0.19 | 1.00 (no CO2 yet) |
| 0.10 Myr | 2.1 | 0.21 | 1.00 |
| 0.25 Myr | 2.1 (CO2 front arrives) | 0.26 | 0.83 |
| 0.50 Myr and later | 0 | 1.00 | 1.00 |

Static model: 0.38 inside, 0.895 outside. The 1.5 solar-mass star behaves alike (0.17 to 0.21 inside during its 0.63 Myr pebble epoch, 0.73 outside in the CO2 window); the 0.3 solar-mass star's pebble epoch is over before 0.02 Myr, so its gas is pebble-free from the first 0.02 Myr on and sits at 1.00 everywhere it is observable.

What this says, before the build:

- **Inside the water line the drift model is firmly lower than the static one during the pebble epoch (0.2 against 0.38)** and firmly higher after it (1.0). That is the direction of the brief and of [R5210] (inner disc: low C/O early, high late).
- **Outside the water line the effect is modest and goes the other way for most of the epoch.** The water is frozen either way. What differs is the CO2 and CH4 the gas holds: pebbles with those ices arrive at the CO2 line (35 AU) only when the front gets there, near the end of the epoch for the Sun, and the gas needs a Myr to carry them in. For most of a disc's life the gas between the lines is the **supply gas, C/O about 1.0**, not the static 0.895; only the lagged CO2 window brings it down (0.83 at `E = 2`, 0.73 for the 1.5 solar-mass star). The brief expected gas outside the water line to run the other way from the gas inside (high early, low late). The model gives high early (1.0) and, after the lag, lower while CO2 arrives (0.8), then high again once the pebbles stop. The honest reading: the gas side outside the water line is a weak lever and may push the over-0.8 share of giants **up**, not down (GEN.1471 left it at 22 to 26 percent against the target of 4 to 9). The build measures it and reports the sign. If the measure goes the wrong way the levers left are the share `m` of the envelope that is metals (`ENVELOPE_METAL_SHARE`, one re-check allowed by GEN.1473) and what a giant swallows (section 2.5, the solids); the gas model is not tuned to hit the target.
- **The cold finger rises the water ice of the solids just outside the snow line**, which lowers the C/O of the solids a giant swallows there. That is the lever that points down for giants that form near the water line.

## 4. Which parts are secondary or speculative

| Part | Basis | Standing |
|---|---|---|
| Ice moves in with pebbles, evaporates at its line, enriches the gas inside; solids enhanced outside | [R5208, R5209, R5210] abstracts read | secondary, read |
| `E = 1` is the static budget | derived from `_gas_and_solid_atoms` | exact |
| Front `r_g(t)`, flux `Mdot_F(t)` formulas and numbers | [R5211] form recalled, numbers recalled, paper not opened | recalled, T3 |
| Viscous `Mdot_g(t)`, `Md_0`, `t_nu`, `alpha` | standard disc theory, observed accretion rates; the numbers are mine | educated guess |
| Pebbles end at the disk's edge as a step (replaced by a taper, GEN.1514, section 10) | my simplest form | speculative |
| The lag as a boxcar of the viscous time | my simplest form | speculative |
| Cold finger `leak = 0.9`, band width | [R5209] abstract's order of magnitude; the width is mine | secondary and speculative |
| The giant's gas-gathering age from the cumulative pebble mass | my derivation to avoid a draw | speculative |
| Pebble flux cut by a giant already formed (isolation mass), planets accreting pebbles on the way | not modelled | open point |

## 5. Wiring for GEN.1473

- New module `physics/disk_drift.py` (the clocks, `enrichment_ratio(species, t, star)`, `gas_atoms(T, t, star)`); `disk_condensation` keeps the static budget and gets `gas_and_solid_c_over_o(..., enrichment=None)` so `None` is today's behavior. `envelope_c_over_o` takes the gas atoms from `disk_drift` when it is given an age. `mixture_giant._disk_c_over_o` passes the giant's gas-gathering age. No shim: the static call is the `E = 1` case, not a second code path.
- Reuse: `LINES_K`, `line_au`, `_gas_and_solid_atoms`, `inventory`, `masses_percent` (the budget), `snow_line_solid_jump`, `formation.draw_snow_line_au` (the stored snow line), `DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER` (the edge), `DISK_MASS_STELLAR_MASS_EXPONENT` (cancels, documented), `scatter.DISK_LIFETIME_YR`, the star's luminosity track, the giant's formation distance and core mass. Nothing here needs a look-alike of an existing function.
- Tests: `E = 1` reproduces the static gas and C/O exactly at every temperature; mass conservation (the C and O in the gas plus the ice left behind equals the budget); `E = 0` gives the supply gas; the lag window never exceeds the age.
- Measure: the same 1,000 G2V systems and 331 fresh giants as the GEN.1395 report, reporting the share over 0.8 and under 0.35 and the median, before and after, and the stars that are not G2V (M dwarf, 1.5 solar masses) for the sign of the epoch lengths.

## 6. Open points

- The numbers of `Mdot_F` and `r_g` are recalled. The build reads [R5211] (the PDF through the web-fetch route) and replaces the recalled numbers with the read ones; DE1 and DE2 move from T3 to T1.
- The gas side outside the water line may not lower the share of giants over 0.8; section 3. The result decides what GEN.1473 reports; it does not decide the model.
- The step end of the pebble flux makes M dwarfs pebble-free almost from birth (0.017 Myr), while observed M-dwarf discs hold pebbles for Myr. A smoother decay (one growth time at the edge) is the first thing to try if the M dwarf rows look wrong.
- `leak = 0.9` is read off an abstract's order of magnitude; the net flux does not depend on it, the band does.
- Not modelled: a giant that has formed blocks the pebbles for the inner disc (the isolation mass); planets that eat pebbles on the way; ice in the soot line; photoevaporation of the gas. Each would lower `E` further in; none is needed to set the sign of the first build.
- The moons' circumplanetary discs (`moon_disk_zone`) use their own lines; they are out of this note.

## 7. What the build changed and measured (GEN.1473)

**Built.** `physics/disk_drift.py` holds the clocks, `enrichment_ratio`, `ice_weights`, `gathering_age_myr` and `giant_ice_weights`; `disk_condensation._gas_and_solid_atoms` and `envelope_c_over_o` take `ice_weights` (`None` and weights of 1 are the static budget exactly, a test); `mixture_subject` gives a giant its stored snow line (`disk_snow_line_au`); `mixture_giant._disk_c_over_o` passes the weights. No new random draw. Stand-ins DE1 to DE10 are rows of section 6.242, and DE1 and DE2 are now read (T1): the pebble flux and the front are Eq. 14 and Eq. 10 of [R5211], and the growth efficiency of 0.05 that gives the 122 AU is the paper's sticking efficiency of 0.5 over its grain growth ratio of 10. The paper does not state when the flux ends (DE3 stays my step at the disc's edge).

**Differences from the note.** (1) The star-mass term of the pebble flux (the paper's 1/3 power) is replaced by the disc-mass scale (the 1.8 power of `tuning.DISK_MASS_STELLAR_MASS_EXPONENT`), so the disc mass cancels in the ratio as section 2.1 says; with the 1/3 power it would not cancel and the ratio of a 1.5 solar-mass star would be half the Sun's. (2) The nitrogen, sulfur and noble gases of `envelope_enrichments` keep their static values; the first build moves the carbon and oxygen only. (3) The doc's 49 AU edge is 18 snow lines of a 2.7 AU snow line; a 2.3 AU snow line (the Sun's) gives an edge of 41 AU and a pebble epoch of 0.20 Myr.

**Measured** (reports/gen1473-drift-measure-2026-11-09.md, 2,787 generated giants of all stars): for G and K stars the share of giants over C/O 0.8 falls from 18 to 12 and 11 percent and the median is unchanged (0.42); the share under 0.35 rises by 2 points; F stars do not move. M dwarfs move a lot (median 0.62 to 0.44, over 0.8: 39 to 12 percent) because 99 percent of their giants gather their gas as the pebbles end, within the lag from the line; that row stands on DE3 and DE5, which are the two stand-ins to revisit first.


## 8. The start state of the gas, and the noble gases (GEN.1517, GEN.1516)

**The start state (GEN.1517).** Section 2 started the gas at the disc's edge, where every ice whose line is warmer than the edge is frozen, and gave a place only the vapor that pebbles had carried in, lagged by the gas flow. That is the state of a disc that has been drained of its ices, not of a new one: at time zero the gas inside each line is the solar gas with that line's vapor in it (the static budget of GEN.1363), and it stays so until gas from outside the line has flowed in. The viscous flow carries gas inward at about 20 AU per Myr at the lag of DE5, so the gas at a place at age `t` began `t / LAG_MYR_PER_AU` AU farther out. `disk_drift.enrichment_ratio` now returns `E = 1` for a species whose line the gas has not crossed (`age < lag_myr(line, place)`); past it, the gas crossed at `age - lag` and the ratio of section 2 applies. Stand-in DE11.

Measured on a grid of stars (0.2 to 1.4 solar masses, [Fe/H] -0.3 to 0.3, formation temperatures 60 to 500 K): before the fix a Sun-like giant's gas held no argon, krypton, xenon or CO2 at all (every weight 0); with it, 18 of 24 cells of a 1 solar mass star and 20 of 24 of a 1.4 solar mass star are the static disc exactly, and the dwarfs keep more drift (12 of 24 at 0.2, 9 of 24 at 0.4). The cause of the rest being static is the age: the gas-gathering age of DE9 is 0.04 Myr at the median, shorter than the lag from any line at least a few AU away, so the drift reaches only giants close to a line or around dwarfs. Whether 0.04 Myr is a plausible age for a giant's gas (a core of 15 Earth masses built from every pebble that passes, with no accretion efficiency) is the next question for the drift model.

**The noble gases, nitrogen and sulfur (GEN.1516).** `disk_condensation.envelope_enrichments` takes the same `ice_weights` as `envelope_c_over_o`: the vapor of N2, NH3, H2S, Ar, Kr and Xe that evaporated on the way is `E` times the static vapor, the solids stay the static ones, and `sigma` is found from the weighted carbon and oxygen as in section 7. With weights of 1 (or none) the factors are the static ones exactly. A giant that formed colder than a line has that ice frozen and takes it with the planetesimals as before; the weight acts only on a gas warmer than the line.

## 9. The age a giant gathers its gas, from the efficiency of pebble accretion (GEN.1518)

Section 7 built the core from every pebble that passes its place (100 percent efficiency), which put the median gas-gathering age of a giant at 0.04 Myr, shorter than the gas needs to carry the vapor of nearly any line to the place. A growing core accretes only a share of the pebbles that drift past it: [R5211: Eq. 33] gives `f = 0.034 (Mc / M_E)^(2/3) (r / 10 AU)^(-1/2)` for pebbles of Stokes number 0.1 (the Hill-sphere accretion of Eq. 28), and a seed of 1e-3 Earth masses (section 3.3). Because `f` is a power of `Mc`, `Mc^(1/3)` grows in proportion to the pebble mass that has passed: `Mc^(1/3) = Mseed^(1/3) + (f_0 / 3) (r / 10 AU)^(-1/2) M_pebbles`. `disk_drift.pebble_mass_for_core` gives the pebble mass a core needs (above `f = 1`, close in and for small cores, every pebble is kept); `gathering_age_myr` takes the age at which that mass has passed since the front arrived. Stand-ins DE9 (now through DE12) and DE12.

The 15 Earth mass critical core needs 66 Earth masses of pebbles at 1 AU, 148 at 5 AU and 209 at 10 AU, about 1 Myr of the pebble flux of the Sun at 5 to 10 AU, which is the paper's "within about 1 Myr" ([R5211] sections 5 and 7, read only in summary here). In the generator that time is longer than the pebbles last for nearly every star: the pebble epoch ends when the front leaves the disk's edge (DE3), 0.2 Myr for the Sun and 0.004 Myr for a 0.2 solar mass dwarf (an edge of 18 snow lines), so on a grid of stars most giants now gather their gas at the end of the pebbles (24 of 24 cells for a 0.2 and a 0.4 solar mass star, 16 of 24 for the Sun). What the drift model then does for a Sun-like giant is set by when the pebbles end, not by the core: the next stand-in to revisit is DE3 (GEN.1514, a smooth decay of the pebble flux instead of a step at the edge).

## 10. The pebble flux tapers instead of stopping (GEN.1514)

Section 4 ended the pebble flux as a step to zero when the front reached the disk's edge (DE3). A front that reaches the exponential cut-off of the gas disk meets less dust, so the flux falls as the gas surface density does at the front's place: `exp(-front / edge)` (`disk_drift.pebble_taper`). The taper is 1 for a young front, 0.37 when the front is at the edge (0.2 Myr for the Sun) and 0.002 at the end of the disk's life, and it multiplies the pebble flux in the enrichment ratio (`tapered_ratio`) and in the pebble mass that builds a core (the integral has a closed form in the age to the 2/3). The pebble epoch no longer has an end age (`pebble_end_myr` is gone); a species is carried until the disk is, `disk_end_myr` (3 Myr).

Measured on the grid of stars (24 cells each): a core of 15 Earth masses at 5 AU needs 148 Earth masses of pebbles (section 9) and the Sun's tapered flux inside its edge brings 48, so a giant of a dwarf or Sun-like star cannot build its core from pebbles inside the generator's disk and gathers its gas as the disk goes (24 of 24 cells for 0.2 and 0.4 solar mass stars, 16 of 24 for the Sun, 9 of 24 for 1.4 solar masses). The gas that reaches it by then crossed the lines after the pebbles had stripped the ices: its weights are 0.02 to 0.05 for CO2, H2S and the noble gases, and the water weight is under 1. Effects on a giant at 3 times its star's heavy elements: the median carbon over oxygen is unchanged (0.49 to 0.48, the planetesimals' share makes up for the gas), the nitrogen and sulfur stay near the Galileo probe's values (3.0 and 2.8 against 2.96 and 2.75 measured), and the argon, krypton and xenon fall from 1 times the star's to 0.03 to 0.05. That last row is a miss against the probe's 2 to 3 times solar (the model's earlier miss was 1): in this model a giant gets its noble gases only from solids that froze them (a giant that formed colder than 48 K), and the planetesimals of a giant at 100 to 200 K carry none. The follow-up that would repair it is a provenance for the planetesimals a giant accretes, not a change to this taper. The disk edge itself (18 snow lines; a 0.2 solar mass star's is 0.9 AU) is the largest lever in these numbers.

## 11. Where the planetesimals formed (GEN.1519)

With the pebble flux tapered (section 10) the gas of most giants has almost none of the noble gases, and the planetesimals a giant accretes were taken to be the ones that froze at its own place, which for a giant warmer than 48 K have none either. A giant does not accrete only what formed where it sits. `disk_condensation.provenance_condensed` takes the mean of the frozen species over the annuli between the giant's place and the disk's edge (equal per AU, 64 annuli), and `envelope_c_over_o` and `envelope_enrichments` take that mix for the solids (`solids_edge_k`, which `mixture_giant` sets to the disk's edge temperature) while the gas keeps the ices of its own place and the drift weights. Stand-ins PRV1 to PRV3.

Measured on the grid (24 cells per star, envelopes at 3 times the star's heavy elements, median over the cells): argon 0.95 to 1.9, krypton 1.65 to 2.5, xenon 2.9 to 3.6 times the star's from dwarfs to 1.4 solar masses (before: 1 with the static solids, 0.03 to 0.05 with the tapered gas alone); the Galileo probe measured 2.5, 2.7 and 2.6. Nitrogen is 2.3 (the probe's 2.96) and sulfur 3.1 to 3.3 (2.75 plus or minus 0.66); the carbon over oxygen is unchanged at the median (0.49). For a giant at 110 K around the Sun at 4 times the heavy elements, the pattern is Ar 1.4, Kr 2.4, Xe 4.3, N 3.0, S 4.3. The noble gases rise in the order of their lines (xenon most) because the share of the disk colder than a line is larger for a warmer line; the probe's equal enrichment would need solids from colder than argon's 48 K, and the generator's edge is at 43 K, so none exist. Whether the edge should be larger is GEN.1520. Sulfur is high because the H2S ice of the planetesimals beyond 80 K adds to a giant that formed warmer: a miss against the probe the item does not try to tune away.
