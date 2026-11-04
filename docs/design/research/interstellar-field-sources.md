# Interstellar asteroid fields: what the sources say (GEN.1098)

Status: read report. Written 5:36 PM EDT 2026-10-30 as the lane reported. Nothing built. Reads the sources of `docs/design/systems/interstellar-asteroid-fields.md` section 10 for the numbers the four field modules need (GEN.1099 collision, GEN.1100 dispersal, GEN.1101 rates, GEN.1102 composition). Every number is a row of the registry section 6.103 (rows IF1 to IF12 of `gap-filling-plan.md`). Boss decided physics only (4:54 PM EDT 2026-10-30): no game rate and no floor dial.

## What was read

| Source | How | What it gave |
|---|---|---|
| [R2633] | full text | the disruption threshold and its scaling (IF3), the largest remnant (IF1, IF2), the impact angle (IF4), the fragment tail (IF5) and the fragment speeds (IF6) |
| [R4226] [R4999] | abstract, sections 2 to 4 | where erosion begins, the debris share of giant impacts, and the core and mantle bookkeeping (IF7, IF8) |
| Su et al. 2006, ApJ 653, 675 [R4247] | abstract, sections 2 to 4, Table 4 | the share of young A stars with a belt, by age (IF9) |
| Do, [R1124] | abstract, section 2 | the interstellar background the fields stand against (IF11) |

Not reached, and what that leaves open:

* **[R3191] and 2020, [R2178]**: rogue density; already the root of `tuning.ROGUE_PLANET_RATE_PER_STAR`, and GEN.1109 (rogue abundance) owns them. The collision rate (design section 2) stays a derived figure on that input.
* **Shock vaporisation thresholds** (Pierazzo, Melosh) and **[R0527]** (flyby truncation): not opened. The design's remark that peak shock above about 400 GPa leaves mostly melt and vapour stays [R]. LS12 has no vaporisation at all (see finding 3).
* **[R0763]** (excess rates of B stars): not opened; IF9 uses A stars as the stand-in for the belts of supernova progenitors.
* **[R2867], 2010** (core and mantle in impacts of differentiated bodies): quoted second-hand through [R4226] (IF8).
* **[R1941]** (unbinding by mass loss): the half-mass rule is derived in IF10 in two lines instead.

## What it changes in the design

1. **The rubble window is not 15 to 30 km/s; it scales with the pair.** The disruption speed is `V* = (32 pi c* / 5)^(1/2) (rho1 G)^(1/2) RC1` with `c* = 1.9`, `rho1 = 1000 kg/m3` and `RC1` the radius of the combined mass at 1000 kg/m3 (IF3). Computed here [C]: two Earths 22.6 km/s (`RC1` 14 180 km, `Q*` 6.4e7 J/kg, escape speed of the pair 10.6 km/s), two 0.3-Earth bodies 15.2 km/s, two 0.1-Earth bodies 10.5 km/s. The window of a large remnant plus rubble is `QR / Q*` from 1 to 1.8, i.e. speeds `V*` to `1.34 V*` (22.6 to 30.3 km/s for two Earths). Below `V*` the pair merges or hits and runs; above it the largest remnant falls under a tenth of the mass (0.044 at 40 km/s, 0.013 at 60 km/s for two Earths). The design's window was right for two Earths and wrong for lighter pairs, where it is lower.
2. **The largest-remnant law replaces the shock table.** The design's section 3.1 outcome table (merger, partial, catastrophic by speed alone) becomes `Mlr / Mtot` against `QR / Q'*RD` (IF1, IF2), with the impact angle drawn from `sin theta cos theta` (IF4) and the hit-and-run test `b > R / (R + r)`.
3. **Extrapolation warning.** The law was fitted to N-body simulations of gravity-dominated rubble at speeds up to a few km/s (impact speeds of the table: tens of m/s to km/s). Applying it to 20 to 80 km/s uses the same scaling far outside the fitted range, and it has no melt or vapour: the module must say so in its history sentence and keep the vapour-and-droplet outcome of design section 3.1 for speeds above about `1.34 V*`, where LS12's own 10 percent largest remnant ends the law's region of validity.
4. **Unbound mass is smaller than the design's "near 1".** Read as dm/dv = 10^(A - S v) per unit of the escape speed with `A = -0.3 Mlr / Mtot + 0.3` (IF6), the share of the mass that leaves faster than the escape speed of the combined mass is about 11 percent at the start of the super-catastrophic regime and tends to 13.6 percent of `Mtot` at high speed [C]; the rest stays in a cloud bound to the (tiny) remnant and to itself. That is a reading of equation 43 (the bin width is not given in the text; the reading matches the paper's "lowest bin about 0.1 `Mtot`"), so GEN.1099 checks it against Figure 7 first. If it holds, a collision between planet-size rogues makes a bound debris swarm round the remnant plus a smaller unbound field, not a field of the whole mass.
5. **Composition from two parents.** The escaping fragments come mostly from the silicate mantles of both bodies, the cores going into the largest remnant (IF8, average of the two idealised models, initial core fraction 1/3). So an unbound collision field is silicate-rich with little iron, not iron plus silicate plus ice in proportion; ice shards only from water-rich parents' outer layers. GEN.1102 starts from that.
6. **Supernova belts are better supported than the design's 0.3.** 41 percent of A stars of 30 Myr or less carry a belt at 24 micron (48 percent at 70 micron), falling to 25, 15 and 2 percent in older bins (IF9, Table 4 of Su). Progenitors are 8 solar masses and more and live under 40 Myr, so 0.41 replaces 0.3 as the stand-in (tier T3: A stars stand in for B stars, and the excess is dust, not belt mass). The belt is unbound whenever the remnant is under half the progenitor's mass (IF10): every neutron-star case, and a black hole of less than half.
7. **A field past 1 ly sits under the interstellar background.** `n_IS` about 0.2 per cubic AU for 100 m bodies (IF11) is the floor the retired field falls into; it is one detection, so factor 3.

## What each module takes

| Item | Rows | Notes |
|---|---|---|
| GEN.1099 `rogue_collision.py` | IF1 to IF7, IF12 | `V*`, `QR`, the remnant, the tail, the speed law; the collision rate of design section 2 stays the derived figure (not a row: [C] from the rogue density) |
| GEN.1100 `field_dispersal.py` | IF6, IF11, IF12 | expansion from the speed law (not a flat 5 to 15 km/s), the retirement radius at 1 ly |
| GEN.1101 `field_rates.py` | IF9, IF10 | supernova belt rate from the core-collapse rate times the belt share; no floor dial |
| GEN.1102 composition | IF8 | mantle-dominated shards, core in the remnant |

## Open for the module lanes

* The reading of equation 43 (finding 4) and the unbound share it gives.
* Whether the remnant planet (own rogue planet, Boss decision) is the "largest remnant" of IF1, with its core fraction from IF8.
* Which `c*` to use for a rogue with a thick ice shell (IF3 gives 1.9 for planets, 5 for small ice or rubble bodies; the water-rich rogues sit between).
