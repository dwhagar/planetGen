# Mantle oxides (GEN.1362)

A silicate rock stores seven more oxides than SiO2, MgO, CaO, Al2O3 and FeO: Na2O, K2O, TiO2, Cr2O3, MnO, NiO and P2O5. Without them the CIPW norm of a mantle has no albite, orthoclase, ilmenite, chromite or apatite, and the felsic recipe cannot come out of the derived chain. The model is in `planetgen/physics/mantle_oxides.py`, called by `bulk_draw.dry_partition`. No new random draw is made: the seven follow from the star, the formation zone and the rock's FeO that are already drawn.

## The chain, for one Mg atom of the body

1. **The star.** Atoms of the element per Mg atom of the host (`HostRatios`, from GEN.1360: Na, K, Ti, Cr, Mn, P, with Ni).
2. **The zone.** The share the zone kept against CI chondrites, `volatile_depletion.retained_fraction(T50, cut)` (GEN.1361), with the cut of the zone the body formed in (`cut_temperature_k` of its formation distance and the star's luminosity; a moon takes its parent's; a body with no star or distance takes the Earth's cut). Ti, Cr and Ni are kept whole near the Earth; Na and K a quarter; Mn and P about 0.6.
3. **The core.** The share in the metal is `D c / (D c + 1 - c)`, c the core's share of the dry mass. D is set so the Earth's own core (c = 0.325, FeO 0.0805) has the Earth's share: Cr 0.55, Mn 0.60, Ni 0.93, P 0.94 (the concordance of [R4608]). Ti has none; Na and K have a fixed 5 percent. D falls as the rock's FeO rises, `FeO^(-n/2)` with the valence n of the oxide (Ni, Cr and Mn 2, P 5). The core fraction itself includes the metal the elements add, so `dry_partition` iterates (three passes) from the Earth's core fraction.
4. **The mantle.** What the core did not take is the rock's oxide.

The Cr, Mn, P, Na and K the core took, and the Ni it took, join the metal; a carbon rock has none of the seven.

## What is stored

Seven columns on `planets`, `moons` and `rogue_planets` (migration 0230): `rock_na2o_fraction`, `rock_k2o_fraction`, `rock_tio2_fraction`, `rock_cr2o3_fraction`, `rock_mno_fraction`, `rock_nio_fraction`, `rock_p2o5_fraction`, each a mass fraction of the rock, NULL for a carbon rock or a body saved before v227. `surface_composition._rock_oxides_of` hands them to `rock_types.rock_oxides`, so the norm and the ground types read them.

## Calibration (a test point, never an input)

| Oxide | Model Earth (wt %) | Bulk silicate Earth | Model Mars | Mars [R4833] |
|---|---|---|---|---|
| Na2O | 0.317 | 0.35 | 0.76 | 0.59 |
| K2O | 0.041 | 0.0285 | 0.097 | 0.043 |
| TiO2 | 0.185 | 0.196 | 0.169 | 0.17 |
| Cr2O3 | 0.364 | 0.377 | 0.51 | 0.88 |
| MnO | 0.115 | 0.134 | 0.28 | 0.37 |
| NiO | 0.255 | 0.242 | 0.51 | 0.046 |
| P2O5 | 0.023 | 0.0189 | 0.13 | 0.17 |

The Earth is within 1.5 times in every oxide (the anchor of D). Mars (formation distance 1.524 AU, cut 1066 K, FeO 0.147) is within 2.5 times in all but its nickel. Mars's mantle nickel is the open item: the leading-term redox law cannot reach it, and the pressure, temperature and sulfur of the real metal-silicate partition are not in the model. Mercury and the Moon are not targets (their history is an impact, `volatile_depletion`).

Stand-ins: MO1 to MO4 in section 6.216 of the gap-filling plan.


## The phosphorus budget (GEN.1367)

`physics/phosphorus_budget.py` takes the P2O5 above to the ocean. The crust is the batch melt of the mantle, so its P2O5 is the mantle's times `1 / F` (phosphorus is incompatible), with F the melt fraction of the body's heat flow (`surface_composition.melt_fraction`); the Earth's 0.021 percent makes 0.14. A logistic of the mantle redox about iron-wuestite minus 3 decides how much of it is schreibersite; the rest is apatite where the crust is wet and merrillite where it is dry. The ocean's phosphate is the Earth's 2.3 umol/L times the ratio of the crust's P2O5, the reactive share of its hosts (merrillite weathers three times faster than apatite, schreibersite barely), the erosion of the lid and the rock area the rain reaches, over the volume of the ocean. A soda lake lacks the calcium to precipitate apatite and holds 400 times longer, to a ceiling of 1000 umol/L; a sealed ocean is cut off from rock and holds 0.01. Left out: the burial with iron oxides (inside the Earth's calibration) and the recycling of a biosphere. Stand-ins PB1 to PB10 are in section 6.221 of the gap-filling plan.
