# Rogue planet review: what is physics, what is a stand-in (GEN.1109, R1)

Planets 1 lane, 10-30. Boss (4:54 PM EDT 10-30): "physics is the driver for rogue planets". This reviews every quantity the rogue-planet generator uses, says whether each is physics, a sourced stand-in or a flat rate, names the physical replacement, and checks how many bodies the planet architecture actually ejects against microlensing. It is the R1 item of [interstellar-asteroid-fields.md](../systems/interstellar-asteroid-fields.md) section 13; R2 (GEN.1111) builds the ejection-made rogues and R3 (GEN.1112) replaces the defaults. Nothing in the generator changes in this report.

Source marks: **[read]** a paper whose text was opened for this report or for GEN.1089 and the oddity lanes; **[R]** recalled, to be confirmed by the replacement item before it is relied on.

## Result in one paragraph

The rogue-planet count is a flat rate (5.8 per star, `ROGUE_PLANET_RATE_PER_STAR`), not an output of any system. The architecture removes about **0.3 bodies per star** that could count as ejected (about 0.22 thrown by a crowded giant plus about 0.08 swept by a hot Jupiter, together roughly 9 Earth masses per star, nearly all of it gas giants), against a flat rate that gives about 2.7 per star over the range microlensing probes (0.33 to 6,660 Earth masses) and a microlensing central value of 21 per star over that range. So the architecture alone is about nine times short of our own flat rate and about seventy times short of the survey central value. The gap is not a reason to keep the flat rate. It says where the missing ejection lives: the architecture throws out almost no small bodies, because it deletes sub-minimum bodies (about 125 per star in the oligarch chain) and merges rocky pairs (losing 10 percent of mass each time, about 2.3 Earth masses per star) rather than scattering them. Barclay et al. 2017 [read] find that planets scattered by giants alone eject about 2.5 terrestrial-mass bodies per star, which is the right scale. The physical build is therefore: a scatter step (GEN.1089) that classifies each loss as ejection or collision by the Safronov number, a tally of embryos and planetesimals the architecture currently deletes, and cluster unbinding (GEN.291). The flat rate then becomes a test target.

## 1. What the architecture ejects today (measured)

Method: 1,500 stars drawn from a Kroupa-like mass function (0.1 to 1.5 solar masses), [Fe/H] normal(-0.1, 0.2), snow line 2.7 M^2 AU, disc edge 200 AU, `planet_architecture.build`, one fixed seed per star; the three loss steps wrapped to count what they remove (throwaway script, not committed; numbers will move with later architecture epochs).

| Quantity per star | Value |
|---|---|
| Oligarch chain bodies before thinning | 134.7 (26.4 Earth masses of cores) |
| After `sculpt` merges | 130.3 (24.1 Earth masses); loss 2.3 Earth masses, the 10 percent per merge of `COLLISION_MASS_RETENTION` |
| Dropped as below `PLANET_MIN_MASS_EARTH` (0.1 inner, 0.6016 outer) | about 125 bodies; not tallied as mass or as ejected |
| Planets kept | 5.0 (60.8 Earth masses including gas) |
| Thrown out by `_scatter_crowded_giants` (the lighter of a crowded giant pair) | 0.22 bodies, 6.2 Earth masses |
| Swept by `_clear_around_giants` (inside a hot Jupiter or within the clearing spacing of a moved giant) | about 0.08 bodies, about 3.3 Earth masses (400-star pass) |
| Giants per star | 0.045 (0.1 to 0.3 solar masses), 0.11 (0.3 to 0.7), 0.16 (0.7 to 1.5) |

Reading: the architecture does have the right mechanism (giants scatter and sweep) but it only fires on the 5 to 16 percent of stars with a giant, it never records speed or fate, and a lost body is deleted, not ejected, bound to a star, or collided. It has no Safronov-number test: Theta = (M_p / M_star)(a / R_p) is about 0.07 for an Earth at 1 AU and about 10 for Jupiter at 5.2 AU, so a scatterer ejects only above Theta of about 1; below that, scattering circulates the orbit or ends in collision. That test is the physics that splits "ejected" from "merged".

## 2. Microlensing checks (test targets)

| Source | Result | Our generator |
|---|---|---|
| Sumi et al. 2023, arXiv 2303.08280 [read] | dN/dlogM = 2.18 (+1.0 / -1.4) per dex per star at M/8 Earth masses, slope -0.96; total 21 (+23 / -13) per star for 0.33 to 6,660 Earth masses | flat bins give about 2.7 per star over the same range (about 8 times below central, near the lower 1 sigma of 8); slope 0.65 against 0.96 |
| Mroz et al. 2017, arXiv 1707.07634 [read] | upper limit 0.25 Jupiter-mass free-floating planets per main-sequence star | Jupiter bin 0.025 per star, under the limit |
| Mroz et al. 2020, arXiv 2009.12377 [read] | an Earth-to-Neptune population needs core-accretion ejection of 0.3 to 1 Earth mass bodies | our terrestrial bin is 4.98 per star over 0.1 to 2 Earth masses |
| Johnson et al. 2020, arXiv 2006.10760 [read] | Roman rates are forecast assuming the Cassan et al. 2012 function | the 5.8 and its bins are Boss's 10-01 research; the code comment says 6.5 (see section 4) |
| Veras and Raymond 2012, arXiv 1201.2175 [read] | planet-planet scattering alone cannot supply the giant free-floating planets | agrees with a small Jupiter bin |
| Barclay et al. 2017, arXiv 1704.08749 [read] | about 2.5 terrestrial-mass planets ejected per star when Jupiter and Saturn analogs scatter; about one third of a 5 Earth-mass disc | the scale the scatter step must reproduce |
| Ma et al. 2016, arXiv 1605.08556 [read] | core accretion gives a free-floating event rate about 1.8e-3 of stellar, about 13 times below Sumi 2011 | agrees that ejection alone is short |

Conclusion for tests: the physical build should land inside Sumi's lower 1 sigma (about 8 per star over 0.33 to 6,660 Earth masses) and under Mroz's Jupiter limit. A yield much below 1 per star means the scatter step is missing small-body ejection; a yield above 40 means it is over-ejecting.

## 3. Every rogue quantity: status and physical replacement

Kinds: **physics** (a law or a sourced fit, keep), **sourced stand-in** (a published number used as a default), **flat rate / default** (nothing behind it), **fitted** (tuned by hand).

| # | Quantity | Today | Kind | Physical replacement | Test target |
|---|---|---|---|---|---|
| 1 | Count per star | `ROGUE_PLANET_RATE_PER_STAR` 5.8, same everywhere | flat rate | sum of architecture ejections plus the scatter step (GEN.1089) plus embryo ejection plus cluster unbinding (GEN.291), by star mass, [Fe/H] and age | Sumi 2023 total 21 (+23/-13) per star over 0.33 to 6,660 Earth masses |
| 2 | Count per volume | rate per star times `REFERENCE_STELLAR_DENSITY_PC3` 0.14 (`PHENOMENON_DENSITY_PC3`) | derived from 1 | per sector, the sum over its stars' donations (storage sector decision open) | per-pc^3 figure unchanged as a check |
| 3 | Mass function | `dN/dlogM ~ M^-0.65`, 0.1 Earth to 13 Jupiter | fitted slope, hand bin edges | the ejected population's own mass function (light bodies are ejected more easily; Theta of the scatterer sets the mass cut) | Sumi slope -0.96 (+0.4 / -0.6 range) |
| 4 | Lower mass edge | 0.1 Earth mass | stand-in | ejection of the oligarch chain below `PLANET_MIN_MASS_EARTH`, the debris class below it | none |
| 5 | Gas-giant threshold | `ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER` 0.05 | rule of thumb | envelope fraction of core mass and gas capture; the ejected body keeps its architecture class | none |
| 6 | Ejection survival of air | hydrogen envelope chance 0.1 (terrestrial) / 0.5 (sub-Neptune); pressure 10 to 1,000 / 100 to 10,000 bar | default | the planet's formation outcome (gas_outcome, nebula capture) and loss by photoevaporation before ejection | Barclay: bodies keep formation atmospheres only if ejected early |
| 7 | Ejection survival of moons | rocks "that kept a moon" draw tidal flux 1e-4 to 1e-2 W/m^2 log-uniform | default | whether a moon's orbit sat inside the ejection-encounter Hill fraction [R: Debes and Sigurdsson 2007; Hong et al. 2018]; the moon then keeps its architecture orbit | share of ejected planets keeping a moon |
| 8 | Water-rich | chance 0.35 / 0.8, water share 1e-4 to 0.1 log-uniform | default | formation zone against the system's stored formation snow line | none |
| 9 | Frozen atmosphere | mass cut `ROGUE_FROZEN_ATMOSPHERE_MIN_MASS_EARTH` 0.3 | default | outgassing and volatile inventory from composition and impact history | none |
| 10 | Age | `ROGUE_PLANET_AGE_RANGE_GY` 0.1 to 10, uniform | default | age of the ejecting system at ejection plus time since, from the galaxy star-formation history | age distribution follows the host population |
| 11 | Temperature | heat law anchors: Jupiter flux 5.4 W/m^2, Neptune 53 K, cooling t^-1.3, photospheres 0.3 and 0.1 bar, cap 2,800 K | sourced fits | keep, but feed the real age (row 10) | none |
| 12 | Internal heat | radiogenic offsets (sigma 0.22 dex for Th/U, K range, 26Al factor 0.1 to 10) and `ROGUE_ACTIVE_HEAT_FLUX_W_M2` 0.04 | placeholders | the ejecting star's stored composition (GEN.476) and the body's own mantle/core split | none |
| 13 | Retention test | `ROGUE_MIN_ESCAPE_PARAMETER` 15 | physics (Jeans) | keep | none |
| 14 | Speed | generic galactic orbit; no ejection kick | missing | ejection speed from the encounter (scatterer escape speed at periapsis) added to the host's velocity (GEN.1111 storage research) | velocity dispersion of the free-floating population |
| 15 | Origin / parent | none | missing | record ejecting system and ejection time (or "formed free") so the page gets a history sentence | every rogue has a parent id or "formed free" |
| 16 | Cluster unbinding | not modelled | missing | stellar-encounter stripping of wide orbits in the birth cluster, GEN.291 | fraction unbound by cluster age |

## 4. Small findings in code

* `ROGUE_PLANET_RATE_PER_STAR` is 5.8 but the comment on `PHENOMENON_DENSITY_PC3` (tuning.py near line 3568) and the block comment above it both say 6.5 and the terrestrial share 5.6. The bins sum to 5.8 and give 4.98 terrestrial. The R3 pass fixes the wording when it replaces the rate.
* The terrestrial bin holds 86 percent of rogues by count, so the rate and the lower edge (rows 1 and 4) decide almost everything; a physical replacement of row 4 moves the total the most.

## 5. Items proposed (for TODO to file)

* **RP1 (GEN.1089 extension, non-lock):** classify every architecture loss (thinning, sub-minimum drop, giant scatter, hot-Jupiter sweep) by Safronov number into ejected, collided or bound, and return an ejection record (mass, composition, orbit, scatterer, age); tests against Barclay's 2.5 terrestrial bodies per star and Theta thresholds.
* **RP2 (non-lock):** ejection mass-function module (light-body ejection bias, Theta cut) producing `dN/dlogM`, test against Sumi slope -0.96.
* **RP3 (non-lock):** photoevaporation and volatile loss before ejection (rows 6, 8, 9) from the formation outcome.
* **RP4 (non-lock):** moon survival through the ejection encounter (row 7), Hill-fraction criterion.
* **RP5 (non-lock):** age and time-since-ejection from star-formation history (row 10) feeding the heat law (row 11).
* **RP6 (non-lock):** ejection speed (row 14) from escape speed at periapsis plus host velocity; interstellar-object speed window.
* **RP7 (non-lock):** cluster unbinding (GEN.291, row 16), a rate per cluster age and density.
* **RP8 (GEN.1111, lock):** create the rogue rows from RP1 to RP6, storage sector per the research decision, parent id and rank in the object id; the flat rate stays only as a test target.
* **RP9 (non-lock):** a calibration test of the full yield against Sumi 2023, Mroz 2017 and Barclay 2017, run after RP8.

R3 (GEN.1112) is RP2 to RP7 as the replacement modules; RP1 is the first of them because every other row needs its ejection record.

## 6. Stand-ins this review registers

None added to GEN.725. The audit adds no number to the generator; rows 1 to 12 and 14 to 16 already count as stand-ins in the design table and are replaced by RP1 to RP7.
