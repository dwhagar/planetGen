# Interstellar object rates

**Status (2026-10-15):** built in 7.18.0 (schema v37, PR #136; GEN.1 and GEN.2). The GEN.100 retune (#793) changed the rates afterwards; see [compact-remnant-regions.md](compact-remnant-regions.md).

Reference numbers for GEN.1 (rates) and GEN.2 (rogue planet masses),
from the research Boss supplied on 2026-09-30. Built in version 7.18.0
(schema v37, PR #136): `program_constants.PHENOMENON_DENSITY_PC3`,
`PHENOMENON_RATE_SCALE`, `phenomenon_rate_per_star`,
`ROGUE_PLANET_MASS_BINS`, `generate.generate_sector_phenomena` and
`flag_fast_stars`. Boss chose the full rates, with the table's
hypervelocity figure replaced by 1e-10 (see "Checks"). The figures are
Boss's chosen reference. The ones that don't add up are flagged under
"Checks". On the pages, rogue planet mass classes have reference pages
from 7.46.0 (`/classes/rogue-planet`), and rogue planets, interstellar
comets, black holes and neutron stars get rendered views from 7.57.0.

Nebulae (GEN.10) followed in 7.30.0: molecular clouds were drawn per star at
the `molecular-cloud` rate below and generated as dark-family classes M-Q
(correction, GEN.47: in a galaxy they now come from the seeded
`nebula_field` instead, see the design note);
H II regions and reflection nebulae grow around O, B and A stars
(`NEBULA_HOST_RULES`; the stars a nebula needs are made with it, GEN.320 and GEN.321, see "The stars a nebula needs" below); every planetary nebula gets its own hot white dwarf
system. Diffuse gas (classes A-B) is the background and is not generated.
See `nebula-and-asteroid-field-classes.md`; its "Corrected facts" section
lists what the 2026-10-09 research found wrong in the host rules and the
H II sizes. Anomalies beyond these kinds (magnetars and the rest) are in
`anomalies.md`. Wolf-Rayet stars have no rate in this table: they are a state
of the very massive stars the bright and mass passes already place
(`physics/wolf_rayet.py`, GEN.243), so their number follows from the star
formation history and the metallicity, not from a density of their own; see
"Wolf-Rayet stars and the ring nebula" below.

## How the code uses these numbers

- Every kind in `PHENOMENON_DENSITY_PC3` becomes a rate per star:
  `density / 0.14 * PHENOMENON_RATE_SCALE[kind]`
  (`phenomenon_rate_per_star`). All scales are 1.0.
- `generate_sector_phenomena` counts the sector's stars (a binary counts
  two, `sector_star_count`) and draws each kind's count from a Poisson
  distribution with mean `rate * stars`. This is the "per star" form the
  model below asks for.
- A kind's rate per star is scaled by the place it is rolled in: the sector
  profile's `kind_factors` ("Rates by place" below, GEN.236). The mean number of
  a kind in a sector is `rate * stars * factor`; a kind with no entry in
  `kind_factors` is scaled by star count only. (The research's gas-density
  scaling for supernova remnants, `n_* * rho_gas`, is the `supernova-remnant`
  factor.) Correction (GEN.47,
  checked against the code 2026-10-09): molecular clouds in a galaxy are
  now scaled by gas density, `GMC_GAS_DENSITY_EXPONENT` is read by
  `galaxy/nebula_field.py`; `GMC_ARM_FILLING_FACTOR` is still defined and
  read by nothing.
- Correction (GEN.100): in a galaxy, black holes, neutron stars, planetary
  nebulae, supernova remnants and hypervelocity stars are not rolled per
  sector. They are drawn galaxy-wide at the same per-star rates and placed
  by `generation/phenomenon_scatter.py` right after the plan; a sector's
  fill only builds what was placed there.
- Runaway and hypervelocity stars are flags on ordinary systems
  (`star_systems.runaway_class`), not phenomenon rows. The hypervelocity
  chance scales as `(8000 pc / r_GC)^2`.
- Interstellar comets kept as rows use a design rate, `"comet": 0.007`
  per pc³ (0.05 per star), not the research density. The real population
  (1e15 per pc³, `INTERSTELLAR_DEBRIS_DENSITY_PC3`) is shown only as a
  computed count on the sector page (`queryDb.interstellar_debris_count`).
- Asteroid fields have rate 0; the type stays for hand-made fields and
  facilities.

## The model

Treat each object type as a Poisson process in space. Its local number
density `n_i` (objects per pc³) scales with the local stellar density
`n_*`, normalized to the solar neighborhood value `n_*0 ≈ 0.14 pc⁻³`
(0.10-0.14 in the literature):

    n_i(r) = n_i0 * (n_*(r) / 0.14)

The chance of at least one object of type i in a volume V is
`1 - exp(-n_i V)`, and the count is Poisson with mean `n_i V`.

Because a sector's system count already follows `n_* V`, a density that
scales with `n_*` is the same as a fixed rate **per star**:
`n_i0 / 0.14`. So the existing per-sector Poisson draw in
`generate.generate_sector_phenomena` keeps working; only the rates change
(and "per star system" should become "per star", or be divided by the
mean stars per system, about 1.3 with binaries).

Three types don't scale with `n_*` alone:

- Giant molecular clouds scale with gas density, `∝ ρ_gas^1.4`
  (Schmidt-Kennicutt), and are extended: a filling factor of about
  0.01-0.02 inside spiral arms.
- Supernova remnants scale with `n_* · ρ_gas`.
- Hypervelocity stars scale with `r_GC⁻²` (distance from the galactic
  center), since Sagittarius A* ejects them.

## Densities (solar neighborhood, n_*0 = 0.14 pc⁻³)

| Type | n at R0 (pc⁻³) | Per star | Scales with | Per 4 pc sector* |
|---|---|---|---|---|
| Interstellar comets and debris (about 100 m and up) | 7e14 to 1.8e15 (1e15 central; 0.1 to 0.2 per au³) | ~7e15 | n_* | ~6e16 |
| Terrestrial rogue planets (0.1-2 M⊕) | 0.5-1.4 (0.7 central) | 2-10 (5) | n_* | ~45 |
| Jupiter-mass rogue planets (> 1 Mjup) | 0.035 | ≤ 0.25 | n_* | ~2.2 |
| Rogue brown dwarfs (13-80 Mjup) | 0.025-0.035 (0.03) | ~0.21 | n_* | ~1.9 |
| Runaway stars (> 30 km/s) | 2.1e-3 | 0.015 (1-2%) | n_* | ~0.13 |
| Isolated neutron stars | 5.6e-4 (1e9 over ~2.5e11 stars) | 4e-3 | n_* | ~0.036 |
| Isolated stellar black holes | 7e-5 (1e8 over ~2.5e11 stars) | 5e-4 | n_* | ~4.5e-3 |
| Giant molecular clouds | 1e-6 to 1e-5 (5e-6) | n/a | ρ_gas^1.4 | ~3e-4 (see the cloud-density check below) |
| Planetary nebulae | 1.4e-8 | ~1e-7 | n_* | ~9e-7 |
| Supernova remnants | 1e-8 to 1e-7 | n/a | n_* · ρ_gas | ~6e-7 |
| Hypervelocity stars (> 500 km/s) | 5e-9 at 8 kpc (see Checks) | n/a | r_GC⁻² | ~3e-7 |
| Isolated asteroid fields | ~0 (disperse in 1e6-1e7 yr) | 0 | n/a | 0 |

\* A 4 pc sector is 64 pc³, about 9 stars at local density. The density
model expects 6.3 systems in a local-density sector (`galaxy_shape.
expected_system_count_at_density_1`); with binaries that is about 8 stars.

Galaxy totals behind these (for scale): 25-100 billion brown dwarfs,
1e7-1e8 runaway stars, 1e3-1e4 hypervelocity stars, ~1e9 isolated neutron
stars, ~1e8 isolated black holes, ~20,000 planetary nebulae, ~1,000
detectable supernova remnants. Intracluster ("hostless") supernovae are
5-20% of supernovae in rich galaxy clusters; they happen between
galaxies, so they have no place in one galaxy's generator.

## Rogue planet mass bins (GEN.2, GEN.45)

Low-mass rogues dominate: disk scattering ejects small bodies while
giants stay bound. Draw a bin by its per-star rate, then a mass inside it.

Since GEN.45 both the rates and the draw inside a bin follow one mass
function, `dN/dlogM ∝ M^-0.65` from 0.1 M⊕ to 13 Mjup (Boss's research of
2026-10-01), normalized to the same 6.5 rogues per star
(`ROGUE_PLANET_MASS_FUNCTION_SLOPE`, `ROGUE_PLANET_RATE_PER_STAR`). The
research writes it `dN/dM ∝ M^-0.65`; read per unit mass, that would make
about 87% of rogues gas giants, against its own aim of keeping them rare,
so it is read per logarithm of mass.

(Pre-retune figures; the GEN.100 retune of 2026-10-09 sets terrestrial rogue planets to 5.8 per star.)

| Bin | Mass | Per star | Check against |
|---|---|---|---|
| Terrestrial | 0.1-2 M⊕ | 5.6 | 2-10: [R2178, R3190] (OGLE-2016-BLG-1928, 0.3-2 M⊕) |
| Sub-Neptune / ice giant | 2-20 M⊕ | 0.72 | [R4254, R4255] find Neptune-mass candidates; no rate given |
| Saturn-class | 20 M⊕-1 Mjup | 0.17 | no rate in the research |
| Jupiter-mass | 1-13 Mjup | 0.028 | un[R3191]'s upper limit of 0.25 |

That is about 95% terrestrial and 5% giants (past
`ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER`, 10 M⊕, the core mass that holds an envelope; GEN.1133 moved it from 16). The
hand-set rates before it (5, 1, 0.25 and 0.25 per star, a log-uniform
mass in each bin) gave 91% and 9%. Above 13 Mjup is a brown dwarf, its
own row in the density table.

A steeper law doesn't fit both ends: `dN/dlogM ∝ M^-0.96` (the [R4254, R4255]
slope) normalized to 0.25 per star above 1 Mjup gives about 575 per star
above 0.1 M⊕, a hundred times the terrestrial estimate.

## Checks

- **Comet separation column (GEN.659).** The old research table gave 1e11-1e14 pc⁻³ (1e-4 to 1e-2 AU⁻³, a spacing of
  45 to 4.5 AU). The papers say more: Do, [R1124] find about 0.2 AU⁻³ for objects like '[R2159] about 0.1 AU⁻³ (1 pc³ = 8.8e15 AU³), so 7e14 to 1.8e15 pc⁻³, a spacing `n^(-1/3)` of about
  2.2 to 1.7 AU. The constant takes 1e15, the middle of that range, for objects of about 100 m and up.
- **Hypervelocity stars.** With `n ∝ r_GC⁻²`, the count inside radius R
  is `4π n0 r0² R`. Using 5e-9 pc⁻³ at 8 kpc gives about 30,000 inside
  8 kpc and about 400,000 out to 100 kpc, against the stated 1e3-1e4
  total. A total of 1e4 out to 100 kpc implies about 1e-10 pc⁻³ at 8 kpc.
  Use 1e-10 unless Boss prefers the table's figure.
- **Supernova remnants.** ~1,000 detectable remnants over a disk of about
  7e11 pc³ is about 1.5e-9 pc⁻³; the table's 1e-8 implies about 7,000.
  The generator before 7.18.0 used 2,000 galaxy-wide. The code now uses
  1e-8, which counts faint, undetectable remnants too.
- **Molecular cloud density (recommendation, 2026-10-09).** 5e-6 pc⁻³ is
  the top of the table's range and, with the GEN.47 field, puts 3% to 46%
  of the volume inside a cloud (13% of arm sectors at the solar circle)
  against about 0.5% to 1% observed; the catalogues (about 8,100 to 9,700
  clouds, 1,064 massive ones) give 1.2e-8 to 1.4e-8 pc⁻³. About 1e-7 pc⁻³
  for class M, with the small dark classes as their own rows, would match.
  Boss chose the 5e-6 research value, so this is a recommendation only
  (derivation in `nebula-and-asteroid-field-classes.md`, "The molecular
  cloud field").
- **Neutron star pulsar fraction (recommendation, 2026-10-09; pre-retune: the GEN.100 retune of 2026-10-09 sets 2% of neutron stars pulsing, and `NEUTRON_STAR_PULSAR_CHANCE` no longer exists in `src/`).** The
  catalogue above gives about 5e8 isolated neutron stars.
  `NEUTRON_STAR_PULSAR_CHANCE = 0.7` makes 3.5e8 of them pulsars, against roughly 1e5
  beamed radio pulsars (an active fraction of about 1e-4 to 1e-3), because
  radio emission ends after about 1e7 years and an isolated star has no
  companion to recycle it. The constant is a gameplay choice in the code
  today; the physical value and an age rule are in `anomalies.md`.
- (Pre-retune values; the GEN.100 retune of 2026-10-09 sets terrestrial rogue planets to 5.8 per star, neutron stars to 0.4% and black holes to 0.05%, see "Retune of 2026-10-09" below.) Terrestrial (0.78 pc⁻³ = 5.6 per star), Jupiter-mass (0.004 = 0.028
  per star, under the 0.035 limit), brown dwarfs (0.03 = 1 per 4.7 stars), runaways (1.5%), neutron
  stars and black holes (0.5-0.7% and 0.05-0.07% of stars, matching the
  galaxy totals) are self-consistent. The generator used 0.5% and 0.1%
  before the retune of 2026-10-09 instead, the shares the mass-and-age star model leaves behind (star-fix
  study, 2026-09-30), so isolated remnants and the star census agree.

## Generation cost

Measured on main at 2026-09-30 (default args, 10 sectors): a sector
builds about 10 systems in about 0.1 s, roughly 650 body rows (stars,
planets, moons). Building one rogue planet object takes about 0.13 ms.

At the research rates a local-density sector adds about 50 rogue rows
(45 terrestrial, 2 giant, 2 brown dwarfs, plus ~1 sub-Neptune class per
star under the default bins, so nearer 60). That is about 8-10% more rows
and well under 1% more generation time. Storage and speed are not the
problem. What changes is what a sector looks like: rogues would
outnumber systems about five to one in the Contents table and on the
Sector Map.

Comets can't be rows: ~6e16 per sector. Neutron stars and everything
rarer are fine as rows; most sectors simply get none.

## References (short form)

Rogue planets: [R3191]
(Nature 548:183); [R2178] (AJ 160:123); [R3190]
(ApJL 903:L11); [R4254, R4255] (AJ
166:109); [R2717, R2720] (ApJL 777:L20); [R3095] (Nat.
Astron. 6:89); [R3467].

Other objects: [R2374] (ApJS 253:7) and Mužić et al.
2017 (MNRAS 471:3699), brown dwarfs; [R4310] (MNRAS 448:L6), runaways
and kicks; [R0556] (ARA&A 53:15) and [R2444] (MNRAS
491:2465), hypervelocity stars; [R0040] (MNRAS 334:553),
[R3348] (A&A 638:A94), [R2538] (ApJL 933:L23) and Sahu
et al. 2022 (ApJ 933:83), isolated compact remnants; Engelhardt et al.
2017 (AJ 153:133), [R4005] (ApJL 896:L8) and Raymond et
al. 2020 (ApJL 894:L22), interstellar comets and the absence of free
asteroid fields; [R2320] (ARA&A 50:531),
nebulae and remnants; [R3907] and McGee &, intracluster supernovae.

## Why it works this way

- **Per star, not per system or per volume.** A sector's system count
  already follows local stellar density, so a fixed rate per star gives the
  density scaling for free and keeps the old per-sector Poisson draw. A
  binary counts as two stars, which is why the rate is per star, not per
  system (section "The model").
- **Full rates.** Boss chose the research rates at full strength on
  2026-09-30. `PHENOMENON_RATE_SCALE` exists so a kind can be dialed down
  (for example 0.1 for rogue planets) without changing its research value.
- **Mass bins for rogue planets.** One power law cannot match both the
  terrestrial and the Jupiter-mass measurements (it overshoots the
  terrestrial count a hundredfold), so the rates are binned
  ("Rogue planet mass bins").
- **Hypervelocity 1e-10.** The table's 5e-9 implies 30 to 40 times more
  hypervelocity stars than its own galaxy total ("Checks").
- **Neutron stars 0.5% and black holes 0.1% of stars (pre-retune; 0.4% and 0.05% after the GEN.100 retune of 2026-10-09, and the black hole density has been metallicity-dependent since GEN.563, black-hole-density.md).** These were the shares
  the mass-and-age star model (7.23.0) leaves behind, so isolated remnants
  and the star census agree.
- **No free asteroid fields.** They disperse in 1e6 to 1e7 years.
- **Comets as a figure, not rows.** About 6e16 per sector cannot be stored.

### Alternatives not taken

- Scaling clouds and remnants by gas density. The constants are in place
  for it, but the reason it was not wired up is not recorded. (Clouds were
  wired up later by GEN.47; remnants still are not.)
- Supernova remnants at the detectable-only density (1.5e-9). The code
  chose 1e-8; the comment says it "counts faint remnants too". No further
  reason is recorded.


## Rates by place: the populations and the kind factors (GEN.236)

The per-star rates above are the solar neighbourhood's. Every sector carries a profile (`galaxy/profile.py`, GEN.235; see
[intergalactic-and-halo-contents.md](../galaxy/intergalactic-and-halo-contents.md) section 5.3) that says which stellar
populations it holds, and `kind_factors_at` turns those populations and the gas into one factor per kind, each in 0 to 1,
that multiplies the kind's rate per star. A halo cell, a cluster or an old satellite therefore makes no H II region and no
core-collapse remnant, and the disk's young, gas-rich places make the most.

The populations (`profile.POPULATIONS`): `young`, `intermediate` and `old` (the thin disk by age, with the young open
clusters of `tuning.OPEN_CLUSTER_POPULATIONS` counted with the disk), `bulge`, `halo`, `cluster`, `satellite` (with
`satellite_young` and `satellite_intermediate`), `nuclear_cluster` and `nuclear_disc`. Each has an age window, [Fe/H]
peaks, a velocity dispersion and a rotation rule.

| Kind | Factor |
|---|---|
| `brown-dwarf` | 1 everywhere (every population makes them) |
| `black-hole` (stellar mass) | the share of the place's massive stars that failed to explode and made their remnant, from the populations' shares, ages and [Fe/H] peaks (`remnant_density.place_factor_from_peaks`, GEN.563): a metal-poor place makes more |
| `neutron-star` | what the black holes leave of the massive stars that made their remnant, less the galaxy's unbound share (GEN.554): a metal-rich place makes more |
| `supernova-remnant` | the gas fraction (`nebula_field.gas_factor`, capped at 1), times the share of stars that are not in a cluster; 0 in the halo and in clusters. Not zeroed in a superbubble: the hot cavity is where the explosions went off. The Type Ia share of the remnants is a profile field, `type_ia_share` (GEN.319): 1 where no young population is |
| `planetary-nebula` | the share of the place's stars that are not young against the solar circle's, capped at 1 (GEN.319): 1 in the disk (within 3%) and in every older place |
| `molecular-cloud` | the cold gas fraction times the young share against the solar circle's |
| `nebula` (emission and reflection) | the cold gas fraction; the host star's own age is the young share here |
| `rogue-planet`, `comet` | the planet occurrence at the place's [Fe/H] against the disk's (`metallicity.planet_occurrence_factor`), averaged over the populations |

A kind not in the table is not scaled by place. The neutron stars and black holes the galaxy lost to a kick
(`ejecta_field`, GEN.238) are not part of the stars' rate: the fill adds them to the cell's own draw at a density that falls
as 1 / r^2. The galaxy-wide scatters (`generation/phenomenon_scatter.py`) draw the black holes, neutron stars, planetary
nebulae and supernova remnants with the same factors, band by band.

## The stars a nebula needs (GEN.320, GEN.321, GEN.241)

A nebula of a class that has stars gets them as real systems, inside its bounds, from seeds that depend on the nebula or the
star alone, never on the order sectors are filled in (`generation/nebula_contents.py`). Whether a hot star has a nebula, and
of which class, is rolled on a stream of the star's own (`hosted_seed`), so a star the scatter pre-placed and the same star
made in a fill agree. A planetary nebula's white dwarf is made with it, and `PLANETARY_NEBULA_COMPANION_CHANCE` of them gain a
close companion. A Bok globule (P) or a star-forming core (Q) gets its protostars from the cloud's own seed, built by whichever
sector holds each point. The cluster behind an H II region of class D or E is its ionizing members, the stars of 5.5 solar
masses and up still alive, placed in a Plummer sphere of half-mass radius 1 pc (D) or 2 pc (E) clipped to the nebula's
radius; only these are real, and the rest of the cluster stays the sector's background density. A supernova remnant keeps its
neutron star or black hole through the remnant object itself.

## Wolf-Rayet stars and the ring nebula (GEN.243)

A star of at least `tuning.WR_MIN_INITIAL_MASS_SOL` (25 solar masses at the Sun's metallicity, rising as Z^-0.5, so 35 at half
and 56 at a fifth of it) spends the last `WR_PHASE_MYR` (0.3 Myr, at most half its supergiant window) of its life as a
Wolf-Rayet star, in place of the end of the supergiant window of `stellar_evolution.evolve_star`: yerkes class `W`, a present
mass of 0.4 to 0.6 of the initial mass, a luminosity log-uniform over 1e5 to 1e6 solar, and a subtype by the share of the
phase gone (WN 60%, WC 30%, WO 10%, the last only from 40 solar masses; `WR_TEMPERATURE_RANGE_K` and `WR_SUBCLASS_RANGE`
give the temperature and the subclass, so the text reads `WN6 Blue Wolf-Rayet Star`). No Wolf-Rayet star is rolled from a rate:
the bright and mass passes place the massive stars and draw their state, which gives about 1,180 in a Milky Way-like disk
against the 1,200 expected, and none in a population older than the phase allows. Their companions follow the subtype
(`WR_COMPANION_CHANCE`: 0.52 for WN, 0.90 for WC and WO, read by the system builder; the `wolf-rayet-binary` row of the phenomena catalogue). A
Wolf-Rayet star can light a **class X ring nebula** (family `emission`, radius 6.5 to 114 ly, a thin swept-up shell around a
hollow interior, `host_only`: only a Wolf-Rayet star makes one, so a nebula drawn with no host is never X); the host
rule in `NEBULA_HOST_RULES` rolls it from the star, and `NEBULA_CLASS_MAX_AGE_MYR["X"]` caps the host's age at 10 Myr. The
rest of the design and the numbers are in the intergalactic design note, section 5.11.

## Retune of 2026-10-09 (Boss: central observed values)

See also [star-types-by-galactic-radius.md](../galaxy/star-types-by-galactic-radius.md)
(GEN.133) for how star-type shares change with radius.

Neutron stars 0.4% and black holes 0.05% of stars; planetary nebulae 1e-7 per star; terrestrial rogue planets 5.8 per star (0.7 per pc³ central); intermediate-mass black holes 0.1% of black holes; accretion disks on 0.1% of black holes; 2% of neutron stars pulsing, 10% of those millisecond pulsars. Regional differences (scale height, radius, type) are in [`compact-remnant-regions.md`](compact-remnant-regions.md) and are not modeled yet.

See also [multistar-and-compact-systems.md](multistar-and-compact-systems.md) for binaries and systems around compact remnants, and [anomalies.md](anomalies.md) for the anomaly classes and rates.
