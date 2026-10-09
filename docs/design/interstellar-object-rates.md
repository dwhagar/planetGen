# Interstellar object rates

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
(`NEBULA_HOST_RULES`); every planetary nebula gets its own hot white dwarf
system. Diffuse gas (classes A-B) is the background and is not generated.
See `nebula-and-asteroid-field-classes.md`; its "Corrected facts" section
lists what the 2026-10-09 research found wrong in the host rules and the
H II sizes. Anomalies beyond these kinds (magnetars, Wolf-Rayet stars and
the rest) are in `anomalies.md`.

## How the code uses these numbers

- Every kind in `PHENOMENON_DENSITY_PC3` becomes a rate per star:
  `density / 0.14 * PHENOMENON_RATE_SCALE[kind]`
  (`phenomenon_rate_per_star`). All scales are 1.0.
- `generate_sector_phenomena` counts the sector's stars (a binary counts
  two, `sector_star_count`) and draws each kind's count from a Poisson
  distribution with mean `rate * stars`. This is the "per star" form the
  model below asks for.
- Every kind is scaled by star count only. The gas-density scaling for
  supernova remnants (`n_* * rho_gas`) is not applied. Correction (GEN.47,
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
  (1e12 per pc³, `INTERSTELLAR_DEBRIS_DENSITY_PC3`) is shown only as a
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
| Interstellar comets and debris (m-km) | 1e11 to 1e14 (1e12 central) | ~7e12 | n_* | ~6e13 |
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

| Bin | Mass | Per star | Check against |
|---|---|---|---|
| Terrestrial | 0.1-2 M⊕ | 5.6 | 2-10: Johnson et al. 2020; Mróz et al. 2020 (OGLE-2016-BLG-1928, 0.3-2 M⊕) |
| Sub-Neptune / ice giant | 2-20 M⊕ | 0.72 | Sumi et al. 2023 find Neptune-mass candidates; no rate given |
| Saturn-class | 20 M⊕-1 Mjup | 0.17 | no rate in the research |
| Jupiter-mass | 1-13 Mjup | 0.028 | under Mróz et al. 2017's upper limit of 0.25 (replacing Sumi et al. 2011's 1.8) |

That is about 96% terrestrial and 4% gas giants (past
`ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER`, about 16 M⊕). The
hand-set rates before it (5, 1, 0.25 and 0.25 per star, a log-uniform
mass in each bin) gave 91% and 9%. Above 13 Mjup is a brown dwarf, its
own row in the density table.

A steeper law doesn't fit both ends: `dN/dlogM ∝ M^-0.96` (the Sumi 2023
slope) normalized to 0.25 per star above 1 Mjup gives about 575 per star
above 0.1 M⊕, a hundred times the terrestrial estimate.

## Checks

- **Comet separation column.** The research table gives a characteristic
  separation of 0.001-0.02 AU for 1e11-1e14 pc⁻³. `n^(-1/3)` of those is
  about 45 AU down to 4.5 AU (1 pc³ = 8.8e15 AU³). The densities
  themselves agree with the ~1e-4 to 1e-2 AU⁻³ quoted beside them.
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
- **Neutron star pulsar fraction (recommendation, 2026-10-09).** The
  catalogue above gives about 5e8 isolated neutron stars.
  `NEUTRON_STAR_PULSAR_CHANCE = 0.7` makes 3.5e8 of them pulsars, against roughly 1e5
  beamed radio pulsars (an active fraction of about 1e-4 to 1e-3), because
  radio emission ends after about 1e7 years and an isolated star has no
  companion to recycle it. The constant is a gameplay choice in the code
  today; the physical value and an age rule are in `anomalies.md`.
- Terrestrial (0.78 pc⁻³ = 5.6 per star), Jupiter-mass (0.004 = 0.028
  per star, under the 0.035 limit), brown dwarfs (0.03 = 1 per 4.7 stars), runaways (1.5%), neutron
  stars and black holes (0.5-0.7% and 0.05-0.07% of stars, matching the
  galaxy totals) are self-consistent. The generator uses 0.5% and 0.1%
  instead, the shares the mass-and-age star model leaves behind (star-fix
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

Comets can't be rows: ~6e13 per sector. Neutron stars and everything
rarer are fine as rows; most sectors simply get none.

## References (short form)

Rogue planets: Sumi et al. 2011 (Nature 473:349); Mróz et al. 2017
(Nature 548:183); Johnson et al. 2020 (AJ 160:123); Mróz et al. 2020
(ApJL 903:L11); Barclay et al. 2023 (AJ 166:162); Sumi et al. 2023 (AJ
166:109); Liu et al. 2013 (ApJL 777:L20); Miret-Roig et al. 2022 (Nat.
Astron. 6:89); Pearson & McCaughrean 2023 (arXiv:2310.01231).

Other objects: Kirkpatrick et al. 2021 (ApJS 253:7) and Mužić et al.
2017 (MNRAS 471:3699), brown dwarfs; Tauris 2015 (MNRAS 448:L6), runaways
and kicks; Brown 2015 (ARA&A 53:15) and Koposov et al. 2020 (MNRAS
491:2465), hypervelocity stars; Agol & Kamionkowski 2002 (MNRAS 334:553),
Olejak et al. 2020 (A&A 638:A94), Lam et al. 2022 (ApJL 933:L23) and Sahu
et al. 2022 (ApJ 933:83), isolated compact remnants; Engelhardt et al.
2017 (AJ 153:133), Seligman & Laughlin 2020 (ApJL 896:L8) and Raymond et
al. 2020 (ApJL 894:L22), interstellar comets and the absence of free
asteroid fields; Draine 2011 and Kennicutt & Evans 2012 (ARA&A 50:531),
nebulae and remnants; Graham et al. 2015, Sand et al. 2011 and McGee &
Balogh 2010, intracluster supernovae.

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
- **Neutron stars 0.5% and black holes 0.1% of stars.** These are the shares
  the mass-and-age star model (7.23.0) leaves behind, so isolated remnants
  and the star census agree.
- **No free asteroid fields.** They disperse in 1e6 to 1e7 years.
- **Comets as a figure, not rows.** About 6e13 per sector cannot be stored.

### Alternatives not taken

- Scaling clouds and remnants by gas density. The constants are in place
  for it, but the reason it was not wired up is not recorded. (Clouds were
  wired up later by GEN.47; remnants still are not.)
- Supernova remnants at the detectable-only density (1.5e-9). The code
  chose 1e-8; the comment says it "counts faint remnants too". No further
  reason is recorded.


## Retune of 2026-10-09 (Boss: central observed values)

See also [star-types-by-galactic-radius.md](star-types-by-galactic-radius.md)
(GEN.133) for how star-type shares change with radius.

Neutron stars 0.4% and black holes 0.05% of stars; planetary nebulae 1e-7 per star; terrestrial rogue planets 5.8 per star (0.7 per pc³ central); intermediate-mass black holes 0.1% of black holes; accretion disks on 0.1% of black holes; 2% of neutron stars pulsing, 10% of those millisecond pulsars. Regional differences (scale height, radius, type) are in [`compact-remnant-regions.md`](compact-remnant-regions.md) and are not modeled yet.

See also [multistar-and-compact-systems.md](multistar-and-compact-systems.md) for binaries and systems around compact remnants, and [anomalies.md](anomalies.md) for the anomaly classes and rates.
