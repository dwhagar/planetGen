# Interstellar object rates

Reference numbers for TODO items 5 (rates) and 6 (rogue planet masses),
from the research Boss supplied on 2026-09-30. Built in schema v37
(`program_constants.PHENOMENON_DENSITY_PC3`, `phenomenon_rate_per_star`,
`ROGUE_PLANET_MASS_BINS`; `generate.generate_sector_phenomena` and
`flag_fast_stars`), at the full rates Boss chose, with the table's
hypervelocity figure replaced by 1e-10 (see "Checks"). Molecular clouds
(filling factor) and the other nebula kinds wait for TODO item 27; until
then the only generated nebulae are planetary ones. The figures are Boss's
chosen reference. The ones that don't add up are flagged under "Checks".

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
| Isolated neutron stars | 1e-3 | ~7e-3 | n_* | ~0.064 |
| Isolated stellar black holes | 1e-4 | ~7e-4 | n_* | ~6e-3 |
| Giant molecular clouds | 1e-6 to 1e-5 (5e-6) | n/a | ρ_gas^1.4 | ~3e-4 |
| Planetary nebulae | 3e-8 | ~2e-7 | n_* | ~2e-6 |
| Supernova remnants | 1e-8 to 1e-7 | n/a | n_* · ρ_gas | ~6e-7 |
| Hypervelocity stars (> 500 km/s) | 5e-9 at 8 kpc (see Checks) | n/a | r_GC⁻² | ~3e-7 |
| Isolated asteroid fields | ~0 (disperse in 1e6-1e7 yr) | 0 | n/a | 0 |

\* A 4 pc sector is 64 pc³, about 9 stars at local density. The generator
makes about 10 systems per default sector today.

Galaxy totals behind these (for scale): 25-100 billion brown dwarfs,
1e7-1e8 runaway stars, 1e3-1e4 hypervelocity stars, ~1e9 isolated neutron
stars, ~1e8 isolated black holes, ~20,000 planetary nebulae, ~1,000
detectable supernova remnants. Intracluster ("hostless") supernovae are
5-20% of supernovae in rich galaxy clusters; they happen between
galaxies, so they have no place in one galaxy's generator.

## Rogue planet mass bins (item 6)

Low-mass rogues dominate: disk scattering ejects small bodies while
giants stay bound. Draw a bin by its per-star rate, then a mass
log-uniformly inside it.

| Bin | Mass | Per star | Source |
|---|---|---|---|
| Terrestrial | 0.1-2 M⊕ | 5 (range 2-10) | Johnson et al. 2020; Mróz et al. 2020 (OGLE-2016-BLG-1928, 0.3-2 M⊕) |
| Sub-Neptune / ice giant | 2-20 M⊕ | 1 (default, not in the research) | Sumi et al. 2023 find Neptune-mass candidates; no rate given |
| Saturn-class | 20 M⊕-1 Mjup | 0.25 (default, not in the research) | fills the gap between the two constrained bins |
| Jupiter-mass | 1-13 Mjup | 0.25 (upper limit) | Mróz et al. 2017; replaces Sumi et al. 2011's 1.8 per star |

Above 13 Mjup is a brown dwarf, its own row in the density table.

A single power law doesn't fit both ends: `dN/dlogM ∝ M^-0.96` (the
Sumi 2023 slope) normalized to 0.25 per star above 1 Mjup gives about
575 per star above 0.1 M⊕, a hundred times the terrestrial estimate.
Hence the bins.

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
  Today's generator uses 2,000 galaxy-wide. Pick 1e-8 only if faint,
  undetectable remnants count too.
- Terrestrial (0.7 pc⁻³ = 5 per star), Jupiter-mass (0.035 = 0.25 per
  star), brown dwarfs (0.03 = 1 per 4.7 stars), runaways (1.5%), neutron
  stars and black holes (0.5-0.7% and 0.05-0.07% of stars, matching the
  galaxy totals) are self-consistent.

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
