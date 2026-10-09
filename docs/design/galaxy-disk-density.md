# Galaxy density model

**Status:** built (`src/planetgen/galaxy/density.py`). The model dates
from the "revision 2" design pass of 2026-09 (Track C). That revision's main
proposal, a stored plan of every qualifying sector, was never built; it is
kept in `archive/galaxy-disk-density-rev2.md`. This file describes the model
as the code uses it now, on the cylindrical grid
(`galaxy-coordinate-system.md`).

## 1. The model

A thin disk carrying the spiral arms, a thick disk and a boxy bar bulge,
each fitted to the Milky Way's published structure (GEN.118, GEN.119):

```
R     = sqrt(x^2 + y^2)
theta = atan2(y, x)
v(z, h) = sech^2(z / 2h)            # isothermal sheet; tail exp(-|z|/h)

thin(R, z)  = exp(-R / L) * v(z, h)                       # L, h: the shape's
                                                          # disk scale length/height
thick(R, z) = a_T * exp(-R / L_T) * v(z, h_T)             # L_T = 0.77 L, h_T = 3 h
bulge(x, y, z) = bulge_amplitude * exp(-r_s^2 / 2)        # in the bar's frame:
r_s^4 = ((x'/x0)^2 + (y'/y0)^2)^2 + (z/z0)^4              # x0 = bulge_scale_radius_pc,
                                                          # y0 = 0.39 x0, z0 = 0.27 x0
theta_arm(R)      = spiral_reference_angle_rad + ln(R / spiral_reference_radius_pc) / tan(pitch_angle_rad)
arm_factor(R, th) = 1 + arm_amplitude * cos(arm_count * (th - theta_arm))

relative_density = k_norm * (bulge + thin * arm_factor + thick)
```

- **Thin disk.** `disk_scale_height_pc` is the exponential scale height
  surveys quote: `sech^2(z/2h)` falls off as `exp(-|z|/h)` away from the
  plane. Only the thin disk carries the arms.
- **Thick disk** (`tuning.THICK_DISK_*`). Scale height 3 times and scale
  length 0.77 times the thin disk's, and 4% of the thin disk's density in
  the plane at the Sun, so `a_T = 0.04 * exp(R0/L_T - R0/L)`. It holds
  about 18% as many stars as the thin disk.
- **Bulge** (`tuning.BULGE_*`). The bar's boxy/peanut core as COBE/DIRBE
  sees it edge on: Dwek et al. 1995's G2 model. `bulge_amplitude` is its
  central density over the thin disk's. Its long axis leads the Sun by 27
  degrees in the direction of rotation (counterclockwise from galactic
  north), so its near end is at positive galactic longitude.
- **The Sun** sits at `R0 = 3.15 L` (8.2 / 2.6 kpc), on the inter-arm
  minimum. The bar's angle and the thick disk's amplitude are measured
  from there. `model_terms(shape)` returns these derived values; the
  Galaxy Map's prisms read them with the shape.

`sech^2` is computed as `4e / (1 + e)^2` with `e = exp(-2|x|)`
(`_sech_squared`), because `1 / cosh(x)^2` overflows once `|x|` passes about
710.

**Normalization.** `build_galaxy_shape` computes `k_norm` so
`relative_density = 1.0` in the plane, at the inter-arm minimum, at
`calibration_radius_pc`. That defaults to `R0` (8,200 pc at the default
shape). The model's fixed parts are ratios to the shape's own lengths,
so the same formulas work for a galaxy of any size.

**Halo floor.** The density never falls below
`tuning.MIN_RELATIVE_DENSITY` (GEN.78).

## 2. Default shape (`planetgen plan`)

| Parameter | Default | Basis |
|---|---|---|
| `--disk-scale-length-pc` | 2,600 | Thin disk, 2.6 +- 0.5 kpc (BHG16) |
| `--disk-scale-height-pc` | 300 | Thin disk, 300 +- 50 pc (BHG16) |
| `--bulge-scale-radius-pc` | 1,580 | Dwek et al. 1995 G2 fit to COBE/DIRBE |
| `--bulge-amplitude` | 3.11 | Bulge 31% of the stars: 1.4-1.7e10 of 5e10 M_sun (BHG16) |
| `--arm-count` | 2 | Grand-design spiral |
| `--pitch-angle-deg` | 15 | Typical grand-design pitch (10-25 degrees) |
| `--arm-amplitude` | 0.4 | Arm to inter-arm contrast of 2.33 |
| `--calibration-radius-pc` | 3.15 x scale length | The Sun's radius, 8.2 kpc |

The shape is stored in `galaxy_shape` with `k_norm`.

Sources: Bland-Hawthorn & Gerhard 2016, ARA&A 54:529 ("BHG16"); Dwek et
al. 1995, ApJ 445:716; Wegg & Gerhard 2013, MNRAS 435:1874.
`tests/test_galaxy_milky_way.py` and the math check
(`galaxy_bulge_to_total`, `galaxy_thick_to_thin_disk`) compare the
defaults with them.

## 3. Which sectors exist

```
E = expected_system_count_at_density_1   # 6.31 systems for a 4 pc sector
predicted_star_count(center) = E * relative_density(center)
a sector qualifies  <=>  predicted_star_count >= 1
```

The threshold is deterministic: a cell's fate is a fact of its position, not
a coin flip. At a given ring and layer the bulge is largest along the bar
and `arm_factor` is at most `1 + arm_amplitude`, so each term at its own
maximum bounds the density over angle (exactly, away from the bar), and
that bound falls as `R` and `|z|` grow. So `galaxySkeleton.build_layer_extents` finds
each layer's outer ring in one pass, and that outline is the galaxy's edge
(`galaxy-coordinate-system.md`, section 4). A sector inside the outline is
checked exactly when it is visited.

When a qualifying sector is generated, its own `relative_density` is the
`--density` multiplier, so its system count follows the local richness.

## 4. Stellar populations (7.38.0)

`population_densities` splits the density at a point into young,
intermediate, old disk and bulge populations. The bulge term is the bulge
population; the thick disk and the halo floor are old. The thin disk is shared by each population's share of star
formation, its own scale height (`STELLAR_POPULATION_SCALE_HEIGHT_RATIO`,
thinner for younger stars: 50, 150 and 300 pc exponential at the default
shape) and its own arm contrast
(`STELLAR_POPULATION_ARM_AMPLITUDE`, strongest for young stars). The parts
always sum to `relative_density`. Each system in a galaxy sector draws its
star's population from this mix, so O and B stars and supergiants gather in
the arms near the plane (CHANGELOG 7.38.0).

The bright-star scatter and backfill draw from these populations. How to
draw them without visiting every cell (an exact block-first draw from the
ring-and-layer density bound), why dropping sectors and boosting the rest is
rejected, and what the backfill and a fill of the sparse outer rim cost, are in
[sampling-backfill-and-resume.md](sampling-backfill-and-resume.md).

## 5. Why it works this way

- **Exponential disk, bulge and log spiral.** The brief was "assume a
  spiral galaxy similar to the Milky Way", and each parameter has a real
  Milky Way basis (section 2). Why this form was picked over other profiles
  is not recorded.
- **Calibrated between the arms.** The Sun sits between two arms (the
  Local Arm or Orion Spur), and calibrating there means typical
  solar-neighborhood space reads as 1.0. Calibrating on an arm would make
  everything off-arm look thinner than the real neighborhood.
- **`k_norm` computed, never typed in.** It must change whenever any shape
  parameter changes.
- **A Milky Way bulge and thick disk (GEN.118, GEN.119).** Revision 2
  shrank the bulge to an amplitude of 1 and a 200 pc sphere so the old
  shell build stayed tractable. That left it 0.6% of the stars against
  the Milky Way's ~31%, and edge on the galaxy was a flat disk. The
  bulge is now the COBE bar and the disk heights are the published ones,
  so a side view looks like the real galaxy and the model can be tested
  against published numbers. The cost: at the default shape the outline
  reaches layer 1,020 (4.1 kpc) instead of 317, and a full galaxy holds
  about 3.4 times as many systems (see the changelog for the numbers),
  close to the Milky Way's real 100-400 billion stars.
- **A hard threshold instead of a random occupancy roll.** Revision 1 used a
  hash-based Bernoulli roll; revision 2 replaced it with "a sector exists
  where it expects at least one star", which needs no randomness and
  defines the edge from the model.
- **Nothing stored per sector.** Revision 2 measured a full plan table at
  about 10.5 billion rows and about 1 TB. Since position and density are
  pure functions, only the outline is stored: one `galaxy_layer` row per
  layer (635) and one `galaxy_column` row per ring (about 3,900).
- **Populations by position.** Before 7.38.0 every sector drew stars from
  the same mix; tying ages to place puts young massive stars in the arms,
  as in a real disk. The reason in the release note is that "star ages
  follow where a sector sits".

### Alternatives considered and rejected

- The `galaxy_sector_plan` table and `galaxyPlan.py` resumable build
  (revision 2, sections 4-8 of the archive copy).
- A probabilistic occupancy roll per sector (revision 1, described only in
  revision 2; revision 1 itself is not in the repository).
- Metallicity gradients: listed as an open question in revision 2 and
  never built. (Its other open question, thin and thick disks as separate
  components, is built: GEN.119.)
- A spherical or oblate bulge (McMillan 2017's axisymmetric model was
  the other candidate for GEN.118): the observed bulge is a bar, boxy
  edge on, so the COBE fit was used. Not modelled: the long thin bar
  beyond the bulge, the nuclear stellar disk inside 200 pc and the
  stellar halo beyond the floor. With no central cusp the density is
  nearly flat through the core (a core sector holds about 880 to 1,060
  systems, far fewer than the real nuclear disk); see
  [fill-order-curves-and-core.md](fill-order-curves-and-core.md), section 4.

## Corrections made to the previous version

- It described the shell grid (`shell_index`, `shell_slot_index`,
  Fibonacci slots), which 6.0.0 replaced with rings, layers and slots.
- It used `E = 4.32` for an 11.5 ly sector; the 4 pc sector gives 6.31.
- It calibrated at `R_sun` from `GALACTIC_CENTER_DISTANCE_LY`; the code
  calibrates at a fixed ratio of the scale length (3.15 since GEN.118).
- It pointed to `galaxyGen.py`, `galaxyPlan.py`, `sectorGen.py` and SQLite,
  none of which are used now.
