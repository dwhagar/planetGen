# Galaxy density model

**Status:** built (`src/planetgen/galaxy/density.py`). The model dates
from the "revision 2" design pass of 2026-09 (Track C). That revision's main
proposal, a stored plan of every qualifying sector, was never built; it is
kept in `archive/galaxy-disk-density-rev2.md`. This file describes the model
as the code uses it now, on the cylindrical grid
(`galaxy-coordinate-system.md`).

## 1. The model

An exponential disk, an exponential bulge and a logarithmic spiral:

```
r_cyl = sqrt(x^2 + y^2)
r_3d  = sqrt(x^2 + y^2 + z^2)
theta = atan2(y, x)

rho_bulge(r_3d)          = bulge_amplitude * exp(-r_3d / bulge_scale_radius_pc)
rho_disk_radial(r_cyl)   = exp(-r_cyl / disk_scale_length_pc)
f_z(z)                   = sech^2(z / disk_scale_height_pc)
theta_arm(r_cyl)         = spiral_reference_angle_rad
                           + ln(r_cyl / spiral_reference_radius_pc) / tan(pitch_angle_rad)
arm_factor(r_cyl, theta) = 1 + arm_amplitude * cos(arm_count * (theta - theta_arm))

relative_density = k_norm * (rho_bulge + rho_disk_radial * f_z * arm_factor)
```

`sech^2` is computed as `4e / (1 + e)^2` with `e = exp(-2|x|)`
(`_sech_squared`), because `1 / cosh(x)^2` overflows once `|x|` passes about
710.

**Normalization.** `build_galaxy_shape` computes `k_norm` so
`relative_density = 1.0` in the plane, at the inter-arm minimum, at
`calibration_radius_pc`. That defaults to `2.82 * disk_scale_length_pc`
(the Milky Way's solar-radius to scale-length ratio), about 7,900 pc at the
default shape. The module does not read `GALACTIC_CENTER_DISTANCE_LY`, so
the same formulas work for a galaxy of any size.

## 2. Default shape (`generate.py plan`)

| Parameter | Default | Basis |
|---|---|---|
| `--disk-scale-length-pc` | 2,800 | Real Milky Way scale |
| `--disk-scale-height-pc` | 350 | Between the thin disk (~300 pc) and the thick disk (~900-1,000 pc) |
| `--bulge-scale-radius-pc` | 200 | Shrunk from a first guess of 500 (see section 5) |
| `--bulge-amplitude` | 1.0 | Shrunk from a first guess of 5 (see section 5) |
| `--arm-count` | 2 | Grand-design spiral |
| `--pitch-angle-deg` | 15 | Typical grand-design pitch (10-25 degrees) |
| `--arm-amplitude` | 0.4 | Arm to inter-arm contrast of 2.33 |
| `--calibration-radius-pc` | 2.82 x scale length | As above |

The shape is stored in `galaxy_shape` with `k_norm`.

## 3. Which sectors exist

```
E = expected_system_count_at_density_1   # 6.31 systems for a 4 pc sector
predicted_star_count(center) = E * relative_density(center)
a sector qualifies  <=>  predicted_star_count >= 1
```

The threshold is deterministic: a cell's fate is a fact of its position, not
a coin flip. Because `arm_factor` is at most `1 + arm_amplitude`, the
density's maximum over angle at a given ring and layer is exact, and it
falls as `R` and `|z|` grow. So `galaxySkeleton.build_layer_extents` finds
each layer's outer ring in one pass, and that outline is the galaxy's edge
(`galaxy-coordinate-system.md`, section 4). A sector inside the outline is
checked exactly when it is visited.

When a qualifying sector is generated, its own `relative_density` is the
`--density` multiplier, so its system count follows the local richness.

## 4. Stellar populations (7.38.0)

`population_densities` splits the density at a point into young,
intermediate, old disk and bulge populations. The bulge term is the bulge
population. The disk term is shared by each population's share of star
formation, its own scale height (`STELLAR_POPULATION_SCALE_HEIGHT_RATIO`,
thinner for younger stars) and its own arm contrast
(`STELLAR_POPULATION_ARM_AMPLITUDE`, strongest for young stars). The parts
always sum to `relative_density`. Each system in a galaxy sector draws its
star's population from this mix, so O and B stars and supergiants gather in
the arms near the plane (CHANGELOG 7.38.0).

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
- **Bulge shrunk to amplitude 1 and 200 pc.** With amplitude 5 and 500 pc
  the bulge alone cleared the one-star threshold out to about 2,950 pc in
  every direction, a solid ball that made the shell-era pruning useless.
  The smaller pair still gives a bright core. The design notes say these
  values made the build tractable and are not a fitted astrophysical
  result.
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
- Thin and thick disk as separate components, and metallicity gradients:
  listed as open questions in revision 2 and never built. The population
  split of section 4 covers part of the first.

## Corrections made to the previous version

- It described the shell grid (`shell_index`, `shell_slot_index`,
  Fibonacci slots), which 6.0.0 replaced with rings, layers and slots.
- It used `E = 4.32` for an 11.5 ly sector; the 4 pc sector gives 6.31.
- It calibrated at `R_sun` from `GALACTIC_CENTER_DISTANCE_LY`; the code
  calibrates at 2.82 scale lengths.
- It pointed to `galaxyGen.py`, `galaxyPlan.py`, `sectorGen.py` and SQLite,
  none of which are used now.
