# Galaxy density model

**Status:** built (`src/planetgen/galaxy/density.py`). The model dates
from the "revision 2" design pass of 2026-09 (Track C). That revision's main
proposal, a stored plan of every qualifying sector, was never built; it is
kept in `archive/galaxy-disk-density-rev2.md`. This file describes the model
as the code uses it, on the cylindrical grid until 2026-10-15 (the grid is shell, band and wedge since then,
`galaxy-coordinate-system.md`).

## 1. The model

A thin disk carrying the spiral arms, a thick disk and a boxy bar bulge,
each fitted to the Milky Way's published structure (GEN.118, GEN.119):

```
R     = sqrt(x^2 + y^2)
theta = atan2(y, x)
v(z, h) = exp(-|z| / h)              # the form surveys fit (GEN.544); column 2 h

thin(R, z)  = exp(-R / L) * v(z, h)                       # L, h: the shape's
                                                          # disk scale length/height
thick(R, z) = a_T * exp(-R / L_T) * v(z, h_T)             # L_T, h_T: the drawn length and height
                                                          # ratios times L and h (GEN.532; height
                                                          # ratio 2.5 to 4.2); 0.77 and 3 are the
                                                          # Milky Way reference check
bulge(x, y, z) = bulge_amplitude * exp(-r_s^2 / 2)        # in the bar's frame:
r_s^4 = ((x'/x0)^2 + (y'/y0)^2)^2 + (z/z0)^4              # x0 = bulge_scale_radius_pc,
                                                          # y0/x0, z0/x0 drawn per galaxy (GEN.533);
                                                          # 0.39 and 0.27 are the Dwek reference
theta_arm(R)      = spiral_reference_angle_rad + ln(R / spiral_reference_radius_pc) / tan(pitch_angle_rad)
arm_factor(R, th) = 1 + arm_amplitude * cos(arm_count * (th - theta_arm))

relative_density = k_norm * (bulge + thin * arm_factor + thick)
```

- **Thin disk.** `disk_scale_height_pc` is the exponential scale height
  surveys quote, and the profile is that exponential, `exp(-|z|/h)` (GEN.544): its midplane density is the one the
  fits give, and its column is `2 h` times it. It replaced the isothermal sheet `sech^2(z/2h)`, whose tail is the same
  but whose column is `4 h` times its midplane, so for one mass its midplane was half the fits'. Only the thin disk
  carries the arms. The nuclear stellar disc keeps the isothermal sheet (`_sech_vertical`).
- **Thick disk** (`shape.structure`, `galaxy/structure_draw.py`, GEN.532). Its mass, scale length and scale height are
  ratios to the thin disk's drawn from the galaxy seed: the mass (0.14 to 0.99 at the circular speeds of 3e9 to 2e11
  solar masses of stars) between the fits of Yoachim and Dalcanton 2006 and Comeron et al. 2018, the length by Comeron
  2018 and the height from 2.5 to 4.2. The amplitude is derived from the mass ratio:
  `a_T = ratio * arm_level / (length_ratio^2 * height_ratio)`. The Milky Way's published values (900 pc, 2.0 kpc, 4% of
  the thin disk's density at the Sun, 18% of its mass; BHG16) are `structure.REFERENCE_STRUCTURE`, a check target.
- **Bulge** (`tuning.BULGE_*`). The bar's boxy/peanut core as COBE/DIRBE
  sees it edge on: Dwek et al. 1995's G2 model. `bulge_amplitude` is its
  central density over the thin disk's. Its long axis angle is drawn per
  galaxy, uniform on [0, 180) degrees, with the axis ratios (GEN.533,
  `galaxy/structure_draw.py` `draw_bar`); the Milky Way's 27 degrees leading
  the Sun, and its 0.39 and 0.27, are the Dwek reference, a check only.
- **The Sun** is not placed: no Sun position is read and no point of the
  galaxy is the calibration point (GEN.404; Boss, 2026-10-16, 9:43 PM EDT:
  the galaxy stores no home or Sun position). `R0 = 3.15 L` (8.2 / 2.6 kpc)
  is a test reference only. `model_terms(shape)` returns the derived values; the
  Galaxy Map's prisms read them with the shape.

The nuclear disc's `sech^2` is computed as `4e / (1 + e)^2` with `e = exp(-2|x|)`
(`_sech_squared`), because `1 / cosh(x)^2` overflows once `|x|` passes about
710. The disks' `exp(-|z|/h)` is clamped at 700 scale heights.

**Normalization (GEN.404).** `build_galaxy_shape` computes `k_norm` so the density, in units of the real local
stellar density (`LOCAL_STELLAR_DENSITY_LY3`), integrates to the systems the galaxy's stellar mass makes
(`structure.stellar_mass_solar` over `density.mean_system_mass_solar`, GEN.545): `density.total_mass` is the integral of the
unnormalized model (the thin and thick disks, the bulge, the nuclear cluster and disc and the stellar halo, each in
closed form). No point of the galaxy is the calibration point and no Sun position is read. A galaxy of the reference
structure (5e10 solar masses) held 1.25e11 systems at the recalled 0.4 solar masses a system (1.03e11 at the derived mean, below), and the density at the reference radius was 0.38 of the local
density the earlier calibration pinned there (2.6 times sparser): the model's disc was a `sech^2` sheet whose midplane
density is half the exponential fits quote for the same column. **GEN.544** made the vertical profile the exponential
the scale heights are fitted with, with the column (the mass) unchanged: the density at the reference radius is now 0.75
of the old pinned value (twice the GEN.404 value, 0.376 to 0.750), the thin disc's column is the same, and the
reference galaxy's bulge amplitude is 3.11/2 because BHG16's 3.11 was measured against the `sech^2` disc and the bulge
keeps its 31% of the mass (`mathcheck.REFERENCE_BULGE_AMPLITUDE`). The model's fixed parts are ratios to the shape's own lengths,
so the same formulas work for a galaxy of any size.

**The mean system mass (GEN.545).** The galaxy's number of systems is its stellar mass over the mean mass of a system, and
that mean is no longer a recalled constant (`tuning.MEAN_SYSTEM_MASS_SOLAR` was 0.4, "within a factor of 1.5"); it is worked
out from the generator's own draws (`generation/system_mass.py`) for stars formed at one age:

- the **placed star**: an initial mass from the initial mass function (Kroupa, 0.08 to 150 solar masses), at its birth mass
  until its last giant phase ends (`stellar_evolution.youngest_progenitor_sol`), then its remnant: a white dwarf by the
  initial-final relation, a neutron star or a black hole at the mean of the kind's mass range. The remnants are systems too;
- the **other stars** of its system: the `star-count` table's mean number by the placed star's class with the 8 to 12 tail
  (`multistar.expected_star_count`: 1.31 stars for a 0.1 solar mass primary, 1.62 for 1, 2.37 for 5, 3.08 for 50), each
  with a birth mass of the primary's times a uniform 0.1 to 1, taken to its own state at the same age.

The mean over a population's ages weights each by the star formation there (the disk's three populations by their history,
`galaxy.age_law`; the others uniformly over their windows), and the mean over a galaxy or a cell weights each population by
its systems (`density.population_units`: the thin disk shared among its populations, the thick disk in `old`, the bulge, the
nuclear cluster and disc, the halo). By systems and not by mass because the density counts systems. For the reference
galaxy: young 0.72, intermediate 0.57, old 0.51, bulge 0.50, halo 0.49, so **0.508 solar masses** a system overall, a
primary drawn from the primary mass function (GEN.928: the Chabrier function less the companions of heavier primaries; the
Kroupa function gave 0.483 with these words: a placed star of 0.35 times 1.42 stars a system) times 1.43 stars a system. The
reference galaxy holds 0.98e11 systems (1.03e11 on the Kroupa function, 1.25e11 at 0.4) and its density at the reference
radius is 0.59 of the pinned local density (0.62 on the Kroupa function, 0.750 at 0.4): the gap
that remains is real if the pinned density (Mamajek's 0.0984 stars per pc^3) counts the stars of a system singly, as it
probably does, since a system here holds 1.42 of them. The compact companions of massive stars' trees add under 0.2% (not
counted), and so do the companions of the remnants the scatter places (a few percent of those systems). A test draws the
systems with the generator's own samplers and compares (`test_gen_system_mass.py`). A cluster's or satellite's mean (once the recalled
`CLUSTER_MEAN_SYSTEM_MASS_SOL`, also 0.4) is derived the same way from the component's own age (GEN.598,
`Component.mean_system_mass_solar`): 0.47 for an 11 Gy globular cluster, 0.49 for the Large Cloud, 0.9 for a cluster a few
thousand years old.

**Stellar halo** (GEN.234). A spheroid with a broken power law, `A * (max(m, core) / break) ** -n`, where
`m = sqrt(R^2 + (z / q)^2)`: inner slope 2.5, outer slope 4.5, break at 25 kpc (9.6 disk scale lengths), axis ratio
`q = 0.7`, flat inside 500 pc (`tuning.HALO_*`; sources and the numbers it gives are in
[intergalactic-and-halo-contents.md](intergalactic-and-halo-contents.md) sections 2.1 and 4). The radii are ratios to the
disk scale length, and the amplitude is solved like `k_norm`: the halo holds the galaxy's own share (`structure.stellar_halo_system_fraction`, GEN.532: its drawn
halo mass over its stellar mass, scaled to the model's counts; 1.17% for the Milky Way's 1.4e9 solar masses of 5e10 at 0.4 solar
masses a system against 3.0e11 systems) of the visible components'
stars (`density.halo_terms`). It replaces the flat floor GEN.78 put under the density (0.1% of the local density
everywhere, which was half the real halo at the Sun and 200 to 340,000 times too much at 50 to 262 kpc). The halo is
0.19% of the local density at the Sun and 3e-9 at 262 kpc, and no cell beyond 8.7 kpc from the center expects more than
about 0.01 systems from it. The alternative normalisation to the local density (Juric et al. 2008, 0.5%, about 2.7 times
higher) is `tuning.HALO_LOCAL_DENSITY_FRACTION`; `density.halo_local_density_scale` gives the factor. The outline
(`galaxy/skeleton.py`) reads the disk, the bulge and the components (the clusters and satellites of `galaxy_components`,
added to the bound by `ComponentField.bound`, GEN.160 and GEN.237) but not the halo, so the halo adds no cell to it. The
halo still fills every cell of the reach: a cell outside the outline is generated from its own density, halo term
included, and its sector profile says `halo` (`galaxy/profile.py`, GEN.235; see
[intergalactic-and-halo-contents.md](intergalactic-and-halo-contents.md) section 5.3). The halo is the `halo` entry of the
profile's `mix`, with the age window of the halo population (`tuning.STELLAR_POPULATION_AGE_RANGES_GY["halo"]`) and its [Fe/H]
(`metallicity.halo_mean_fe_h`, `tuning.METALLICITY_HALO`).

## 2. The shape (`planetgen plan`)

GEN.401: no shape parameter has a default value any more. `planetgen plan` draws each from the galaxy seed
(`galaxy/shape_draw.py`, the laws in `intergalactic-and-halo-contents.md` 5.13, rows 4 to 6) unless the option is
given; an option given replaces its draw, and the bulge's amplitude is solved for the disk actually used.

| Parameter | Where it comes from |
|---|---|
| `--disk-scale-length-pc` | Size-mass relation of late-type galaxies (Shen et al. 2003), `R_d = R_50 / 1.678`, log-normal scatter |
| `--disk-scale-height-pc` | Length over `q = 7.3` with a 0.13 dex scatter (Kregel et al. 2002) |
| `--bulge-scale-radius-pc` | A third of the bar's length, which is drawn in disk scale lengths by stage (Erwin 2005) |
| `--bulge-amplitude` | Solved so the bulge's mass is the drawn bulge-to-total of the bulge and disks (Graham and Worley 2008, by stage) |
| `--arm-count` | Shares of Hart et al. 2017 (two arms in 64.5 percent) |
| `--pitch-angle-deg` | Falls with stellar mass between Yu and Ho 2019's two binned means (earlier notes say 2018; Crossref gives ApJ 871, 194, 2019 for the pitch-angle paper), 5.5 degrees of scatter |
| `--arm-density`, `--interarm-density` | 1.4 and 0.6, the arm contrast (not drawn; no source) |

The Milky Way's measured shape (BHG16: 2,600 pc, 300 pc, 1,580 pc, 3.11, two arms at 15 degrees) is kept only as the
reference the math check measures the model's machinery against (`physics/mathcheck.py`).

The shape is stored in `galaxy_shape` with `k_norm`.

Sources: Bland-Hawthorn & Gerhard 2016, ARA&A 54:529 ("BHG16"); Dwek et
al. 1995, ApJ 445:716; Wegg & Gerhard 2013, MNRAS 435:1874.
`tests/test_galaxy_milky_way.py` and the math check
(`galaxy_bulge_to_total`, `galaxy_thick_to_thin_disk`) compare the
model with them on the reference shape.

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
intermediate, old disk, bulge, halo and cluster populations. The bulge term is the bulge
population, the halo term the halo population (10 to 13 Gyr, [Fe/H] about -1.6 inside 25 kpc and -2.2 outside) and the
thick disk is old. The cluster population (10 to 13 Gyr) is 0 until the cluster components (GEN.237) add their
density. The thin disk is shared by each population's share of star
formation, its own scale height (`structure.population_height_ratios`, GEN.532: the squared vertical
dispersion of the population over the old disc's, thinner for younger stars) and its own arm contrast
(`STELLAR_POPULATION_ARM_AMPLITUDE`, strongest for young stars: 0.75, 0.42 and 0.10 for the young, intermediate and old,
GEN.931; the old disk is the observed contrast of about 10% (Eilers et al. 2020; Drimmel and Spergel 2001), the young one the
gas ratio of 4 to the Kennicutt-Schmidt index 1.4, and the intermediate one falls with the log of the age between them;
registry 6.65). The parts
always sum to `relative_density`. Each system in a galaxy sector draws its
star's population from this mix, so O and B stars and supergiants gather in
the arms near the plane (CHANGELOG 7.38.0).

The bright-star scatter and backfill draw from these populations. How to
draw them without visiting every cell (an exact block-first draw from the
ring-and-layer density bound), why dropping sectors and boosting the rest is
rejected, and what the backfill and a fill of the sparse outer rim cost, are in
[sampling-backfill-and-resume.md](../generation/sampling-backfill-and-resume.md).

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
  beyond the bulge and the stellar halo's inner structure. The nuclear star
  cluster and nuclear stellar disc are modelled (GEN.280, schema v99,
  [realistic-galactic-core.md](realistic-galactic-core.md) section 3); see
  [fill-order-curves-and-core.md](fill-order-curves-and-core.md), section 4.

## Corrections made to the previous version

- It described the shell grid (`shell_index`, `shell_slot_index`,
  Fibonacci slots), which 6.0.0 replaced with rings, layers and slots.
- It used `E = 4.32` for an 11.5 ly sector; the 4 pc sector gives 6.31.
- It calibrated at `R_sun` from the Sun's distance; the code calibrated at a fixed ratio of the scale length (3.15
  since GEN.118) until GEN.404 normalized the whole galaxy to its stellar mass.
- It pointed to `galaxyGen.py`, `galaxyPlan.py`, `sectorGen.py` and SQLite,
  none of which are used now.

## ADM.49: user-set arm, inter-arm, core and bulge densities

The thin disk's arm factor is `arm_level * (1 + arm_amplitude * cos(m (theta - theta_arm)))`. The user sets the density on an arm crest `a` and between arms `i` (the disk's own units, usual 1.4 and 0.6), and `density.arm_terms(a, i)` gives `arm_amplitude = (a - i) / (a + i)` and `arm_level = (a + i) / 2`. The Gaussian core this item added was replaced by the nuclear star cluster and nuclear stellar disc (GEN.280, schema v99, [realistic-galactic-core.md](realistic-galactic-core.md) section 3); `core_amplitude` became `nuclear_amplitude`. The bulge density is `bulge_amplitude`. Superseded (GEN.404; Boss, 2026-10-16, 9:43 PM EDT): no point is the calibration point and no Sun position is read; the total is set by the galaxy's stellar mass, so a denser arm or core raises those places against the rest of the galaxy (8.2 kpc is a test reference only). The skeleton and object-first bounds, the stellar populations, the component masses and the Galaxy Map's prisms all read the same terms.

## GEN.625: the all-Auto window of the size and the flatness

A galaxy with the mass, size (or compactness), flatness, bulge share and bulge size all on Auto draws its disc size and flatness through `shape_draw._auto_size`: the seed's draws of the size-mass relation (Shen et al. 2003, median half-light radius over 1.678 for the scale length, scatter 0.47 falling to 0.34 in ln) and of the flatness (Kregel et al. 2002, median 7.3, 0.13 dex), and then again from the same two streams, one more draw each, until the shape raises none of `dials.shape_warnings` (the thin disc's crest at least 20 times the one-star density, a scale height of at least 8 pc, a bulge of at least 4 pc) and the flatness is inside the Structure dial's 3 to 20. The crest tested is the plan's own (`plan_figures.crest`, GEN.626) of the shape built at every dial's Auto value and the default sector edge, which costs one density normalization a draw (the disc history is cached, so the accepted draw's build is free). At most 32 draws are made; the shape with the best crest is kept if none passes (it does not happen at the Auto masses). The first draw is the unchanged one, so a seed that did not warn keeps its galaxy; one seed in five drew a shape that warned and gets a later draw. Pin any of the five and the draws are the seed's own first ones, as before, and the grid's warnings fire for what the pin makes.

Why a window and not a shifted law: the crest goes as the mass over the scale length cubed times the flatness, and the size-mass scatter is a log-normal with no upper bound, so a bounded sourced law cannot keep every draw clear of the warning; the Auto mass is already a window (3e9 to 6e9, Boss, 10-21) and this is the same kind of statement, the galaxies the 4 pc plan represents well. The sourced medians and scatters are untouched; the window removes the large, thick tail the plan cannot hold.

Measured over 1,000 seeds before and after, against the plan's own crest (with the nuclear disc of GEN.628): the share of Auto draws that warn 17.5 percent to 0, the mean disc scale length 1,470 pc to 1,224 pc (47 percent of the 2.6 kpc a Milky Way disc has, Bland-Hawthorn and Gerhard 2016, used only as a test yardstick), the mean stellar mass 4.3e9 solar masses (8 percent of a Milky Way's 5.4e10), the mean flatness 7.6 to 7.8, the 5th percentile of the crest 6.1 to 23.6.