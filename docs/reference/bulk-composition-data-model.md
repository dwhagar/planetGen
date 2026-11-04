# The bulk composition data model

What a planet or moon is made of in bulk, which stored values record it, and how its mass, radius, density and gravity follow
from it. The code is the authority; each module named here opens with a docstring that gives its sources and limits, and the
designs are [bulk-composition.md](../design/composition/bulk-composition.md) and
[composition-spine-contracts.md](../design/composition/composition-spine-contracts.md). The air is in
[atmosphere-data-model.md](atmosphere-data-model.md); column types and versions are in [database-schema.md](database-schema.md).

**The answer to "does the composition drive the mass, density and gravity?"** Yes, for every body the generator makes. Before the
composition build (main at 8.0.1635) a rocky body's density was a uniform draw for its class, its composition one fixed text per
class, and no rock, ice, metal or gas fractions existed. Now the composition is drawn first, from the host star's abundances and the
body's origin, and the radius, density and gravity follow from it. No class letter is read by the draw.

## 1. The order of the work

`physics/planets.generate_planet_properties`:

1. **Kind and mass.** A planet draws its kind (terrestrial or giant) and its mass from the formation law at its distance
   (`planet_kind.draw_formation`); a moon's mass comes from its parent's moon size law; a caller may give a mass, a radius or both
   (kept, and checked against `tuning.GIVEN_MASS_KG_BOUNDS` and `GIVEN_RADIUS_KM_BOUNDS`).
2. **A terrestrial body** (one the formation law does not make a giant, up to `bulk_draw.MAX_MASS_EARTH` = 20 Earth masses) draws its composition and takes its radius from it
   (`bulk_draw.draw_terrestrial`).
3. **A giant** draws its heavy elements, radius and interior layers (`giant_draw.draw_giant`).
4. The volume, the surface gravity (`planets.calculate_surface_gravity`, `G M / R^2` in Earth gravities), the radiogenic content and
   the interior heat follow. The air, the water and the ground come after (see the atmosphere page and
   `physics/surface_composition.py`).

Mass, radius and gravity are never drawn apart: the mass and the composition say the radius (or a given radius says the mass),
the density is the mass over the volume of the radius, and the gravity is `G M / R^2`.

## 2. What is stored

On `planets` and `moons` alike (schema v144, DB.69; every column is `NULL` on a row saved before it, which `bulk_fractions` reads
back as `None`):

| Column | Unit | What it is |
|---|---|---|
| `metal_mass_fraction` | share of the body's mass | The iron-rich core, light elements included |
| `rock_mass_fraction` | share | The silicate mantle and crust, or silicon carbide and graphite when `carbon_rock` |
| `water_mass_fraction` | share | H2O in every phase |
| `ice_mass_fraction` | share | Every other ice (CO2, NH3, CH4, CO, N2, H2S, CH3OH, organics) as one |
| `envelope_mass_fraction` | share | Hydrogen and helium gas |
| `rock_feo_fraction`, `rock_al2o3_fraction` | share of the rock | FeO and Al2O3 |
| `rock_mg_si_ratio` | molar ratio | The rock's Mg/Si |
| `envelope_helium_fraction`, `envelope_z_fraction` | share of the envelope | Helium (Y) and heavy elements (Z) |
| `porosity` | share of volume | The void share of a small body (`NULL` where the structure solver decides) |
| `carbon_rock` | flag | The rock is silicon carbide and graphite |
| `composition_origin` | `drawn` or `inverse` | Composition first, or worked back from a given mass and radius |
| `composition` | text | The phrase written from the fractions (section 6) |
| `mass_kg`, `radius_km`, `density_g_cm3`, `gravity_g` | kg, km, g/cm3, Earth gravities | The results; `density = mass / (4/3 pi radius^3)` |

The five shares of the whole body (metal, rock, water, ice, envelope) are not negative and add up to 1. A giant's interior adds
`interior_kind` (`gas_giant`, `ice_giant`, `brown_dwarf`), `core_mass_earth`, `dilute_core_mass_earth`, `metallization_pressure_gpa`,
`helium_depletion` and `has_helium_shell` (v154, `physics/giant_interior.py`); the layer masses add up to the body's mass and the layers
hold exactly the heavy-element mass and natal helium (tested to 1e-9). Where a body formed relative to its star's snow line is
`born_beyond_snow_line` (v132), and the star's own snow line `stars.snow_line_km`. The ground (v159) is a separate set of
`ground_*` columns, drawn from the bulk fractions and described in `physics/surface_composition.py`.

## 3. How a terrestrial body's composition is drawn

`physics/bulk_draw.py` (GEN.428), in this order (the order is part of the generator epoch):

1. **The star sets the dry ratios.** The host star's abundances (`stellar_composition`) give Fe/Mg, Si/Mg, Al/Mg, Ca/Mg, Ni/Mg and
   C/O. A rocky body takes all the star's Mg, Si, Ca and Al into its rock as oxides and the iron its mantle does not hold as FeO into
   its metal core (`dry_partition`). An Earth-like star gives a core mass fraction near 0.3 with no draw; an iron-rich or iron-poor
   star moves it, the observed trend [R1146, R3916].
2. **What the star cannot say is drawn**, from `RANGES` and the tables, each with its source in the module: the kind of rock
   (`ROCK_KIND_TABLE`: silicate, or silicon carbide and graphite for C/O of 0.8 or more, [R0455, R4469]); the
   water share, log-uniform over the range of the body's origin (`tuning.DRY_WATER_MASS_FRACTION` 1e-6 to 1e-3; a water world,
   `WATER_WORLD_SHARE` 0.03 of bodies, 0.01 to 0.1; an icy body beyond the snow line `ICY_WATER_MASS_FRACTION` 0.05 to 0.5); for an
   icy body the other ices from the comet volatiles; the rock's FeO (2 to 18 percent), the core's light elements (3 to 12 percent), a
   small icy body's porosity (0 to 0.3), a carbon body's graphite over rock (5 to 40 percent); then the iron branch of a rocky planet
   (`IRON_BRANCH_TABLE`: its iron follows the star's with a scatter, or is enriched to a core mass fraction of 0.5 to 0.8).
3. **Radius from mass.** `bulk_structure.structure` gives the radius the composition holds the mass at. A caller that gave a radius
   gets the mass the composition needs for it, by bisection; a body too dense to reach its radius below 20 Earth masses is set to
   the radius of 20 (`MAX_MASS_EARTH`). A body below `SMALL_BODY_MASS_KG` (1e21 kg) is a heap with a porosity and its mass is its
   volume times the mixing density.
4. **A body given both a mass and a radius** goes the other way (`composition_origin = inverse`): they say the core mass fraction
   (`bulk_structure.cmf_from_mass_radius`, at most `MAX_CMF` 0.9), as an observer's route does.
5. **A moon that holds its origin** (`moon_composition`, GEN.466) takes its shares from it: a disk moon of a giant (ice rising with the
   formation distance past the disk's ice line, [R3922]; the moon system of [R0643]), an impact moon (little
   metal, no ice, [R2724]), a captured body, a ring moon or the partner of a dwarf planet. Every weight of the origin table is
   open (no source), so the rows a moon supports share the probability equally. Below the taper's end of the compaction law
   (`compaction.TAPER_END_KG`, GEN.442) its mass follows from its radius by the tapered brittle law and its grain density; heavier, the
   structure solver decides.

Small bodies (belts, comets, icy bodies, rings) draw their composition first from rows of the materials table and derive their
density from it with `physics/compaction.py` (`small_body_composition.py`, GEN.442: the Bus and DeMeo taxa of the Solar System's one
belt are flagged stand-ins).

## 4. How the radius follows: the structure model

`physics/bulk_structure.py` (GEN.429), pure functions with no draws. The radius is at the 1 bar level (Boss, 1:24 AM EDT
2026-10-17). One reference solver and three fast forms:

- **`layered_structure`**, the reference: a two-layer integrator extended with a water layer on top whose density is the 300 K
  isotherm of liquid water (IAPWS-95 through the NIST WebBook, to 950 MPa), the ice VI, VII and X reference points of the materials
  table, and the polytrope of [R3992] above 1 TPa. Hydrostatic balance in `s = m^(1/3)`, fourth-order Runge-Kutta,
  shooting on the central pressure to 1 bar at the surface. Checked against a scipy reference built from the full mantle tables
  (0.20 percent in radius, 0.28 in the moment of inertia, 0.93 in the core-mantle pressure over 90 cases) and, with water, within 1.5
  percent of [R4881] from 0.5 to 50 Earth masses. Validated 0.1 to 20 Earth masses.
- **Dry rock and iron**, `radius_fit_earth`: 0.1 to 20 Earth masses, core mass fraction 0 to 0.8; at most 0.40 percent (rms 0.10)
  against the solver.
- **Water-bearing**, `water_world_radius_earth`: 0.3 to 20 Earth masses, water fraction 0 to 1; rms 1.2 and at most 3.6 percent.
- **Small differentiated bodies**, `small_body_radius_earth`: 1.5e-4 to 0.3 Earth masses; rms 0.57, at most 2.0 percent.
- **An envelope on a core**, `envelope_radius_earth`: [R2741] Eqs. 1 and 3 at the 20 mbar level plus ln(50) scale
  heights for 1 bar; 1 to 20 Earth masses, envelope 0.01 to 20 percent; rms 6.3 percent against their tables at 1 Gyr. The valley of
  the radius distribution (`physics/radius_valley.py`, GEN.509) tunes the starting envelope against the Kepler population.

A body whose envelope is more than `GIANT_ENVELOPE_FRACTION` (0.5) of its mass, or that has over 20 Earth masses of core, is a giant
and outside this module.

**A giant** (`giant_radius.py`, `giant_draw.py`, GEN.460): the envelope's metallicity is a multiple of the star's, the heavy-element
mass follows [R0710] (`M_Z = M_core + f_Z (M - M_core)` with its scatter, at most `MAX_HEAVY_FRACTION` 0.9 of the mass), and the
radius is interpolated from the cool-start tables of [R1429] (17 to 3,594 Earth masses) plus the inflation of
[R4029] above a threshold flux (an empirical median, a stand-in for the thermal evolution the project does
not have). Where the tables do not hold the body (an ice giant, a giant under 17 Earth masses, a moon) the anchored mass-radius
relation of `planets.sample_giant_radius_km` is used. The density of a giant is bounded by `tuning.GIANT_DENSITY_BOUNDS_G_CM3` (0.05 to
60) and a terrestrial body's by `bulk_draw.DENSITY_BOUNDS_G_CM3` (0.6 to 25).

## 5. The materials table

`physics/materials.py` (GEN.427): the one table of every element, mineral, ice, organic and alloy a body can contain, 161 rows from
`physics/data/materials.json`, which `scripts/build_materials_table.py` builds from the research parts' tables (E materials, B's
phases and iron rows, C's ice phases) and `scripts/materials_data.py` ([R3992: Tables 1 to 3], the Handbook of Mineralogy,
PubChem). A row has its formula, molar mass (kg/mol, from the element table), grain density (g/cm3 at ambient conditions), phase
points, latent heats, band albedos and, where it has one, the equation-of-state constants of the model GEN.426 chose, each number with
its source and a source state (`verified`, `partly verified` with its limit, or not reached). A value with no source is `None`; none is
kept from memory. The `gate` says which classes can hold a row: `any` for every planet and moon class and small-body kind (135 rows),
`planetary` for the planet and moon classes only (26 rows: the high-pressure silicate phases, high-pressure ices, core iron alloys and
hydrogen-helium envelopes). The loader refuses a file without its sources, a row without a source and state, a formula it cannot
read, an unknown gate, and a missing equation-of-state constant. The 84 asteroid components, the taxa and the class composition
texts all resolve to rows of it, and a test fails when one does not. The licence of the data is in `DATA-LICENSE.md`.

## 6. The words

`physics/composition_text.py` (GEN.430) writes the `composition` phrase of every planet and moon from the stored fractions, so the
words and the numbers cannot disagree: each material the body holds at least 1 percent of is named with the level its share falls in
(`mostly` from 0.5, `a large share of` from 0.2, `a moderate share of` from 0.05, `a small share of` from 0.01), largest first, for
example "mostly silicate rock, with a large share of iron-rich metal and a small share of water". A porosity of 0.1 or more adds
"loosely packed".

## 7. The migration and the sweep

Migration 0144 adds the columns, all `NULL`: a row made before it has no composition to read back (it needs its star's abundances and
the draws), so the regenerated galaxy of the generator epoch replaces it. A stopped run can be started again.

`tests/test_sweep_bulk.py` (TEST.210) runs the rules of `tests/audit/sweep_bulk.py` over bodies of every class: `fractions-sum`,
`detail-ranges`, `porosity-range`, `density-from-mass-radius`, `gravity-from-mass-radius`, `radius-from-structure` (to 1e-6 for a
planet and 4 percent for a moon, whose mass follows the tapered brittle law; a body with an envelope is left out because its flux and
age are not stored), `denser-than-grains` and `text-materials`. Its last test breaks each rule and checks the sweep reports it.

## 8. Limits

The Fortney tables do not inflate, so the inflation is the Sestovic median (a stand-in, marked in the gap plan); the water-world share
and the origin weights of moons are open; the small-body taxa kernel is the Solar System's one belt and is flagged. These stand-ins are
registered with their tier and tweak level in [stand-in-values.md](../design/research/stand-in-values.md).
