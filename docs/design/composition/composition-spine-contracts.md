# Composition spine: interface contracts (lane I, 2026-10-23)

Owner: lane I (thread "Composition spine for habitability"). Items: GEN.565 (PlanetaryObject part), GEN.557, GEN.428 with DB.69, GEN.429, GEN.452 part 2, GEN.411 with DB.71 and GEN.424. Everything after GEN.428 (GEN.430, 441, 442, 460, 461, 466 and the rest) belongs to the lanes named in the plan, not to this lane. Designs: docs/design/composition/composition-model.md, bulk-composition.md, atmosphere-composition.md, atmosphere-colour.md section 7.

Status: the accessors of section 1 are in PR 1 (no behaviour change, no epoch); sections 2.1 and 2.2 describe the code as built in the GEN.429 and GEN.428 PR.

## 0. Order of work (and one change to the order the plan gives)

The plan lists GEN.565 then GEN.557 then GEN.428. GEN.428, GEN.429 and GEN.411 do not need the PlanetaryObject class or the wanderer draw: they are functions of a plain `BulkSubject` / `AirSubject` (below), called from `generate_planet_properties` today and from the one body draw later. GEN.557 also waits on GEN.553 (lane A) and DB.73 (lane B). So the order is:

1. PR "spine-contracts": this document in docs/design, plus `physics/bulk_composition.py` (data class + accessor) and `physics/mixture.py` (accessors). No epoch.
2. PR GEN.429 kernels: `physics/bulk_structure.py`, pure functions (reference solver, closed forms, small-body mixing). No epoch.
3. PR GEN.428 + DB.69 (+ the derived radius, density and gravity of GEN.429 in the generator): epoch bump, migration, golden update.
4. PR GEN.452 part 2 (#2613 rebased) + GEN.411 + DB.71 + GEN.424: epoch bump, migration.
   Built in three merges (2026-10-24): 4a GEN.452 part 2 (#3041), 4b DB.71 + GEN.424 (#3110, epoch 95), 4c GEN.411 + GEN.422 (epoch 103): the air is the fixed point of `physics/surface_air.py` (`surface-pressure-supply-law.md` section 6), which replaced the class density, molar density and greenhouse ranges, `CLASS_MIXES` and the mixing-ratio draw (GEN.683); the supply law (GEN.726, #3111) came between 4b and 4c.
5. GEN.565 PlanetaryObject and GEN.557 once lanes A and B land.

## 1. Accessors (stable from PR 1)

### 1.1 Bulk composition: `planetgen.physics.bulk_composition`

```python
@dataclass(frozen=True)
class BulkFractions:
    metal: float      # iron-rich core incl. light elements, share of the WHOLE body's mass
    rock: float       # silicate mantle and crust (SiC and graphite when carbon_rock)
    water: float      # H2O in every phase (ice, ocean, mantle-bound is NOT included)
    ice: float        # every other ice: CO2, NH3, CH4, CO, N2, H2S, CH3OH, organics (one share)
    envelope: float   # H/He gas (with envelope_z of heavy elements inside it)
    # metal + rock + water + ice + envelope == 1 to 1e-9
    rock_feo: float | None         # FeO, mass fraction of the rock
    rock_al2o3: float | None       # Al2O3, mass fraction of the rock (Earth bulk silicate 0.045)
    rock_mg_si: float | None       # molar Mg/Si of the rock
    envelope_helium: float | None  # Y of the envelope (mass fraction of the envelope)
    envelope_z: float | None       # heavy elements, mass fraction of the envelope
    porosity: float | None         # small bodies only
    carbon_rock: bool
    origin: str                    # "drawn" (route 2) or "inverse" (route 1)

def bulk_fractions(body) -> BulkFractions | None   # None for a body saved before GEN.428 and a body with no composition
```

Body attributes and columns (planets and moons tables, DB.69): `metal_mass_fraction`, `rock_mass_fraction`, `water_mass_fraction` (exists today in hydrosphere; from GEN.428 it is the drawn value, not a separate draw), `ice_mass_fraction`, `envelope_mass_fraction`, `rock_feo_fraction`, `rock_al2o3_fraction`, `rock_mg_si_ratio`, `envelope_helium_fraction`, `envelope_z_fraction`, `porosity`, `carbon_rock` (bool), `composition_origin` (`drawn` | `inverse`). All dimensionless; declared in `db/column_units.py` in the DB.69 migration. The detailed mixture (minerals, oxides, ices by name, the rows of the shared materials table GEN.427) is a child table `planet_composition` / `moon_composition` (body id, component = materials key, mass_fraction) written in the same PR; summary columns are the contract for everyone else, the child table is for the description writer (GEN.430) and the surface (GEN.441).

For lane H (answers to its questions): silicate_mass_fraction = `rock_mass_fraction`; refractory_enrichment = `rock_al2o3_fraction / 0.045`; core mass fraction of the dry part = `metal_mass_fraction / (metal_mass_fraction + rock_mass_fraction)`; water = `water_mass_fraction`. Until PR 3 lands `bulk_fractions(body)` returns None: use your stand-ins when it does.

### 1.2 Atmosphere mixture: `planetgen.physics.mixture`

```python
def partial_pressures_pa(body) -> dict[str, float]    # gas-table formula ("CO2", "N2", "He", "NH3", "Na") -> Pa, zeros dropped
def total_pressure_pa(mix) -> float                   # sum
def mole_fractions(mix) -> dict[str, float]
def mean_molar_mass_kg_mol(mix) -> float              # sum x_i M_i, M_i from physics.gases
```

Keys are exactly the `physics/gases.py` formulas (case as the table writes them). Since DB.71 (PR 4b, migration 0146) the air is the child table `planet_gases` / `moon_gases` (body id, gas formula, `partial_pa`), so ANY gas of the table can be stored, and the balance gas is a named stored gas like the rest (GEN.424: helium, ammonia and sodium are rows, no unstored remainder). The body attribute is `gas_partials_pa` (formula to Pa); a column-named dict (a row) carries the same key once `db.store.with_gases` has attached it. The ten `p_<gas>_kpa` columns are gone (no shim). Lane G, F and the radiation code read only through `mixture.partial_pressures_pa(body)` / `atmosphere.partial_kpa(body, gas)` / `mixture.mean_molar_mass_kg_mol`. The air column leaves out a gas with no Rayleigh cross section (sodium) until the gas table has one.

Stored consistency (tested in PR 4b): `sum(partial_pressures_pa) == atmospheric_pressure` to 1e-9 relative; `atm_molar_density == mean_molar_mass_kg_mol` (name and unit kept, kg/mol); `scale_height` follows from R T_e / (M g). GEN.412's `bond_albedo`, `greenhouse_factor`, `albedo_basis` are lane D's and are read, not written, by this lane. Haze, dust and cloud columns (GEN.414, DB.68) are lane G's: use the design names `tau_haze`, `tau_dust`, `tau_cloud` (550 nm), `cloud_cover`, `aerosol_draw` (atmosphere-colour.md section 5), not `cloud_tau`.

## 2. Step contracts

### 2.1 GEN.429 kernels: `planetgen.physics.bulk_structure` (pure; as built)

All SI unless named. No draws, no state.

```python
structure(bulk, mass_kg, *, flux_earth=1.0, age_gy=4.5, temperature_k=None, x_fe=None, core_scale=1.0, prefer_solver=False) -> Structure
    # radius_km, density_g_cm3, gravity_g, escape_kms, moment_of_inertia, cmb_pressure_gpa, method
    # method: "closed-dry" | "closed-water" | "closed-small" | "closed-envelope" | "solver" | "mixing"
layered_structure(mass_earth, metal, rock, water, x_fe, core_scale, ...) -> dict          # the reference solver (iron, silicate, water)
radius_fit_earth(mass_earth, cmf, x_fe, core_scale) -> float                              # part B closed form, 0.1 to 20 Earth masses
water_world_radius_earth(mass_earth, cmf, water) -> float                                 # part C closed form, 0.3 to 20
small_body_radius_earth(mass_earth, cmf, water) -> float                                  # fit over 1.5e-4 to 0.3 Earth masses
envelope_radius_earth(mass_earth, envelope, flux_earth, age_gy, ...) -> float             # [R2741], 1 to 20
mixing_density_kg_m3(metal, rock, water, ice, porosity) -> float                          # below 1e21 kg
cmf_from_mass_radius(mass_earth, radius_earth, ...)                                  # route 1 inverse (the water fraction's inverse was deleted unused, TEST.525)
x_fe_from_feo_mass_fraction(feo) -> float                                                 # the mantle's Fe/(Fe+Mg) from FeO
```

Water and ice together under 0.1 percent of the mass are solved as dry. A giant (envelope over 0.5) raises `OutsideRange`: the giant radius law of GEN.34 is `physics/planets.py`'s until GEN.460.

### 2.2 GEN.428: `planetgen.physics.bulk_draw` (as built)

```python
host_ratios(star) -> HostRatios     # fe, si, al, ca, ni (atoms per Mg), c_over_o, y_envelope, z; from the star's abundance key
dry_partition(ratios, rock_feo, light_elements, carbon_rock, graphite_over_rock=0.0) -> DryPartition   # metal, rock, rock_feo, rock_al2o3, rock_mg_si
draw_terrestrial(body, radius_given, mass_given)   # sets the BulkFractions attributes, mass, radius, density
ROCK_KIND_TABLE                                    # gated GEN.472 table "bulk-rock-kind": silicate | carbide-graphite | carbide-graphite-by-class
```

A giant is drawn by `planetgen.physics.giant_draw.draw_giant` (GEN.460, GEN.642), not here. The body needs `planet_class`, `radius` (or `mass`), `star` (with `metallicity_fe_h`, `age`, `mass`, `galactic_center_dist_ly`) and `is_moon`. Draw order (a change is an epoch change): rock kind, water share, other ices (class D), FeO, core light elements, porosity (class D), graphite (carbon rock). The mass is then solved from the class radius; class S and a body given both a mass and a radius are worked back from them (`origin` is `inverse`). Columns (DB.69, schema 144): `metal_mass_fraction`, `rock_mass_fraction`, `ice_mass_fraction`, `envelope_mass_fraction`, `rock_feo_fraction`, `rock_al2o3_fraction`, `rock_mg_si_ratio`, `envelope_helium_fraction`, `envelope_z_fraction`, `porosity`, `carbon_rock`, `composition_origin`, with the existing `water_mass_fraction`; read them through `bulk_composition.bulk_fractions(body)`, which returns `None` for a row saved before schema 144.

Hook points in `generate_planet_properties` after PR 3 (so other lanes know where to add draws): (a) bulk composition (`bulk_draw`), (b) lane H `radiogenic.draw_content(planet, silicate_mass_fraction=..., refractory_enrichment=...)` immediately after (a), (c) `calculate_atmospheric_conditions` (lane D's albedo and greenhouse, then GEN.411), (d) hydrosphere (reads the water share the bulk draw made). Draws come from the ambient bound stream in that order; each module draws in a fixed order of its own.

Extension coming: giant interiors (GEN.641, GEN.642, DB.89, lane J) add layers and a core to the giant composition; the five whole-body fractions stay the totals.

### 2.3 GEN.452 part 2 + GEN.411 (PR 4)

`atmosphere.generate_atmosphere(planet)` is replaced by: recipes `RECIPES.draw_all(subject)` (GEN.478 tables) -> partial pressures at the surface (saturation caps from `physics.saturation`, cold-end freeze-out, GEN.422 order) -> total pressure = sum -> mean molar mass -> scale height -> text from the stored mixture (GEN.415 is later). `atm_density` and the molar-density range draw of today are gone. Epoch bump with GEN.452 part 2's stored freeze-out values.

## 3. Rules for everyone building against this

- Read composition only through the two accessors. No lane reads or writes `p_<gas>_kpa` after PR 4, or `composition` text for logic.
- New draws: after the hook they belong to, from the ambient stream; never reorder existing draws without telling this lane.
- New numeric columns: declare in `db/column_units.py` in the same PR. Epochs and migrations are numbered at merge (first to merge takes the next free).
- Files this lane edits (avoid parallel edits, tell me first): `physics/planets.py` (`generate_planet_properties`, `calculate_atmospheric_conditions`, `reconcile_zone_and_class`), `physics/atmosphere.py`, `generation/planet.py` (SERIALIZABLE_FIELDS), `db/store.py` (planet and moon field lists), `db/models.py`, `db/schema.sql`, `db/column_units.py`, `physics/hydrosphere.py` (the water draw line only), `physics/habitability_world.py` (molar mass table only), `tuning.py` (class density and atmosphere ranges). New: `physics/bulk_composition.py`, `physics/bulk_structure.py`, `physics/mixture.py`.

## 4. End-to-end draw order: the star to the PHI-6 (DOC.78)

One table for every lane. "Stream" says where the draws come from; "Draws" says whether the block consumes random numbers at all. A block that must not move later draws has a stream of its own, named by a key.

**The reproducibility rule.** Every draw goes through `util/draw.py`. Inside a unit of work (a sector fill, a scatter run, a backfill block) there is one ambient stream, so the order of draws inside it is part of the output: inserting, removing or reordering one draw shifts every later draw and is an epoch change. Draws inside a block run in a fixed order. A block either has a stream of its own or appends after the last draw of the block before it. Derived values draw nothing. A stored formation branch is never redrawn, and a recalculation after an orbit validation or an edit never redraws a stored composition (F10, GEN.422). A table's version is part of the order: adding, removing or reweighting a row changes every seed in its classes.

### 4.1 Per star (once)

| Step | What | Stream | Draws |
|---|---|---|---|
| 1 | Initial mass (Kroupa IMF, sector truncation) | ambient | yes |
| 2 | Age from the population's star formation history | ambient | yes |
| 3 | Redraw while mass and age give a collapsed object | ambient | yes |
| 4 | Evolve to the present state (`evolve_star`) | none | no |
| 5 | Classify: spectral letter, luminosity class | none | no |
| 6 | Habitable zone, Hill sphere, orbit, spin, activity | none | no |
| 7 | [Fe/H] (a horizontal-branch star keeps the one it drew for its colour) | `metallicity:{mass}:{age}:{distance}` (bright stars `metallicity:bright:{seed}`) | yes |
| 8 | Snow line | `snow_line:{mass}:{age}:{distance}` | yes |
| 9 | Disc of a pre-main-sequence star | `disk:{mass}:{age}` | yes |
| 10 | Star composition (`star_makeup`) | `abundance:{key}` | yes, its own stream |

Planned: GEN.553 moves [Fe/H] before the state (step 7 before 4) and gives each step of the star its own stream; GEN.651 wires [Fe/H] into the evolution; GEN.652 settles the two main-sequence tables. When GEN.553 lands, this table is replaced by its five steps (origin, birth, state, composition, motion).

### 4.2 Per body (planet, moon), in this order

| Block | Draws in it | Stream | Source of the order |
|---|---|---|---|
| Count and slots | how many bodies, spacing, belts, zone | ambient | `system.py` `_generate_planets` |
| Formation | born beyond or migrant, formation distance and temperature, ice set | ambient | section 2.2, bulk-composition.md section 7 |
| Bulk composition | water share, dry part, mass (decision 28), envelope fraction and loss | ambient (hook a) | section 2.2 above |
| Structure | radius, density, gravity, escape speed | none | GEN.429 |
| Radiogenic content | K, U, Th, Al-26 | ambient (hook b, straight after a) | lane H |
| Atmosphere | `RECIPES.draw_all` in the tables' order, then saturation and freeze-out (no draws) | ambient (hook c) | section 2.3 |
| Albedo and greenhouse | the two draws GEN.412 records | ambient (hook c) | lane D |
| Temperature | fixed point of pressure, greenhouse and temperature | none | GEN.422 |
| Hydrosphere, ozone, aerosols, magnetism, spin | each in its own module's order | ambient (hook d) | lanes F, J, L |
| Radiation | surface dose from the stored inputs | none | lanes G, H |
| Moons | recursion per orbit slot | planned `Stream(f"{host_seed}:orbit:{slot}")` (GEN.556) | GEN.556 |
| After the body is complete | `classify` (GEN.646), PHI-6 (GEN.643), the four life scores (GEN.649) | none, except the evolutionary-stage draw for the life scores | lane D, class lane |

Until GEN.646 lands the class letter is still an input at the count, formation and bulk blocks; afterwards it is only an output of `classify`.

### 4.3 Rules for changing the order

1. A new draw goes after the hook it belongs to, or on a stream of its own if it must not move later draws. Say so in the PR.
2. Moving a draw, or changing a table's rows or weights, is an epoch bump with a migration, numbered at merge.
3. Whoever changes a row in this table updates it in the same PR.

## 5. Boundaries agreed with other lanes

### L boundary (lane L, 2026-10-23): who writes the mixture

Agreed as lane L proposed (/mnt/project-files/reports/lane-l-boundary-2026-10-23.md). Lane L writes `planetgen/physics/mixture_draw.py`: `draw_mixture(subject) -> MixtureDraw(mole_fractions: {gas-table formula: fraction, sum 1}, regimes drawn, redox_diw)`, a pure function of a `recipes.Subject` (inputs named in `atmosphere_recipes.INPUTS`), with no pressure logic and no storage. Lane I (PR 4) keeps the total pressure (retention row d), the saturation caps and freeze-out (GEN.452, GEN.422 order), the sum of partial pressures, the mean molar mass, the storage and DB.71. PR 4 calls `draw_mixture` in place of `CLASS_MIXES` and deletes `CLASS_MIXES` in the same PR (GEN.479's done-line), so GEN.479 merges before PR 4 or in step with it. Lane L does not edit `planets.py`, `atmosphere.py` or `mixture.py`; until `draw_mixture` is on main PR 4 builds against its signature above.

### K boundary (lane K, 2026-10-23): the small-body density

Agreed. Lane K writes `physics/compaction.py` (`grain_density_kg_m3(mass_fractions_by_component)` (times one minus the porosity; `small_body_density_kg_m3` was deleted unused, TEST.525), `compacted_density_kg_m3(radius_km, rock_fraction)`, the 1e19 to 1e21 kg taper). `bulk_structure.structure` keeps the rocky, water and envelope branches and the solver above 1e21 kg; its `mixing` branch (below 1e21 kg) is the zero-pressure mixing rule with porosity (`mixing_density_kg_m3`) until `compaction.py` is on main, and the first PR after that replaces the branch with a call to lane K's functions and deletes `mixing_density_kg_m3`. Lane K's comet, belt and icy-body draws return `BulkFractions` with `porosity` set.
