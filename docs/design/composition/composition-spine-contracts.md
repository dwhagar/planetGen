# Composition spine: interface contracts (lane I, 2026-10-23)

Owner: lane I (thread "Composition spine for habitability"). Items: GEN.565 (PlanetaryObject part), GEN.557, GEN.428 with DB.69, GEN.429, GEN.452 part 2, GEN.411 with DB.71 and GEN.424. Everything after GEN.428 (GEN.430, 441, 442, 460, 461, 466 and the rest) belongs to the lanes named in the plan, not to this lane. Designs: docs/design/composition/composition-model.md, bulk-composition.md, atmosphere-composition.md, atmosphere-colour.md section 7.

Status: DRAFT 1. Names below are final unless a downstream lane objects in its first reply. The accessors of section 1 land in PR "spine-contracts" first (no behaviour change, no epoch), so lanes can import them today from that branch and from main once it merges.

## 0. Order of work (and one change to the order the plan gives)

The plan lists GEN.565 then GEN.557 then GEN.428. GEN.428, GEN.429 and GEN.411 do not need the PlanetaryObject class or the wanderer draw: they are functions of a plain `BulkSubject` / `AirSubject` (below), called from `generate_planet_properties` today and from the one body draw later. GEN.557 also waits on GEN.553 (lane A) and DB.73 (lane B). So the order is:

1. PR "spine-contracts": this document in docs/design, plus `physics/bulk_composition.py` (data class + accessor) and `physics/mixture.py` (accessors). No epoch.
2. PR GEN.429 kernels: `physics/bulk_structure.py`, pure functions (reference solver, closed forms, small-body mixing). No epoch.
3. PR GEN.428 + DB.69 (+ the derived radius, density and gravity of GEN.429 in the generator): epoch bump, migration, golden update.
4. PR GEN.452 part 2 (#2613 rebased) + GEN.411 + DB.71 + GEN.424: epoch bump, migration.
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

Keys are exactly the `physics/gases.py` formulas (case as the table writes them). Today it reads the ten `p_<gas>_kpa` columns (and maps `o2` to `O2`, and so on); from PR 4 it reads the child table `planet_gases` / `moon_gases` (body id, gas formula, partial_pa), so ANY gas of the table can be stored (DB.71) and the balance gas is a named stored gas like the rest (GEN.424: helium, ammonia and sodium become rows, no unstored remainder). The ten columns are dropped in PR 4 (no shim); lane G, F and the radiation code should read only through `partial_pressures_pa(body)` / `mixture.mean_molar_mass_kg_mol`, never `p_<gas>_kpa`.

Stored consistency (tested in PR 4): `sum(partial_pressures_pa) == atmospheric_pressure` to 1e-9 relative; `atm_molar_density == mean_molar_mass_kg_mol` (name and unit kept, kg/mol); `scale_height` follows from R T_e / (M g). GEN.412's `bond_albedo`, `greenhouse_factor`, `albedo_basis` are lane D's and are read, not written, by this lane. Haze, dust and cloud columns (GEN.414, DB.68) are lane G's: use the design names `tau_haze`, `tau_dust`, `tau_cloud` (550 nm), `cloud_cover`, `aerosol_draw` (atmosphere-colour.md section 5), not `cloud_tau`.

## 2. Step contracts

### 2.1 GEN.429 kernels: `planetgen.physics.bulk_structure` (PR 2, pure)

All SI unless named. No draws, no state.

```python
rocky_structure(mass_earth, cmf, x_fe, core_scale=1.0, p_surface_pa=1e5) -> RockyStructure
    # radius_earth, central_pressure_gpa, cmb_pressure_gpa, gravity_m_s2, escape_kms, moment_of_inertia_factor
    # (part B 3.8 verbatim, standard library; valid 0.1 to 20 Earth masses)
radius_fit_earth(mass_earth, cmf, x_fe, core_scale=1.0) -> float            # B closed form
water_world_radius_earth(mass_earth, cmf, water_fraction) -> float          # C closed form
envelope_radius_earth(mass_earth, envelope_fraction, flux_earth, age_gy) -> float   # Lopez and Fortney, tables from lf14.csv
small_body_density_kg_m3(fractions, mass_kg) -> float                       # 1/sum(w_i/rho_i), (1 - porosity), brittle compaction below the taper
structure(fractions, mass_kg, *, flux_earth, age_gy, temperature_k) -> Structure
    # radius_km, density_g_cm3, gravity_g, escape_kms, moment_of_inertia_factor, branch ("rocky"|"water"|"envelope"|"giant"|"small")
```

`structure` picks the branch from the fractions and the mass and stays inside each form's stated range; outside it uses the solver. Inverse for route 1: `fractions_from_mass_radius(mass_kg, radius_km, ...)` returns `BulkFractions` with `origin="inverse"`.

### 2.2 GEN.428: `planetgen.physics.bulk_composition.draw_bulk(subject, class_gate)` (PR 3)

```python
@dataclass
class BulkSubject:
    planet_class: str | None; is_moon: bool; zone: str
    formation_au: float; snow_line_au: float; born_beyond_snow_line: bool
    star_fe_h: float; star_age_gy: float; star_mass_sol: float; star_makeup: StarMakeup  # physics.stellar_composition
    host_mass_kg: float | None       # moon: its planet's mass
    mass_kg: float | None            # given or None (then the class mass range, flagged stand-in, decision 28)
    radius_km: float | None          # given or None
```

Draw order (a change is an epoch change; fixed): star ratios from `star_makeup` (once per star) -> born beyond or migrant (stored `born_beyond_snow_line`) -> water share -> dry part (Fe/Mg, Si/Mg, branch, w_FeO, core O and S) -> mass (decision 28) -> envelope fraction and loss -> derive structure (GEN.429). Rows come from gated GEN.472 tables registered as `bulk-*`; the class letter is an input to the gate while GEN.90 is open, never the cause. A recalculation never redraws a stored composition (F10).

Hook points in `generate_planet_properties` after PR 3 (so other lanes know where to add draws): (a) bulk composition, (b) lane H `radiogenic.draw_content(planet, silicate_mass_fraction=..., refractory_enrichment=...)` immediately after (a), (c) `calculate_atmospheric_conditions` (lane D's albedo and greenhouse, then GEN.411), (d) hydrosphere. Draws come from the ambient bound stream in that order; each module draws in a fixed order of its own.

### 2.3 GEN.452 part 2 + GEN.411 (PR 4)

`atmosphere.generate_atmosphere(planet)` is replaced by: recipes `RECIPES.draw_all(subject)` (GEN.478 tables) -> partial pressures at the surface (saturation caps from `physics.saturation`, cold-end freeze-out, GEN.422 order) -> total pressure = sum -> mean molar mass -> scale height -> text from the stored mixture (GEN.415 is later). `atm_density` and the molar-density range draw of today are gone. Epoch bump with GEN.452 part 2's stored freeze-out values.

## 3. Rules for everyone building against this

- Read composition only through the two accessors. No lane reads or writes `p_<gas>_kpa` after PR 4, or `composition` text for logic.
- New draws: after the hook they belong to, from the ambient stream; never reorder existing draws without telling this lane.
- New numeric columns: declare in `db/column_units.py` in the same PR. Epochs and migrations are numbered at merge (first to merge takes the next free).
- Files this lane edits (avoid parallel edits, tell me first): `physics/planets.py` (`generate_planet_properties`, `calculate_atmospheric_conditions`, `reconcile_zone_and_class`), `physics/atmosphere.py`, `generation/planet.py` (SERIALIZABLE_FIELDS), `db/store.py` (planet and moon field lists), `db/models.py`, `db/schema.sql`, `db/column_units.py`, `physics/hydrosphere.py` (the water draw line only), `physics/habitability_world.py` (molar mass table only), `tuning.py` (class density and atmosphere ranges). New: `physics/bulk_composition.py`, `physics/bulk_structure.py`, `physics/mixture.py`.
