-- stellarObjects/schema.sql
--
-- MySQL (InnoDB) schema for planetGen persistence (TODO.md Phase 2, ported
-- off SQLite for Phase 5's MySQL migration -- see docs/TODO.md). Plain SQL
-- DDL, no ORM. Column lists are verified against the actual generator
-- source (stellarObjects/config.py, starData.py, doubleStar.py,
-- planetData.py, asteroidData.py, spaceSector.py) as of schema version 2 --
-- not against TODO.md's earlier field-count estimates.
--
-- Searchable-field principle: the primary fields a caller will search on
-- are exactly the ones each object already exposes in its own generated
-- text -- the `*_properties` dict built in `to_paragraph_list()` for
-- stars/binaries/planets (`table_*` columns below), and the equivalent
-- prose facts for asteroid belts, which have no properties dict at all
-- (density, the distance range, and composition -- see the asteroid_belts
-- section below). Raw generative scalar columns (mass_kg, radius_km, etc.)
-- are kept alongside for full object-graph fidelity (Phase 1 merge
-- decision), but the `table_*`/prose-derived columns are the ones meant
-- for search, since they're exactly what's published.
--
-- Distance-unit convention -- TWO standard units, by scale:
--   - SECTOR-SCALE PLACEMENT (`_mpc` suffix): star_systems.position_x/y/z
--     and sectors.edge -- where a system sits within its sector. Kept in
--     milliparsecs rather than km specifically because this is the one
--     place km's numbers get unwieldy (an 11.5-ly sector edge is
--     ~1.09e14 km vs. ~3526 mpc) and because sector geometry has no other
--     unit competing for consistency the way orbital distances do (every
--     orbital-scale quantity already shares km with radius/volume, so
--     there's no benefit standardizing sector geometry to km too).
--   - EVERYTHING ELSE (`_km`/`_km3` suffix): every other distance/length
--     column -- both "distance between things" (orbital distances,
--     habitable zones, system perimeter, heliosphere radius, hill radius,
--     scale height, binary separation) AND "size of an object" (radius,
--     volume) -- is kilometers, one standard unit instead of the
--     generator's native mix of km/AU, so search/comparison across
--     system-scale quantities never requires unit-aware query logic.
-- Every `table_*` column (the rendered wiki-table snapshot strings) is
-- unaffected by either convention -- those are copies of already-formatted
-- display text (still AU/ly/km as the wiki shows today), independent of
-- the raw column's storage unit. Conversion helpers live in
-- stellarObjects/utils.py: `ly_to_milliparsecs`/`milliparsecs_to_ly` for
-- the sector-scale columns; AU-to-km needs no helper for the km columns,
-- it's a single multiply by the existing `physical_constants.AU_TO_KM`.
-- Generation/physics code is untouched either way and keeps using its own
-- native km/AU/ly throughout; only the persistence layer converts, at the
-- moment of writing to (or reading from) these columns.
--
-- star_systems.quadrant stores the sector octant label (Roman numerals
-- I-VIII, one of the 8 sign(x)/sign(y)/sign(z) combinations -- see
-- spaceSector.py's `classify_octant`/`program_constants.SECTOR_OCTANT_LABELS`)
-- alongside the raw x/y/z -- purely derived from position, but worth
-- storing (not just recomputing on read) so it's directly
-- queryable/indexable without a UDF or generated-column expression. The
-- column name is kept as `quadrant` for schema stability, but every
-- human-facing label now displays it as "Octant" (html/sector.py,
-- html/system.py, html/static/sectormap.js) so it doesn't collide with
-- the unrelated, galaxy-scale "Quadrant" concept `sectors.center_x/y/z_pc`
-- and html/lib/galaxymap.py introduced later (4 azimuthal regions
-- spanning many sectors, not this column's 8 sign-combination regions
-- within one sector's own cube).
--
-- v3: star_systems.location stores a human-readable "sector name + nearest
-- neighbors" summary, e.g. "Voranthis Kelmoor -- nearest: Alpha Prime
-- (4.2 ly), Beta Cerise (7.8 ly), Gamma Ost (9.1 ly)" -- the same
-- "derive once at write time, persist for queryability" treatment
-- `quadrant` gets above, so a system's neighborhood is directly searchable
-- without recomputing nearest-neighbor distances on every read. NULL under
-- exactly the same condition as `quadrant`/`position_x/y/z_mpc`: a system
-- never placed in a sector. Computed in `stellarObjects/_db.py` from
-- `SpaceSector.nearest_neighbors` (up to 3, nearest first) at
-- `insert_sector` time; `migrate_database`'s `_migrate_v1_to_v2`/
-- `_migrate_v2_to_v3` backfill it for a database migrated from an older
-- version, applying the same nearest-neighbor logic directly to the
-- already-stored rows instead of a live `SpaceSector` object.
--
-- Two independent version numbers, per TODO.md's "Schema versioning"
-- decision:
--   - `schema_migrations` (below) is the DDL-level structure version --
--     replaces SQLite's `PRAGMA user_version` (no MySQL equivalent) with a
--     real tracking table, per TODO.md's Phase 5 MySQL-migration note.
--   - star_systems.schema_version (per row) is the version of the
--     serialized object-graph shape (Phase 1's to_dict()) that produced it.
-- Both started at 1; the DDL version below is CURRENT_SCHEMA_VERSION in
-- `_db.py`.
--
-- v4: added six nullable galaxy-frame placement columns to `sectors` --
-- center_x/y/z_pc (Cartesian center, parsecs, galactic-origin-relative),
-- galactic_radius_pc (sqrt(x^2+y^2+z^2), persisted for the same
-- can't-index-an-expression reason star_systems.quadrant is persisted
-- rather than recomputed), and shell_index/shell_slot_index (the
-- sector's stable address within the shell/Fibonacci-sphere tiling
-- scheme -- see docs/design/galaxy-coordinate-system.md). All six are
-- NULL together for a sector never placed in a galaxy -- every sector
-- generated by sectorGen.py's own standalone (non-galaxy) tooling, and
-- any migrated from an older database -- enforced by the CHECK below on
-- the first four (shell_index/shell_slot_index are deliberately excluded
-- from that CHECK: they're bookkeeping for the placement *algorithm*, not
-- physically implied by a center point, so a sector could in principle
-- have a hand-authored galaxy position without this scheme's own shell
-- addressing).
--
-- v2: moons split out of `planets` into their own `moons` table (with a
-- `planet_id` FK to the planet they orbit), instead of self-referencing
-- via `planets.parent_planet_id`/`is_moon` -- the same "give it its own
-- table" treatment `star_systems`/`planets` already have, so a query or a
-- UI facet can tell a Class M planet from a Class M moon without an
-- `is_moon` filter. `planet_evolutionary_paragraphs`/
-- `planet_reflection_spectrum` gained moon-owned counterparts
-- (`moon_evolutionary_paragraphs`/`moon_reflection_spectrum`) for the same
-- reason. Moons never generate their own moons (`Planet.__init__` only
-- calls `generate_moons` `if not self.is_moon`), so `moons` has no
-- self-reference of its own. `stellarObjects/_db.py`'s `migrate_database`
-- converts an existing v1 database in place (backing up the original
-- first); `migrateDb.py` runs it over every database in a directory, and
-- `install.sh` (and so `update.sh`, which calls it) does this on every
-- deploy.
-- v5: dropped every pre-rendered `table_*`/`binary_table_*` TEXT column
-- (stars.table_type/table_radius/table_mass/table_temp/table_lum/table_hab/
-- table_loc; planets/moons.table_class/table_distance/table_period/
-- table_radius/table_gravity; star_systems.binary_table_*). These held a
-- display string (e.g. "{{Exp|4.20|5}} kg (1.00% of Sol)") rendered once at
-- generation time from `Star`/`Planet`/`BinaryStarProxy.get_table_properties()`,
-- baked in whichever of wikitext-template or HTML form `SystemConfig.MARKDOWN`
-- happened to be at insert time -- wrong for any consumer wanting the other
-- form (the interactive HTML viewer in particular, which was displaying raw
-- unrendered `{{Exp|...}}` wikitext template syntax). No data is lost: every
-- number these strings were derived from already has its own proper
-- numeric/text column on the same row (mass_kg, radius_km, temperature_k,
-- luminosity_w, habitable_zone_inner/outer_km, distance_km, period_years,
-- gravity_g, star_type, name, planet_class, binary_effective_mass_kg,
-- binary_effective_luminosity_w, binary_separation_km, binary_type, ...) --
-- see this file's own header comment above. Consumers now compute display
-- formatting on demand from those columns instead of reading a frozen
-- pre-rendered copy (`html/lib/tabledisplay.py` for the HTML viewer;
-- `get_table_properties()` is still used, unchanged, to build the actual
-- wiki-page text in `wikitext_content`/`markdown_content`).
-- v6/v7: give every galaxy-placed sector exact vertices -- built from a
-- local spherical Voronoi tessellation among its own same-shell
-- neighbors, extruded radially between the shell's inner and outer
-- bounding spheres (`stellarObjects/sectorGeometry.prism_vertices`) --
-- genuinely gap-free against same-shell neighbors (not merely reduced;
-- see that module's docstring for why exact circumcenters, not nudged
-- approximations, make this possible) and area-matched (not vertex-
-- matched) against the shells in front of and behind it. Vertex/face
-- count varies per sector (typically 5-7, not a fixed number) because it
-- has to: a cube tiling of a sphere cannot be gap-free in general (the
-- same reason a soccer ball needs pentagons mixed with hexagons).
-- Supersedes an earlier, never-released "relax a fixed 8-vertex cube
-- toward its neighbors" approach, which only approximately closed gaps.
--   - v6 stored this as `sectors.vertices_pc`, a JSON `{"inner": [...],
--     "outer": [...]}` blob -- reconsidered almost immediately in favor
--     of v7's plain relational table below, so no released version ever
--     depended on the JSON shape.
--   - v7 replaces `sectors.vertices_pc` with the `sector_vertices` table
--     (see that table's own comment) -- one row per vertex, ordinary
--     columns throughout, no serialized blob anywhere in the schema.
--
-- v8: the galaxy-wide density "skeleton" (`galaxyPlan.py`), stored as
-- compact structural facts rather than one row per sector. Position,
-- density, and vertices are pure deterministic functions of `(shell_index,
-- shell_slot_index)` and a handful of galaxy-wide shape parameters
-- (`galaxyDensity.GalaxyShape`, `sectorGeometry.prism_vertices`) -- none
-- of that needs to be stored per candidate sector, since it's cheaper to
-- recompute on demand than to look up. What genuinely needs precomputing
-- is *where the galaxy has any content at all*, which is why this adds:
--   - `galaxy_shape`: one singleton row holding the whole galaxy's shape
--     parameters, its calibration constant, its sector edge length, and
--     its outer edge (the last shell with any qualifying content) -- the
--     handful of numbers `stellarObjects.galaxyDensity`/`galaxySkeleton`
--     need to recompute any sector's exact position/density on demand.
--   - `galaxy_shell_band`: one row per contiguous *candidate* slot-index
--     band per shell (`stellarObjects.galaxySkeleton.find_shell_bands`) --
--     the safe (possibly slightly wider than exact) range of slots that
--     might hold qualifying content in that shell, found via the exact
--     closed-form phi<->slot-index inverse
--     (`galaxyGeometry.slot_index_bounds_for_phi_range`) already used
--     elsewhere in this codebase. Deliberately NOT one row per qualifying
--     sector, and deliberately NOT storing that sector's position/density/
--     vertices -- an individual slot's exact qualification is checked live
--     (cheap: a single `relative_density` evaluation) only when that slot
--     is actually visited (see `galaxyGen.ensure_sector_generated`), which
--     is also the only time its full content, vertices, and everything
--     else about it get generated and stored in `sectors`/`sector_vertices`
--     /the star-system tables -- lazily, once, never recomputed again.
--     Almost every shell has exactly one band (this density model is
--     symmetric about and peaks at the galactic plane for any realistic
--     parameter choice -- verified directly across a full real-Milky-Way-
--     scale build), but the schema allows more than one per shell rather
--     than assuming it, since nothing about the model rules it out for
--     every possible parameter choice.
--   - `sectors` also gains a `UNIQUE (shell_index, shell_slot_index)`
--     constraint (see that table's own comment) -- now that a sector can
--     be lazily generated the moment it's visited
--     (`galaxyGen.ensure_sector_generated`), rather than only ever
--     through a single-process batch script, two concurrent visits to the
--     same never-before-generated address are possible; this constraint
--     turns that race into a clear `IntegrityError` the visit path
--     recovers from (re-fetch and return the sector the other visit just
--     created) instead of a silent duplicate row at the same address.
--
-- v9: orbital motion. `planets`/`moons` each gain three fixed-at-
--   generation-time orbital-orientation columns (`orbital_inclination_deg`,
--   `orbital_ascending_node_deg`) plus one that changes over time
--   (`orbital_phase_deg`, this body's current position angle around its
--   otherwise-circular orbit) and one static descriptive stat
--   (`rotation_period_hours`, axial "day length" -- this generator doesn't
--   track rotational phase, only period). `updateOrbits.py` is a new,
--   separately-run script that advances every body's `orbital_phase_deg`
--   in place based on `period_years` and real elapsed time, using the new
--   `orbit_simulation_state` singleton row (one per database, same
--   pattern as `galaxy_shape` above) to track when it last ran --
--   `stellarObjects.galaxyGeometry`'s "derive, don't store" principle
--   doesn't apply here: the whole point of this feature is a value that
--   changes with real-world time rather than being a pure function of the
--   body's other (fixed) generative facts, so it has to be stored and
--   periodically updated, not recomputed on read.
--
-- v10: galactic orbit. `stars` gains `galactic_orbital_speed_kms` and
--   `galactic_orbital_period_gy` -- a star system's circular orbital speed
--   and period around the galactic center, derived from its (or the fixed
--   `physical_constants.GALACTIC_CENTER_DISTANCE_LY` fallback's) distance
--   from it via a simple rotation-curve model (see
--   `physical_constants.GALACTIC_ROTATION_FLAT_VELOCITY_KMS`'s comment and
--   `utils.calculate_galactic_orbit`) -- the same "fixed at generation
--   time, stored rather than re-derived on read" treatment
--   `system_perimeter_km`/`heliosphere_radius_km` already get on this same
--   table, since both ultimately depend on which sector (if any) the
--   system was placed in. `star_systems` gains the equivalent
--   `binary_galactic_orbital_speed_kms`/`binary_galactic_orbital_period_gy`
--   pair among its `binary_*` columns, NULL under the same
--   binary-only condition as `binary_system_perimeter_km`.
--
-- v11: planet/moon position. `planets`/`moons` each gain
--   `position_x_km`/`_y_km`/`_z_km` (this body's Cartesian position
--   relative to its orbital anchor -- the star, or a binary system's
--   combined center, for a planet; the parent planet for a moon) and
--   `orbital_speed_kms` (its constant circular-orbit speed). Position is
--   derived from `distance_km` and the v9 orbital-motion columns
--   (`orbital_inclination_deg`/`orbital_ascending_node_deg`/
--   `orbital_phase_deg`) via `utils.orbital_position_au` -- the same
--   "each body positioned relative to its immediate primary" convention
--   `docs/design/galaxy-coordinate-system.md` already uses one level up
--   for sectors/systems relative to the galactic center (see v10 above).
--   `updateOrbits.py`/`_db.advance_orbital_phases` recomputes position in
--   lockstep with `orbital_phase_deg` as time passes; `orbital_speed_kms`
--   only changes if `distance_km`/`period_years` themselves do (e.g.
--   `StarSystem.validate_system` resolving an orbital overlap at
--   generation time), never from phase advancing alone.
--
-- v12: floating-point update guard. `planets`/`moons` gain
--   `min_update_interval_years` -- the shortest `elapsed_years` worth
--   calling `_db.advance_orbital_phases` for. Below this, the phase delta
--   `elapsed_years` would add is smaller than `orbital_phase_deg`'s own
--   IEEE 754 double-precision resolution, so `MOD(orbital_phase_deg +
--   delta, 360)` is guaranteed to round right back to the exact value
--   already stored -- a wasted write that changes nothing. Derived purely
--   from `period_years` (`utils.minimum_update_interval_years`:
--   `period_years * math.ulp(360.0) / 360` -- the coarsest representable
--   step anywhere in `orbital_phase_deg`'s `[0, 360)` range, used as a
--   single conservative bound rather than a per-row value that would
--   itself need updating every time phase does), so it's fixed at
--   generation time and only changes if `period_years` itself does (e.g.
--   `StarSystem.validate_system` resolving an orbital overlap post-hoc,
--   the same trigger `orbital_speed_kms` recomputes on); phase advancing
--   alone never touches it. Not a narrative/display stat -- purely a
--   guard value `advance_orbital_phases` reads to skip a row's `UPDATE`
--   entirely when a call's `elapsed_years` wouldn't move it. Scoped to
--   `planets`/`moons` only: `stars`' galactic-orbit values are fixed
--   forever at generation time (no periodic update mechanism exists for
--   them the way `advance_orbital_phases` exists for `orbital_phase_deg`),
--   so there's nothing for this guard to protect there.
--
-- v13: star motion -- galactic orbit phase, plus binary mutual orbit.
--   Supersedes v12's stars scoping note above: stars now DO have a
--   periodic update mechanism, so they need the same guard planets/moons
--   already have.
--     `stars` gains `galactic_orbital_phase_deg` (this star's current
--   angular position around its galactic orbit -- the same role
--   `orbital_phase_deg` plays for a planet/moon, advanced by
--   `_db.advance_orbital_phases` based on `galactic_orbital_period_gy`) and
--   `galactic_min_update_interval_years` (the same floating-point update
--   guard v12 added for planets/moons, `utils.minimum_update_interval_years`
--   applied to `galactic_orbital_period_gy * 1e9` years instead of
--   `period_years`). Both stars of a binary pair -- and `star_systems`'
--   own `binary_galactic_orbital_phase_deg`/
--   `binary_galactic_min_update_interval_years` pair, added alongside the
--   existing `binary_galactic_orbital_speed_kms`/`_period_gy` -- always
--   carry the identical phase: a binary's AU-scale separation is
--   negligible next to its light-year-scale galactic orbit radius, so the
--   pair moves around the galaxy together, not independently
--   (`StarSystem.__init__` rolls the phase once and threads it to both
--   stars and the proxy -- see that method's docstring).
--     `star_systems` also gains the binary pair's own MUTUAL orbit --
--   entirely separate from (and vastly faster than) the galactic orbit
--   above: the two stars circling their common barycenter, NULL under the
--   same binary-only condition as every other `binary_*` column.
--   `binary_mutual_orbital_period_years`/`_speed_kms` are Kepler's third
--   law and the standard circular-orbit speed formula
--   (`planetPhysics.calculate_orbital_period_years`/
--   `utils.circular_orbital_speed_kms`) applied to the pair's already-
--   stored `binary_separation_km`/`binary_effective_mass_kg` -- the same
--   formulas a planet's orbit around its star already uses, just with the
--   combined pair mass standing in for "the primary". `_inclination_deg`/
--   `_ascending_node_deg`/`_phase_deg` fully orient this mutual orbit in
--   3D and track the pair's current position within it (the "orbital
--   direction" a binary's own orbital plane needs, unlike the galactic
--   orbit above, which this generator treats as planar) -- the same
--   `utils.orbital_position_au` orbital-element convention planets/moons
--   already use, just with no small-tilt bias the way
--   `PLANET_ORBITAL_INCLINATION_MAX_DEG` gives a planet's protoplanetary-
--   disk-derived orbit: a binary pair's mutual orbital plane has no
--   preferred alignment, so `_inclination_deg` is drawn from the full
--   `[0, 180)` range. `binary_mutual_min_update_interval_years` is the
--   same guard formula again, applied to
--   `binary_mutual_orbital_period_years`. `_db.advance_orbital_phases`
--   advances `binary_mutual_orbital_phase_deg` the same way it advances
--   `orbital_phase_deg`/`galactic_orbital_phase_deg` elsewhere, guarded by
--   its own interval.
--
-- v14: binary mutual orbit position. `star_systems` gains
--   `binary_mutual_position_x_km`/`_y_km`/`_z_km` -- the secondary's
--   Cartesian position relative to the primary, derived from
--   `binary_separation_km` and the v13 `binary_mutual_orbital_
--   {inclination,ascending_node,phase}_deg` columns via
--   `utils.orbital_position_au`, the same "each body positioned relative
--   to whatever it actually orbits" convention `planets`/`moons.
--   position_x/y/z_km` already use one level down (a planet relative to
--   its star, a moon relative to its parent planet -- see the v11 note
--   above) and `sectors`/`star_systems` use one level up (relative to the
--   galactic center). `_db.advance_orbital_phases` recomputes this in
--   lockstep every time `binary_mutual_orbital_phase_deg` advances, mirroring
--   exactly how it already keeps a planet's/moon's position in lockstep
--   with its own phase -- position has no independent update of its own,
--   it just has to move whenever phase does. NULL under the same
--   binary-only condition as every other `binary_*` column.
--
-- v15: S-type (wide) binary support (planetGen's second real binary
--   configuration alongside the P-type/close pair every `binary_*` column
--   above already covered). `star_systems` gains `binary_configuration`
--   (`'close'` | `'wide'` | NULL for a single star -- the authoritative
--   discriminator going forward; `is_binary` is kept, and now means "this
--   system has two stars", true for either configuration) plus
--   `binary_eccentricity`/`binary_periapsis_km`/`binary_apoapsis_km` (the
--   pair's own orbital eccentricity and periapsis/apoapsis separation --
--   real and used for `utils.holman_wiegert_critical_semimajor_axis`, but
--   NOT used to make the pair's live position/phase-advance tracking
--   eccentric -- see `doubleStar.WideBinaryPair`'s module docstring for why
--   that stays the existing circular approximation). These three are 0/
--   `binary_separation_km` for a `'close'` pair (tidal circularization
--   makes near-zero eccentricity a legitimate simplification there) and
--   real, generally non-zero values for a `'wide'` one. The existing
--   `binary_separation_km`/`binary_mutual_orbital_*`/
--   `binary_mutual_position_*_km` columns are reused UNCHANGED for a wide
--   pair's own mutual orbit -- they already represent exactly what a wide
--   pair's (circularly-approximated) mutual orbit needs, no new columns
--   for that part. NULL for a `'wide'` pair, unlike a `'close'` one (no
--   merged effective star exists to describe): `binary_type` (the
--   pre-existing column -- note this is NOT the same thing as the new
--   `binary_configuration` above; `binary_type` holds the close pair's
--   merged spectral-summary string, e.g. "Binary (G/K)", a name chosen
--   before S-type support existed and kept as-is for that pre-existing
--   meaning rather than reused, to avoid a silent meaning change for
--   existing data), `binary_temperature_k`, `binary_radius_km`,
--   `binary_effective_mass_kg`, `binary_effective_luminosity_w`,
--   `binary_age_gy`, `binary_lifespan_gy`,
--   `binary_habitable_zone_inner_km`/`_outer_km`, `binary_system_perimeter_km`,
--   `binary_heliosphere_radius_km`, `binary_galactic_orbital_*` (four
--   columns) -- a wide binary's two stars already carry their own galactic
--   orbit and habitable-zone columns individually on their own `stars`
--   rows, so nothing here needs a merged/combined equivalent.
--     `stars` gains `wide_binary_a_crit_km` -- this star's own Holman &
--   Wiegert (1999) critical semi-major axis (`utils.
--   holman_wiegert_critical_semimajor_axis`), the maximum orbit distance
--   that stays long-term stable given its companion's perturbation. NULL
--   for a single star or either constituent of a `'close'` pair (whose
--   planets orbit the merged proxy instead, with no per-star limit of this
--   kind), populated for both stars of a `'wide'` pair.
--     `asteroid_belts` gains `star_id`, the same nullable "which specific
--   star this orbits" column `planets`/`moons` already have (see
--   `planets.star_id`'s own comment below) -- needed because a wide
--   binary's two stars can now each have their own asteroid belts, not
--   just their own planets.
--     `planets.star_id`'s own comment (obsolete as of this version) is
--   updated below: a wide binary's planets DO now populate `star_id` with
--   the specific star (primary or secondary) they orbit, unlike a `'close'`
--   pair's (still NULL, since a close pair's planets orbit the merged
--   proxy, never one constituent star).
--
-- v16: exotic stellar phenomena (phenomenonGen.py's separate, rarer
--   generation mode -- see that module's docstring; these are never
--   produced by systemGen.py/sectorGen.py's normal generation odds). Six
--   new tables, none altering any pre-existing table's shape:
--     `black_holes`/`neutron_stars` -- satellite tables extending a
--   `stars` row (`star_id`, nullable) the same way `asteroid_belt_composition`
--   extends `asteroid_belts` -- populated only when `phenomenonGen.py
--   --anchor-system` wraps the compact remnant in a full `StarSystem`
--   (`compactRemnant.py`'s `BlackHole`/`NeutronStar`, which set the
--   owning `stars.yerkes_class` to the literal marker `'BH'`/`'NS'`
--   rather than a real Yerkes class). `star_id` is NULL for a remnant
--   generated standalone (no owning `StarSystem` at all).
--     `nebulae`/`supernova_remnants`/`rogue_planets`/`interstellar_comets`
--   -- always standalone (nothing in this generator places these within a
--   `StarSystem`), with a nullable `sector_id` reserved for a future
--   sector-context encounter (unused by `phenomenonGen.py` today, which
--   never creates or attaches a sector). `supernova_remnants` additionally
--   references at most one of `black_holes`/`neutron_stars` (`compact_remnant_kind`
--   discriminator, mirroring `StarSystem.binary_type`'s own "type string
--   selects which nullable reference is populated" convention) for a
--   core-collapse remnant whose collapsed core is still detectable within
--   it -- always both NULL for a Type Ia remnant, which leaves nothing
--   behind. `interstellar_comets` has its own `interstellar_comet_composition`
--   child table, mirroring `asteroid_belt_composition`'s per-component
--   breakdown (but without a concentration level -- `InterstellarComet`'s
--   own composition list carries no such gradient).
--   These six tables are new tables only -- a database migrated from an
--   older version needs no `ALTER TABLE` for them (see `_db._migrate_v15_to_v16`),
--   since `_ensure_schema`'s `CREATE TABLE IF NOT EXISTS` (run on every
--   connection) already creates them directly from this file regardless
--   of the database's recorded `schema_migrations` version.
--
-- v17: galactic-orbital motion for every standalone exotic phenomenon,
--   plus a seventh phenomenon, standalone asteroid fields. A rogue planet,
--   interstellar comet, nebula, supernova remnant, or asteroid field is
--   unbound from any specific STAR, not from the galaxy itself -- it still
--   orbits the galactic center on the same timescale a lone star does, via
--   the same mass-independent rotation-curve model `stars.
--   galactic_orbital_phase_deg` already uses (see `utils.
--   generate_galactic_orbit_fields`). `black_holes`/`neutron_stars`/
--   `nebulae`/`supernova_remnants`/`rogue_planets`/`interstellar_comets`
--   each gain `galactic_orbital_speed_kms`, `galactic_orbital_period_gy`,
--   `galactic_orbital_phase_deg`, `galactic_min_update_interval_years` --
--   nullable on `black_holes`/`neutron_stars` (only meaningful/populated
--   when `star_id IS NULL`; an anchored remnant's motion already lives on
--   its owning `stars` row instead, avoiding two sources of truth for the
--   same anchored object's position), `NOT NULL` on the other four (always
--   standalone, always populated). Unlike v16's brand-new tables, these
--   are new columns on already-existing tables, so (unlike
--   `_migrate_v15_to_v16`) `_migrate_v16_to_v17` needs real `ALTER TABLE`
--   statements, not just a bookkeeping row.
--     New table `asteroid_fields` (plus child `asteroid_field_composition`,
--   mirroring `asteroid_belt_composition`/`interstellar_comet_composition`)
--   -- `asteroidFieldData.AsteroidField`, physically the same object as
--   `asteroid_belts` (density + mineral composition, generated via the
--   same shared `asteroidData.generate_asteroid_composition`/
--   `format_composition_summary` helpers) but standalone, drifting in open
--   space rather than orbiting a star -- same "always standalone, nullable
--   `sector_id` reserved" shape as `nebulae` above, plus the four
--   `galactic_orbital_*` columns from the start (a new table, so no
--   `ALTER TABLE` needed for it specifically).
--
-- v18: proper two-body (barycentric) trajectories for binary stars,
--   star<->planet, and planet<->moon -- previously each pair modeled only
--   "the lighter body orbits a fixed primary," which is a poor
--   approximation for a binary's secondary (sampled at 0.1-0.8x the
--   primary's mass) and for a moon near the generator's own mass cap
--   (up to 1/10 its parent planet's mass, approaching real "double
--   planet" ratios like Pluto/Charon). Existing "relative position"
--   columns (`planets.position_x/y/z_km`, `moons.position_x/y/z_km`,
--   `star_systems.binary_mutual_position_x/y/z_km`) are UNCHANGED --
--   still the true separation a large amount of existing physics
--   (insolation, Hill sphere, tidal locking) depends on. New columns add
--   the ORBITED body's own small "reflex offset"/"wobble" away from its
--   nominal fixed point instead -- see `utils.calculate_reflex_offset`'s
--   docstring for the formula.
--     `stars` gains `reflex_offset_x/y/z_km` -- a star's own displacement
--   from the combined pull of every planet orbiting it directly
--   (`planets.star_id`); class-blind, so applies identically to an
--   anchored `black_holes`/`neutron_stars` row. NULL/0 with no planets.
--     `planets` gains `reflex_offset_x/y/z_km` -- a planet's own
--   displacement from the combined pull of its own moons (`moons.planet_id`).
--   NULL/0 with no moons; not applicable to `asteroid_belts` (its own
--   table, never modeled with a gravitating mass anywhere in this codebase).
--     `star_systems` gains `binary_primary_position_x/y/z_km` and
--   `binary_secondary_position_x/y/z_km` -- each binary member's own
--   offset from the pair's barycenter, mirroring the pre-existing
--   `binary_mutual_position_*` (still the secondary relative to the
--   primary, unchanged); `secondary_position = primary_position +
--   binary_mutual_position` always holds, but both are stored explicitly
--   the same "ready-to-query, not derived on read" convention
--   `binary_mutual_position_*` itself already set. Also gains
--   `binary_secondary_mass_fraction` (the constant `secondary_mass /
--   (primary_mass + secondary_mass)`, stored so `_db.advance_orbital_phases`
--   never needs to join back to `stars` for either mass) and
--   `binary_planetary_wobble_x/y/z_km` -- the additional pair-wide wobble
--   from CIRCUMBINARY (P-type) planets, which orbit the merged proxy
--   (`star_id IS NULL`), not either individual star; their pull is
--   modeled as one shared wobble applied to the whole pair (splitting it
--   unevenly between primary/secondary would require solving a real
--   3+-body problem, out of scope) -- NULL unless `binary_configuration =
--   'close'`.
--     These are new columns on already-existing tables (`star_systems`,
--   `stars`, `planets`), so `_migrate_v17_to_v18` needs real `ALTER TABLE`
--   statements, the same as v17's column additions.
--
-- MySQL port -- type mapping and idempotency notes (TODO.md Phase 5):
--   - SQLite's `INTEGER PRIMARY KEY` (a 64-bit rowid alias) becomes
--     `BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY` throughout, with every
--     referencing FK column matching that same `BIGINT UNSIGNED` type --
--     InnoDB requires matching column types/signedness for a foreign key,
--     and BIGINT headroom means this never needs revisiting even at
--     Phase 4 galaxy scale (tens of millions of stars/planets/moons).
--   - SQLite's `REAL` becomes `DOUBLE` (MySQL's `REAL` is a synonym for
--     `DOUBLE` under default SQL modes, but this spells it out rather than
--     relying on a server-side mode setting).
--   - SQLite's untyped `TEXT` splits by actual content length: short
--     identifier/enum-like/name fields become `VARCHAR(n)` (indexable,
--     and MySQL's utf8mb4 index prefix limits make an unbounded TEXT index
--     impractical anyway); genuinely long generated prose (descriptions,
--     flavor text, composition summaries, evolutionary paragraphs) stays
--     `TEXT`; the two full-page rendered columns (`wikitext_content`/
--     `markdown_content`) that can exceed `TEXT`'s 64KiB ceiling for a
--     large system become `LONGTEXT`.
--   - `created_at` becomes a native `TIMESTAMP NOT NULL DEFAULT
--     CURRENT_TIMESTAMP` instead of `TEXT DEFAULT CURRENT_TIMESTAMP`.
--   - Every `CHECK` constraint carries over unchanged (MySQL enforces
--     `CHECK` from 8.0.16 -- this schema targets MySQL 8.0.16+ /
--     compatible MariaDB); multi-column CHECKs (`sectors`) are equally
--     supported.
--   - SQLite's column-level `REFERENCES` shorthand is parsed but NOT
--     enforced by MySQL/InnoDB -- every foreign key here is instead a
--     named table-level `CONSTRAINT ... FOREIGN KEY (...) REFERENCES
--     ...` clause, which InnoDB does enforce.
--   - SQLite's separate `CREATE INDEX IF NOT EXISTS` statements don't
--     translate directly: MySQL's `CREATE INDEX` has no `IF NOT EXISTS`
--     form, which would make re-running this script against an existing
--     database fail with "Duplicate key name" on the second run. Every
--     index is instead declared inline as a `KEY` clause inside its
--     table's `CREATE TABLE IF NOT EXISTS` -- since the whole table
--     (indexes included) is only ever created once, this script stays
--     safe to run repeatedly against an existing database, matching the
--     SQLite schema's original idempotency guarantee.
--   - MySQL has no `CREATE VIEW IF NOT EXISTS`; `sector_objects` below
--     uses `CREATE OR REPLACE VIEW` instead -- equally idempotent (and
--     "replace" is the more correct intent for a view, whose definition
--     might legitimately need updating across a migration in a way a
--     physical table's columns wouldn't).
--   - Every table is explicit `ENGINE=InnoDB` (foreign keys require it)
--     and `utf8mb4`/`utf8mb4_unicode_ci` (full Unicode, including
--     generated names/flavor text outside the Basic Multilingual Plane,
--     rather than MySQL's legacy 3-byte `utf8`).

-- ---------------------------------------------------------------------
-- schema_migrations -- DDL-level structure version tracking. Replaces
-- SQLite's `PRAGMA user_version` (see the MySQL port note above);
-- `_db.py`'s `migrate_database` reads `MAX(version)` from this table
-- and inserts a row per migration step it applies.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS schema_migrations (
    version      INT NOT NULL PRIMARY KEY,
    applied_at   TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- sectors
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS sectors (
    id                  BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    name                VARCHAR(255) NOT NULL,
    edge_mpc            DOUBLE NOT NULL,

    -- Galaxy-frame placement (v4) -- see the header comment's "v4" note.
    -- NULL together: a sector never placed in a galaxy (today's
    -- sectorGen.py standalone tooling, or migrated from a pre-v4
    -- database).
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    galactic_radius_pc  DOUBLE,

    -- The sector's stable (shell, slot) address within the shell/
    -- Fibonacci-sphere tiling scheme -- deliberately not part of the
    -- CHECK below (see the header comment's "v4" note).
    shell_index         INT,
    shell_slot_index    INT,

    CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    -- v8: at most one sector per (shell_index, shell_slot_index) address --
    -- see the header comment's "v8" note (galaxyGen.ensure_sector_generated
    -- relies on this to make a lazy-generation race produce a clear
    -- IntegrityError rather than a silent duplicate sector at the same
    -- address). NULL-together sectors (never placed in a galaxy) don't
    -- collide with each other or with a placed sector -- ordinary SQL NULL
    -- semantics for UNIQUE.
    UNIQUE (shell_index, shell_slot_index),

    KEY idx_sectors_galactic_radius_pc (galactic_radius_pc),
    KEY idx_sectors_shell_index (shell_index)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- This sector's own exact prism vertices (v7 -- see the header comment's
-- "v6/v7" note), one row per vertex rather than a serialized blob:
-- `sectorGeometry.prism_vertices` returns a variable number of vertices
-- per sector (typically 5-7, not fixed), split into an "inner" ring (on
-- the shell's inner bounding sphere) and an "outer" ring (on its outer
-- bounding sphere), both in the same cyclic order. `vertex_index` is that
-- cyclic position (0-based) within its own ring, not a global ordering --
-- pairing `(sector_id, vertex_index)` across the two rings gives the
-- lateral edge each pair of inner/outer vertices spans. No row exists for
-- a sector never placed in a galaxy (mirrors, at the application level
-- rather than a cross-table CHECK -- SQLite can't express "rows exist in
-- another table" as a CHECK constraint -- the same NULL-together
-- condition `sectors`'s own galaxy-placement columns enforce directly).
CREATE TABLE IF NOT EXISTS sector_vertices (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id     BIGINT UNSIGNED NOT NULL,
    ring          VARCHAR(8) NOT NULL CHECK (ring IN ('inner', 'outer')),
    vertex_index  INT NOT NULL,
    x_pc          DOUBLE NOT NULL,
    y_pc          DOUBLE NOT NULL,
    z_pc          DOUBLE NOT NULL,

    CONSTRAINT fk_sector_vertices_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    UNIQUE (sector_id, ring, vertex_index),
    KEY idx_sector_vertices_sector_id (sector_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- galaxy_shape / galaxy_shell_band -- the galaxy-wide density "skeleton"
-- (v8 -- see the header comment's "v8" note). `galaxy_shape` is a
-- singleton (`id` pinned to 1, enforced by the CHECK) -- there is exactly
-- one galaxy. Building or rebuilding the skeleton (`galaxyPlan.py`)
-- replaces this row and every `galaxy_shell_band` row wholesale; neither
-- table is ever partially updated.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS galaxy_shape (
    id                          BIGINT UNSIGNED PRIMARY KEY CHECK (id = 1),

    -- galaxyDensity.GalaxyShape fields, verbatim -- see that module for
    -- what each one means and how `k_norm` is calibrated.
    disk_scale_length_pc        DOUBLE NOT NULL,
    disk_scale_height_pc        DOUBLE NOT NULL,
    bulge_scale_radius_pc       DOUBLE NOT NULL,
    bulge_amplitude             DOUBLE NOT NULL,
    arm_count                   INT NOT NULL,
    pitch_angle_rad             DOUBLE NOT NULL,
    arm_amplitude               DOUBLE NOT NULL,
    spiral_reference_radius_pc  DOUBLE NOT NULL,
    spiral_reference_angle_rad  DOUBLE NOT NULL,
    k_norm                      DOUBLE NOT NULL,

    -- The sector edge length this skeleton was built at, in parsecs
    -- (galaxyGeometry's own scope -- every shell/sector-address formula
    -- takes this as a parameter rather than assuming a fixed constant).
    edge_pc                     DOUBLE NOT NULL,

    -- SpaceSector(edge_ly=...).expected_system_count() at relative_density
    -- = 1 -- cached because every qualification check
    -- (predicted_star_count >= 1.0) needs it, and it's a fixed galaxy-wide
    -- fact, not worth re-deriving from edge_pc on every call.
    expected_system_count_at_density_1  DOUBLE NOT NULL,

    -- The last shell index with any qualifying content -- this galaxy's
    -- real edge, discovered by `galaxyPlan.py` (a run of consecutive empty
    -- shells beyond it), not picked as an arbitrary radius.
    outer_shell_index           INT NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- One contiguous candidate slot-index band per shell (almost always
-- exactly one -- see the header comment's "v8" note). `band_index` orders
-- multiple bands within the same shell (0-based); a shell with no
-- qualifying content at all has no rows here.
CREATE TABLE IF NOT EXISTS galaxy_shell_band (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    shell_index      INT NOT NULL,
    band_index       INT NOT NULL,
    slot_index_min   INT NOT NULL,
    slot_index_max   INT NOT NULL,

    UNIQUE (shell_index, band_index),
    KEY idx_galaxy_shell_band_shell_index (shell_index)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Singleton row (same pattern as galaxy_shape above) tracking when
-- updateOrbits.py last advanced every planet's/moon's orbital_phase_deg in
-- this database, so the next run knows how much real time has actually
-- elapsed since then. Absent entirely until updateOrbits.py's first run
-- against a given database (it creates this row itself).
CREATE TABLE IF NOT EXISTS orbit_simulation_state (
    id               BIGINT UNSIGNED PRIMARY KEY CHECK (id = 1),
    last_updated_at  TIMESTAMP NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- system_configs -- one row per SystemConfig "recipe"
-- (stellarObjects/config.py:29-33 SERIALIZABLE_FIELDS)
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS system_configs (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    markdown          TINYINT(1) NOT NULL DEFAULT 0 CHECK (markdown IN (0, 1)),
    habitable_world   TINYINT(1) CHECK (habitable_world IN (0, 1)),
    asteroid_belt     TINYINT(1) CHECK (asteroid_belt IN (0, 1)),
    large_star        TINYINT(1) CHECK (large_star IN (0, 1)),
    moons             TINYINT(1) CHECK (moons IN (0, 1)),
    max_planets       TINYINT(1) CHECK (max_planets IN (0, 1)),
    planets           TINYINT(1) CHECK (planets IN (0, 1)),
    star_type         VARCHAR(64),
    name              VARCHAR(255),
    age               VARCHAR(16) CHECK (age IN ('young', 'old')),
    intelligent_life  TINYINT(1) CHECK (intelligent_life IN (0, 1)),
    binary_system     TINYINT(1) CHECK (binary_system IN (0, 1)),
    num_orbits        INT
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Child table for the variable-length SLOTS recipe list
-- (config.py:135-151).
CREATE TABLE IF NOT EXISTS system_config_slots (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    config_id     BIGINT UNSIGNED NOT NULL,
    orbit_index   INT NOT NULL,
    type          VARCHAR(16) CHECK (type IN ('planet', 'asteroid_belt')),
    planet_class  VARCHAR(16),
    moons         INT,

    CONSTRAINT fk_system_config_slots_config
        FOREIGN KEY (config_id) REFERENCES system_configs(id) ON DELETE CASCADE,
    KEY idx_system_config_slots_config_id (config_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- star_systems
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS star_systems (
    id                     BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id              BIGINT UNSIGNED,
    system_config_id       BIGINT UNSIGNED NOT NULL,
    name                   VARCHAR(255) NOT NULL,

    -- Position within its sector (spaceSector.py SectorSystemEntry.position),
    -- milliparsecs, relative to the sector's cubic center. NULL iff this
    -- system was never placed in a sector.
    position_x_mpc          DOUBLE,
    position_y_mpc          DOUBLE,
    position_z_mpc          DOUBLE,

    -- Sector octant label derived from the position above (Roman numeral
    -- I-VIII, see the header comment) -- NULL iff position is NULL.
    quadrant                VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- Human-readable "sector name + nearest neighbors" summary, see the
    -- header comment's "v3" note -- NULL iff position is NULL.
    location                TEXT,

    -- v15: now true for EITHER binary configuration (see binary_configuration
    -- below), not just a 'close' (P-type) pair.
    is_binary              TINYINT(1) NOT NULL DEFAULT 0 CHECK (is_binary IN (0, 1)),

    -- v15: 'close' (P-type/circumbinary, doubleStar.BinaryStarProxy) |
    -- 'wide' (S-type, doubleStar.WideBinaryPair) | NULL for a single star --
    -- see the header comment's "v15" note for why this is a distinct column
    -- from the pre-existing binary_type below, not a repurposing of it.
    binary_configuration   VARCHAR(8) CHECK (binary_configuration IN ('close', 'wide')),

    -- BinaryStarProxy-derived fields (doubleStar.py) -- all NULL for a
    -- single-star system or a 'wide' binary (see the header comment's "v15"
    -- note), stored rather than re-derived since _effective_mass/
    -- _effective_luminosity are computed once at generation time.
    binary_separation_km            DOUBLE,
    -- v15: the pair's own orbital eccentricity and periapsis/apoapsis
    -- separation -- 0/binary_separation_km for a 'close' pair, real
    -- (generally non-zero) values for a 'wide' one. See the header
    -- comment's "v15" note.
    binary_eccentricity              DOUBLE,
    binary_periapsis_km              DOUBLE,
    binary_apoapsis_km               DOUBLE,
    binary_type                     VARCHAR(64),
    binary_temperature_k            DOUBLE,
    binary_radius_km                DOUBLE,
    binary_effective_mass_kg        DOUBLE,
    binary_effective_luminosity_w   DOUBLE,
    binary_age_gy                   DOUBLE,
    binary_lifespan_gy              DOUBLE,  -- NULL = float('inf')
    binary_habitable_zone_inner_km  DOUBLE,
    binary_habitable_zone_outer_km  DOUBLE,
    binary_system_perimeter_km      DOUBLE,
    binary_heliosphere_radius_km    DOUBLE,
    binary_galactic_orbital_speed_kms   DOUBLE,
    binary_galactic_orbital_period_gy   DOUBLE,
    binary_galactic_orbital_phase_deg          DOUBLE,  -- v13, see header comment
    binary_galactic_min_update_interval_years  DOUBLE,  -- v13, see header comment

    -- v13 (see header comment): the pair's own mutual orbit -- separate
    -- from the galactic_* columns above.
    binary_mutual_orbital_period_years         DOUBLE,
    binary_mutual_orbital_speed_kms            DOUBLE,
    binary_mutual_orbital_inclination_deg      DOUBLE,
    binary_mutual_orbital_ascending_node_deg   DOUBLE,
    binary_mutual_orbital_phase_deg            DOUBLE,
    binary_mutual_min_update_interval_years    DOUBLE,
    -- v14 (see header comment): the secondary's position relative to the
    -- primary, kept in lockstep with binary_mutual_orbital_phase_deg --
    -- same "position relative to whatever this orbit is around" convention
    -- as planets/moons' own position_x/y/z_km.
    binary_mutual_position_x_km   DOUBLE,
    binary_mutual_position_y_km   DOUBLE,
    binary_mutual_position_z_km   DOUBLE,
    -- v18 (see header comment): each binary member's own offset from the
    -- pair's barycenter -- secondary_position = primary_position +
    -- binary_mutual_position above always holds, but both are stored
    -- explicitly rather than derived on read.
    binary_primary_position_x_km    DOUBLE,
    binary_primary_position_y_km    DOUBLE,
    binary_primary_position_z_km    DOUBLE,
    binary_secondary_position_x_km  DOUBLE,
    binary_secondary_position_y_km  DOUBLE,
    binary_secondary_position_z_km  DOUBLE,
    -- v18: secondary_mass / (primary_mass + secondary_mass), constant --
    -- stored so advance_orbital_phases never needs to join back to `stars`.
    binary_secondary_mass_fraction  DOUBLE,
    -- v18: additional pair-wide wobble from circumbinary (P-type) planets
    -- (star_id IS NULL) -- NULL unless binary_configuration = 'close'.
    binary_planetary_wobble_x_km    DOUBLE,
    binary_planetary_wobble_y_km    DOUBLE,
    binary_planetary_wobble_z_km    DOUBLE,

    system_flavor_text   TEXT,
    schema_version       INT NOT NULL DEFAULT 1,

    -- Fully rendered page text, ready to paste/upload. Both must be
    -- rendered from the SAME already-generated StarSystem object
    -- (toggle system_config.MARKDOWN, render, toggle back, render again)
    -- -- generation mixes unseedable `secrets` with seedable `random`
    -- (see spaceSector.py's module docstring), so the other format can
    -- never be faithfully regenerated later. LONGTEXT, not TEXT -- see
    -- the MySQL port note above.
    wikitext_content     LONGTEXT,
    markdown_content     LONGTEXT,

    -- One system = one wiki page (stars/planets/moons are sections within
    -- it, per StarSystem.__str__), so exactly one URL per wiki target.
    mediawiki_url        VARCHAR(2048),
    wikijs_url           VARCHAR(2048),

    created_at           TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT fk_star_systems_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE SET NULL,
    CONSTRAINT fk_star_systems_config
        FOREIGN KEY (system_config_id) REFERENCES system_configs(id),
    KEY idx_star_systems_sector_id (sector_id),
    KEY idx_star_systems_system_config_id (system_config_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- stars -- one row per single star, two rows (primary/secondary) per
-- binary. Never a row for the BinaryStarProxy itself (its snapshot lives
-- on star_systems above). Fields verified against starData.py:399-445 /
-- Star.__init__.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS stars (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_system_id            BIGINT UNSIGNED NOT NULL,
    role                      VARCHAR(16) NOT NULL CHECK (role IN ('primary', 'secondary', 'single')),
    name                      VARCHAR(255) NOT NULL,
    star_type                 VARCHAR(64) NOT NULL,   -- Star.type, e.g. "G2V Yellow Main Sequence Star"
    yerkes_class              VARCHAR(16) NOT NULL,
    mass_kg                   DOUBLE NOT NULL,
    radius_km                 DOUBLE NOT NULL,
    temperature_k             DOUBLE NOT NULL,
    luminosity_w              DOUBLE NOT NULL,
    age_gy                    DOUBLE NOT NULL,
    lifespan_gy               DOUBLE,            -- NULL = float('inf'), white dwarfs
    habitable_zone_inner_km   DOUBLE NOT NULL,
    habitable_zone_outer_km   DOUBLE NOT NULL,
    system_perimeter_km       DOUBLE NOT NULL,
    heliosphere_radius_km     DOUBLE NOT NULL,
    galactic_orbital_speed_kms    DOUBLE NOT NULL,
    galactic_orbital_period_gy    DOUBLE NOT NULL,
    galactic_orbital_phase_deg          DOUBLE NOT NULL,  -- v13, see header comment
    galactic_min_update_interval_years  DOUBLE NOT NULL,  -- v13, see header comment
    -- v15: this star's own Holman & Wiegert (1999) critical semi-major axis
    -- -- NULL except for a constituent of a 'wide' binary. See the header
    -- comment's "v15" note.
    wide_binary_a_crit_km               DOUBLE,
    -- v18 (see header comment): this star's own reflex-offset "wobble"
    -- from the combined pull of every planet orbiting it directly --
    -- NULL/0 with no planets. Class-blind: applies identically to an
    -- anchored black_holes/neutron_stars row.
    reflex_offset_x_km       DOUBLE,
    reflex_offset_y_km       DOUBLE,
    reflex_offset_z_km       DOUBLE,

    CONSTRAINT fk_stars_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    KEY idx_stars_star_system_id (star_system_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- planets -- top-level planets only as of schema v2 (see the header
-- comment's "v2" note); moons live in the `moons` table below instead of
-- self-referencing here. Fields verified against planetData.py:117-200 /
-- Planet.__init__ (27 scalar fields; excludes system_config/star
-- back-refs, the `moons` list [-> the moons table], `evolutionary_data`
-- [-> planet_evolutionary_paragraphs], and reflection_spectrum_visible/
-- non_visible [-> planet_reflection_spectrum]).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS planets (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_system_id            BIGINT UNSIGNED NOT NULL,
    -- The specific star this planet orbits, when that's a real stored
    -- `stars` row -- true for every single-star system and, as of v15,
    -- every 'wide' binary's planets too (each orbits one specific
    -- constituent star). Still NULL for a 'close' (P-type) binary's
    -- planets: systemData.py generates those against the merged
    -- `BinaryStarProxy`, never one individual constituent star, and the
    -- proxy is deliberately not stored as its own `stars` row (see
    -- star_systems.binary_* above). `star_system_id` above is always the
    -- reliable owning link regardless of star_id.
    star_id                   BIGINT UNSIGNED,
    orbital_index             INT NOT NULL,
    body_type                 VARCHAR(4) NOT NULL CHECK (body_type IN ('t', 'g')),
    name                      VARCHAR(255) NOT NULL,
    planet_class              VARCHAR(16),
    distance_km               DOUBLE NOT NULL,   -- from the star
    radius_km                 DOUBLE NOT NULL,
    mass_kg                   DOUBLE NOT NULL,
    volume_km3                DOUBLE NOT NULL,   -- derived, stored not recomputed
    period_years              DOUBLE NOT NULL,   -- derived, stored not recomputed
    zone                      VARCHAR(4) CHECK (zone IN ('h', 'e', 'c')),
    description               TEXT,
    gravity_g                 DOUBLE,
    surface_temperature_k     DOUBLE,
    density_g_cm3             DOUBLE,
    atmosphere                TEXT,
    atm_density               DOUBLE,
    atm_molar_density         DOUBLE,
    atmospheric_pressure_pa   DOUBLE,
    composition               TEXT,            -- descriptive string (contrast asteroid_belt_composition)
    scale_height_km           DOUBLE,
    hill_radius_km            DOUBLE,
    min_orbit_distance_km     DOUBLE,
    habitable_zone_inner_km   DOUBLE NOT NULL,   -- copied from host star at generation time
    habitable_zone_outer_km   DOUBLE NOT NULL,
    life_chemical             VARCHAR(64),
    evolutionary_speed        VARCHAR(64),
    flavor_text               TEXT,
    flavor_text_count         INT NOT NULL DEFAULT 0,
    -- Orbital motion (v9, see header comment) -- inclination/ascending
    -- node are fixed at generation time; phase changes over time, advanced
    -- in place by updateOrbits.py. rotation_period_hours is a separate,
    -- static "day length" stat.
    orbital_inclination_deg     DOUBLE NOT NULL,
    orbital_ascending_node_deg  DOUBLE NOT NULL,
    orbital_phase_deg           DOUBLE NOT NULL,
    -- Position/speed (v11, see header comment) -- Cartesian position
    -- relative to this planet's orbital anchor (the star, or the
    -- BinaryStarProxy's combined center for a binary system), derived from
    -- distance_km and the orbital-motion columns above; changes in
    -- lockstep with orbital_phase_deg as updateOrbits.py advances it.
    -- orbital_speed_kms is constant around a circular orbit -- it only
    -- changes if distance_km/period_years do.
    position_x_km               DOUBLE NOT NULL,
    position_y_km               DOUBLE NOT NULL,
    position_z_km               DOUBLE NOT NULL,
    orbital_speed_kms           DOUBLE NOT NULL,
    min_update_interval_years   DOUBLE NOT NULL,  -- v12, see header comment
    rotation_period_hours       DOUBLE NOT NULL,
    -- v18 (see header comment): this planet's own reflex-offset "wobble"
    -- from the combined pull of its own moons -- NULL/0 with no moons.
    reflex_offset_x_km       DOUBLE,
    reflex_offset_y_km       DOUBLE,
    reflex_offset_z_km       DOUBLE,

    CONSTRAINT fk_planets_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_planets_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,
    KEY idx_planets_star_system_id (star_system_id),
    KEY idx_planets_star_id (star_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Child table for the variable-length evolutionary narrative list
-- (Planet.evolutionary_data).
CREATE TABLE IF NOT EXISTS planet_evolutionary_paragraphs (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    planet_id   BIGINT UNSIGNED NOT NULL,
    position    INT NOT NULL,
    paragraph   TEXT NOT NULL,

    CONSTRAINT fk_planet_evolutionary_paragraphs_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    KEY idx_planet_evolutionary_paragraphs_planet_id (planet_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Child table replacing reflection_spectrum_visible/non_visible -- no
-- JSON column, same position+value shape as the paragraphs/composition
-- child tables (this schema is a pure relational design, not a
-- JSON-blob-in-a-column hybrid).
CREATE TABLE IF NOT EXISTS planet_reflection_spectrum (
    id             BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    planet_id      BIGINT UNSIGNED NOT NULL,
    spectrum_type  VARCHAR(16) NOT NULL CHECK (spectrum_type IN ('visible', 'non_visible')),
    position       INT NOT NULL,
    value          VARCHAR(255) NOT NULL,

    CONSTRAINT fk_planet_reflection_spectrum_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    KEY idx_planet_reflection_spectrum_planet_id (planet_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- moons -- one row per moon, as of schema v2 (see the header comment's
-- "v2" note). Exactly the same column shape as `planets` (moons are
-- `Planet` instances too, with `is_moon=True`), except `planet_id`
-- replaces `star_system_id`/`star_id`'s role as "what this body directly
-- orbits" -- `star_system_id`/`star_id` are still carried here too
-- (redundant with the owning planet's own columns) purely so a moon can
-- be queried/joined to its system or star without an extra hop through
-- `planets`. No self-reference: moons never generate their own moons
-- (`Planet.__init__` only calls `generate_moons` `if not self.is_moon`).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS moons (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    planet_id                 BIGINT UNSIGNED NOT NULL,
    star_system_id            BIGINT UNSIGNED NOT NULL,
    star_id                   BIGINT UNSIGNED,
    orbital_index             INT NOT NULL,   -- position in the parent planet's moons list
    body_type                 VARCHAR(4) NOT NULL CHECK (body_type IN ('t', 'g')),
    name                      VARCHAR(255) NOT NULL,
    planet_class              VARCHAR(16),
    distance_km               DOUBLE NOT NULL,   -- from the parent planet
    radius_km                 DOUBLE NOT NULL,
    mass_kg                   DOUBLE NOT NULL,
    volume_km3                DOUBLE NOT NULL,
    period_years              DOUBLE NOT NULL,
    zone                      VARCHAR(4) CHECK (zone IN ('h', 'e', 'c')),
    description               TEXT,
    gravity_g                 DOUBLE,
    surface_temperature_k     DOUBLE,
    density_g_cm3             DOUBLE,
    atmosphere                TEXT,
    atm_density               DOUBLE,
    atm_molar_density         DOUBLE,
    atmospheric_pressure_pa   DOUBLE,
    composition               TEXT,
    scale_height_km           DOUBLE,
    hill_radius_km            DOUBLE,
    min_orbit_distance_km     DOUBLE,
    habitable_zone_inner_km   DOUBLE NOT NULL,
    habitable_zone_outer_km   DOUBLE NOT NULL,
    life_chemical             VARCHAR(64),
    evolutionary_speed        VARCHAR(64),
    flavor_text               TEXT,
    flavor_text_count         INT NOT NULL DEFAULT 0,
    -- Orbital motion (v9, see header comment) -- inclination/ascending
    -- node are fixed at generation time; phase changes over time, advanced
    -- in place by updateOrbits.py. rotation_period_hours is a separate,
    -- static "day length" stat.
    orbital_inclination_deg     DOUBLE NOT NULL,
    orbital_ascending_node_deg  DOUBLE NOT NULL,
    orbital_phase_deg           DOUBLE NOT NULL,
    -- Position/speed (v11, see header comment) -- Cartesian position
    -- relative to this moon's orbital anchor, its parent planet; derived
    -- from distance_km and the orbital-motion columns above the same way
    -- planets.position_x/y/z_km are (see that table's comment).
    position_x_km               DOUBLE NOT NULL,
    position_y_km               DOUBLE NOT NULL,
    position_z_km               DOUBLE NOT NULL,
    orbital_speed_kms           DOUBLE NOT NULL,
    min_update_interval_years   DOUBLE NOT NULL,  -- v12, see header comment
    rotation_period_hours       DOUBLE NOT NULL,

    CONSTRAINT fk_moons_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    CONSTRAINT fk_moons_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_moons_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,
    KEY idx_moons_planet_id (planet_id),
    KEY idx_moons_star_system_id (star_system_id),
    KEY idx_moons_star_id (star_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Moon counterpart of planet_evolutionary_paragraphs.
CREATE TABLE IF NOT EXISTS moon_evolutionary_paragraphs (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    moon_id     BIGINT UNSIGNED NOT NULL,
    position    INT NOT NULL,
    paragraph   TEXT NOT NULL,

    CONSTRAINT fk_moon_evolutionary_paragraphs_moon
        FOREIGN KEY (moon_id) REFERENCES moons(id) ON DELETE CASCADE,
    KEY idx_moon_evolutionary_paragraphs_moon_id (moon_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Moon counterpart of planet_reflection_spectrum.
CREATE TABLE IF NOT EXISTS moon_reflection_spectrum (
    id             BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    moon_id        BIGINT UNSIGNED NOT NULL,
    spectrum_type  VARCHAR(16) NOT NULL CHECK (spectrum_type IN ('visible', 'non_visible')),
    position       INT NOT NULL,
    value          VARCHAR(255) NOT NULL,

    CONSTRAINT fk_moon_reflection_spectrum_moon
        FOREIGN KEY (moon_id) REFERENCES moons(id) ON DELETE CASCADE,
    KEY idx_moon_reflection_spectrum_moon_id (moon_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- asteroid_belts -- belts have no properties-dict data table
-- (asteroidData.py:93-141 is prose only), so their searchable columns are
-- the facts that prose always states instead: density, the distance range
-- from the star, and composition.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS asteroid_belts (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_system_id       BIGINT UNSIGNED NOT NULL,
    -- v15: same "which specific star this orbits" semantics as
    -- planets.star_id above -- NULL for a single star's or a 'close'
    -- binary's belt, set to the owning star for a 'wide' binary's.
    star_id              BIGINT UNSIGNED,
    orbital_index        INT NOT NULL,
    distance_km          DOUBLE NOT NULL,
    lower_limit_km       DOUBLE NOT NULL,   -- distance, low (asteroidData.py's "distance_text" range)
    upper_limit_km       DOUBLE NOT NULL,   -- distance, high
    density              VARCHAR(16) NOT NULL CHECK (density IN ('dense', 'sparse', 'typical')),

    -- Searchable composition summary, built the same way as the prose
    -- sentence in AsteroidBelt.to_paragraph_list() (asteroidData.py:119-137),
    -- e.g. "high concentrations of iron, moderate concentrations of nickel,
    -- and trace amounts of platinum" -- a single column so a belt's
    -- composition can be searched without a join, alongside the
    -- structured per-component breakdown in asteroid_belt_composition
    -- below for queries that need one specific component.
    composition_summary  TEXT NOT NULL,

    CONSTRAINT fk_asteroid_belts_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_asteroid_belts_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,
    KEY idx_asteroid_belts_star_system_id (star_system_id),
    KEY idx_asteroid_belts_star_id (star_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Structured per-component detail behind composition_summary above --
-- the belt's (component, concentration) list.
CREATE TABLE IF NOT EXISTS asteroid_belt_composition (
    id             BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    belt_id        BIGINT UNSIGNED NOT NULL,
    position       INT NOT NULL,
    component      VARCHAR(64) NOT NULL,
    concentration  VARCHAR(16) NOT NULL CHECK (concentration IN ('high', 'moderate', 'small', 'trace')),

    CONSTRAINT fk_asteroid_belt_composition_belt
        FOREIGN KEY (belt_id) REFERENCES asteroid_belts(id) ON DELETE CASCADE,
    KEY idx_asteroid_belt_composition_belt_id (belt_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- black_holes / neutron_stars -- v16 exotic phenomena, see this file's
-- header comment's "v16" note. Satellite tables extending a `stars` row
-- (nullable `star_id`) when a compact remnant anchors a full StarSystem
-- (`phenomenonGen.py --anchor-system`); standalone (`star_id` NULL)
-- otherwise.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS black_holes (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_id                   BIGINT UNSIGNED,
    name                      VARCHAR(255) NOT NULL,
    mass_solar                DOUBLE NOT NULL,
    event_horizon_radius_km   DOUBLE NOT NULL,
    spin                      DOUBLE NOT NULL,
    has_accretion_disk        TINYINT(1) NOT NULL CHECK (has_accretion_disk IN (0, 1)),
    temperature_k             DOUBLE NOT NULL,
    luminosity_w              DOUBLE NOT NULL,
    age_gy                    DOUBLE NOT NULL,
    -- v17: NULL when star_id is set (an anchored remnant's galactic motion
    -- lives on its own stars row instead); populated only for a standalone
    -- black hole. See this file's "v17" header note.
    galactic_orbital_speed_kms           DOUBLE,
    galactic_orbital_period_gy           DOUBLE,
    galactic_orbital_phase_deg           DOUBLE,
    galactic_min_update_interval_years   DOUBLE,

    CONSTRAINT fk_black_holes_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE CASCADE,
    KEY idx_black_holes_star_id (star_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS neutron_stars (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_id                   BIGINT UNSIGNED,
    name                      VARCHAR(255) NOT NULL,
    mass_solar                DOUBLE NOT NULL,
    radius_km                 DOUBLE NOT NULL,
    spin_period_ms            DOUBLE NOT NULL,
    magnetic_field_gauss      DOUBLE NOT NULL,
    pulsar_type               VARCHAR(16) NOT NULL CHECK (pulsar_type IN ('young', 'millisecond', 'non-pulsing')),
    surface_temperature_k     DOUBLE NOT NULL,
    luminosity_w              DOUBLE NOT NULL,
    age_gy                    DOUBLE NOT NULL,
    -- v17: same nullable-when-anchored convention as black_holes above.
    galactic_orbital_speed_kms           DOUBLE,
    galactic_orbital_period_gy           DOUBLE,
    galactic_orbital_phase_deg           DOUBLE,
    galactic_min_update_interval_years   DOUBLE,

    CONSTRAINT fk_neutron_stars_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE CASCADE,
    KEY idx_neutron_stars_star_id (star_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- nebulae -- v16 exotic phenomenon, always standalone. `sector_id` is
-- reserved for a future sector-context encounter; unused (always NULL)
-- by phenomenonGen.py today.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS nebulae (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id         BIGINT UNSIGNED,
    name              VARCHAR(255) NOT NULL,
    nebula_type       VARCHAR(16) NOT NULL CHECK (nebula_type IN ('emission', 'reflection', 'planetary', 'dark')),
    radius_ly         DOUBLE NOT NULL,
    composition       TEXT NOT NULL,
    formation_cause   TEXT NOT NULL,
    -- v17: always populated (a nebula is always standalone). See this
    -- file's "v17" header note.
    galactic_orbital_speed_kms           DOUBLE NOT NULL,
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,

    CONSTRAINT fk_nebulae_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    KEY idx_nebulae_sector_id (sector_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- supernova_remnants -- v16 exotic phenomenon, always standalone. At most
-- one of compact_remnant_black_hole_id/compact_remnant_neutron_star_id is
-- ever set (see compact_remnant_kind), and only for a core-collapse
-- progenitor whose collapsed core is still detectable -- never for a
-- Type Ia progenitor, which leaves nothing behind.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS supernova_remnants (
    id                                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id                         BIGINT UNSIGNED,
    name                              VARCHAR(255) NOT NULL,
    morphology                        VARCHAR(16) NOT NULL CHECK (morphology IN ('shell', 'plerion', 'composite')),
    age_years                         DOUBLE NOT NULL,
    radius_ly                         DOUBLE NOT NULL,
    progenitor_type                   VARCHAR(16) NOT NULL CHECK (progenitor_type IN ('Type Ia', 'core-collapse')),
    compact_remnant_kind              VARCHAR(16) CHECK (compact_remnant_kind IN ('black_hole', 'neutron_star')),
    compact_remnant_black_hole_id     BIGINT UNSIGNED,
    compact_remnant_neutron_star_id   BIGINT UNSIGNED,
    -- v17: always populated (a supernova remnant is always standalone).
    galactic_orbital_speed_kms           DOUBLE NOT NULL,
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,

    CONSTRAINT fk_supernova_remnants_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    CONSTRAINT fk_supernova_remnants_black_hole
        FOREIGN KEY (compact_remnant_black_hole_id) REFERENCES black_holes(id) ON DELETE SET NULL,
    CONSTRAINT fk_supernova_remnants_neutron_star
        FOREIGN KEY (compact_remnant_neutron_star_id) REFERENCES neutron_stars(id) ON DELETE SET NULL,
    KEY idx_supernova_remnants_sector_id (sector_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- rogue_planets -- v16 exotic phenomenon, always standalone (a rogue
-- planet is by definition unbound from any star).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS rogue_planets (
    id                  BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id           BIGINT UNSIGNED,
    name                VARCHAR(255) NOT NULL,
    planet_type         VARCHAR(4) NOT NULL CHECK (planet_type IN ('t', 'g')),
    mass_kg             DOUBLE NOT NULL,
    radius_km           DOUBLE NOT NULL,
    composition         TEXT NOT NULL,
    has_internal_heat   TINYINT(1) NOT NULL CHECK (has_internal_heat IN (0, 1)),
    has_moons           TINYINT(1) NOT NULL CHECK (has_moons IN (0, 1)),
    -- v17: always populated (a rogue planet is always standalone).
    galactic_orbital_speed_kms           DOUBLE NOT NULL,
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,

    CONSTRAINT fk_rogue_planets_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    KEY idx_rogue_planets_sector_id (sector_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- interstellar_comets -- v16 exotic phenomenon, always standalone (an
-- interstellar comet is, by definition, unbound from any star).
-- `interstellar_comet_composition` mirrors `asteroid_belt_composition`'s
-- per-component breakdown, minus a concentration level (InterstellarComet's
-- own composition list carries no such gradient).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS interstellar_comets (
    id                     BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id              BIGINT UNSIGNED,
    name                   VARCHAR(255) NOT NULL,
    nucleus_diameter_km    DOUBLE NOT NULL,
    velocity_kms           DOUBLE NOT NULL,
    is_active              TINYINT(1) NOT NULL CHECK (is_active IN (0, 1)),
    composition_summary    TEXT NOT NULL,
    -- v17: always populated (an interstellar comet is always standalone).
    -- Its own hyperbolic velocity_kms above is a separate, non-advancing
    -- descriptive stat -- see this file's "v17" header note.
    galactic_orbital_speed_kms           DOUBLE NOT NULL,
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,

    CONSTRAINT fk_interstellar_comets_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    KEY idx_interstellar_comets_sector_id (sector_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS interstellar_comet_composition (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    comet_id    BIGINT UNSIGNED NOT NULL,
    position    INT NOT NULL,
    component   VARCHAR(64) NOT NULL,

    CONSTRAINT fk_interstellar_comet_composition_comet
        FOREIGN KEY (comet_id) REFERENCES interstellar_comets(id) ON DELETE CASCADE,
    KEY idx_interstellar_comet_composition_comet_id (comet_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- asteroid_fields -- v17 exotic phenomenon (seventh phenomenon type),
-- always standalone (a field drifting in open space, as opposed to
-- asteroid_belts, which always orbits a star). `sector_id` is reserved
-- for a future sector-context encounter; unused (always NULL) by
-- phenomenonGen.py today, the same convention `nebulae` uses.
-- `asteroid_field_composition` mirrors `asteroid_belt_composition`
-- exactly (same per-component/concentration shape -- both are generated
-- via the same shared `asteroidData.generate_asteroid_composition`).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS asteroid_fields (
    id                    BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id             BIGINT UNSIGNED,
    name                  VARCHAR(255) NOT NULL,
    density               VARCHAR(16) NOT NULL CHECK (density IN ('dense', 'sparse', 'typical')),
    radius_ly             DOUBLE NOT NULL,
    composition_summary   TEXT NOT NULL,
    galactic_orbital_speed_kms           DOUBLE NOT NULL,
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,

    CONSTRAINT fk_asteroid_fields_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    KEY idx_asteroid_fields_sector_id (sector_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS asteroid_field_composition (
    id             BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    field_id       BIGINT UNSIGNED NOT NULL,
    position       INT NOT NULL,
    component      VARCHAR(64) NOT NULL,
    concentration  VARCHAR(16) NOT NULL CHECK (concentration IN ('high', 'moderate', 'small', 'trace')),

    CONSTRAINT fk_asteroid_field_composition_field
        FOREIGN KEY (field_id) REFERENCES asteroid_fields(id) ON DELETE CASCADE,
    KEY idx_asteroid_field_composition_field_id (field_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- sector_objects -- unified "every stellar object in a sector" search.
--
-- Standard relational choice here: a VIEW that UNIONs the four typed
-- tables above, not a fifth physical table duplicating their rows. A
-- real table would need to be kept in sync on every insert/update/delete
-- to the tables it mirrors (or drift out of sync); a view has no storage
-- and no sync problem -- MySQL resolves it against current data on every
-- query, and each row still traces back to its real table via
-- (object_type, object_id). Query with e.g.
-- `SELECT * FROM sector_objects WHERE sector_id = %s`.
-- ---------------------------------------------------------------------
CREATE OR REPLACE VIEW sector_objects AS
    SELECT
        'star'                          AS object_type,
        s.id                            AS object_id,
        s.star_system_id                AS star_system_id,
        ss.sector_id                    AS sector_id,
        s.name                          AS name,
        s.star_type                     AS summary,
        NULL                            AS orbital_index
    FROM stars s
    JOIN star_systems ss ON ss.id = s.star_system_id

    UNION ALL

    SELECT
        'planet'                        AS object_type,
        p.id                            AS object_id,
        p.star_system_id                AS star_system_id,
        ss.sector_id                    AS sector_id,
        p.name                          AS name,
        COALESCE(p.planet_class, p.body_type) AS summary,
        p.orbital_index                 AS orbital_index
    FROM planets p
    JOIN star_systems ss ON ss.id = p.star_system_id

    UNION ALL

    SELECT
        'moon'                          AS object_type,
        m.id                            AS object_id,
        m.star_system_id                AS star_system_id,
        ss.sector_id                    AS sector_id,
        m.name                          AS name,
        COALESCE(m.planet_class, m.body_type) AS summary,
        m.orbital_index                 AS orbital_index
    FROM moons m
    JOIN star_systems ss ON ss.id = m.star_system_id

    UNION ALL

    SELECT
        'asteroid_belt'                 AS object_type,
        ab.id                           AS object_id,
        ab.star_system_id               AS star_system_id,
        ss.sector_id                    AS sector_id,
        'Asteroid Belt'                 AS name,
        ab.density                      AS summary,
        ab.orbital_index                AS orbital_index
    FROM asteroid_belts ab
    JOIN star_systems ss ON ss.id = ab.star_system_id;
