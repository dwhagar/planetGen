-- planetgen/db/schema.sql
--
-- MySQL (InnoDB) schema for planetGen persistence (TODO.md Phase 2, ported
-- off SQLite for Phase 5's MySQL migration -- see docs/TODO.md). Plain SQL
-- DDL, no ORM. Column lists are verified against the actual generator
-- source (planetgen/generation/config.py, starData.py, doubleStar.py,
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
--     place km's numbers get unwieldy (a 4-pc sector edge is
--     ~1.23e14 km vs. 4000 mpc) and because sector geometry has no other
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
-- planetgen/physics/units.py: `ly_to_milliparsecs`/`milliparsecs_to_ly` for
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
-- and planetgen/web/maps/galaxymap.py introduced later (4 azimuthal regions
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
-- never placed in a sector. Computed in `planetgen/db/store.py` from
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
-- `store.py`.
--
-- v4: added six nullable galaxy-frame placement columns to `sectors` --
-- center_x/y/z_pc (Cartesian center, parsecs, galactic-origin-relative),
-- galactic_radius_pc (sqrt(x^2+y^2+z^2), persisted for the same
-- can't-index-an-expression reason star_systems.quadrant is persisted
-- rather than recomputed), and shell_index/shell_slot_index (the
-- sector's stable address within the shell/Fibonacci-sphere tiling
-- scheme -- see docs/design/galaxy/galaxy-coordinate-system.md). All six are
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
-- self-reference of its own. `planetgen/db/store.py`'s `migrate_database`
-- converts an existing v1 database in place (backing up the original
-- first); `planetgen.cli.migrate` runs it over every database in a directory, and
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
-- pre-rendered copy (`planetgen/web/lib/tabledisplay.py` for the HTML viewer;
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
--     handful of numbers `planetgen.galaxy.density`/`galaxySkeleton`
--     need to recompute any sector's exact position/density on demand.
--   - `galaxy_shell_band`: one row per contiguous *candidate* slot-index
--     band per shell (`planetgen.galaxy.skeleton.find_shell_bands`) --
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
--   track rotational phase, only period). `planetgen.cli.orbits` is a new,
--   separately-run script that advances every body's `orbital_phase_deg`
--   in place based on `period_years` and real elapsed time, using the new
--   `orbit_simulation_state` singleton row (one per database, same
--   pattern as `galaxy_shape` above) to track when it last ran --
--   `planetgen.galaxy.geometry`'s "derive, don't store" principle
--   doesn't apply here: the whole point of this feature is a value that
--   changes with real-world time rather than being a pure function of the
--   body's other (fixed) generative facts, so it has to be stored and
--   periodically updated, not recomputed on read.
--
-- v10: galactic orbit. `stars` gains `galactic_orbital_speed_kms` and
--   `galactic_orbital_period_gy` -- a star system's circular orbital speed
--   and period around the galactic center, derived from its (or the fallback's, the radius that holds half of the galaxy's disc stars) distance
--   from it via a simple rotation-curve model (see
--   `physical_constants.GALACTIC_ROTATION_FLAT_VELOCITY_KMS`'s comment and
--   `galactic_orbit.calculate_galactic_orbit`) -- the same "fixed at generation
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
--   `orbital_phase_deg`) via `orbits.orbital_position_au` -- the same
--   "each body positioned relative to its immediate primary" convention
--   `docs/design/galaxy/galaxy-coordinate-system.md` already uses one level up
--   for sectors/systems relative to the galactic center (see v10 above).
--   `planetgen.cli.orbits`/`_db.advance_orbital_phases` recomputes position in
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
--   from `period_years` (`orbits.minimum_update_interval_years`:
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
--   guard v12 added for planets/moons, `orbits.minimum_update_interval_years`
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
--   `orbits.circular_orbital_speed_kms`) applied to the pair's already-
--   stored `binary_separation_km`/`binary_effective_mass_kg` -- the same
--   formulas a planet's orbit around its star already uses, just with the
--   combined pair mass standing in for "the primary". `_inclination_deg`/
--   `_ascending_node_deg`/`_phase_deg` fully orient this mutual orbit in
--   3D and track the pair's current position within it (the "orbital
--   direction" a binary's own orbital plane needs, unlike the galactic
--   orbit above, which this generator treats as planar) -- the same
--   `orbits.orbital_position_au` orbital-element convention planets/moons
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
--   `orbits.orbital_position_au`, the same "each body positioned relative
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
--   real and used for `orbits.holman_wiegert_critical_semimajor_axis`, but
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
--   sector-context encounter (unused by `phenomenonGen.py` at the time,
--   which never created or attached a sector -- see this file's "v21"
--   header note for `sectorGen.py`'s own later per-sector generation,
--   which does). `supernova_remnants` additionally
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
-- v18: galaxy-frame placement for nebulae/asteroid fields, the first real
--   use of the `sector_id` column v16/v17 left "reserved for a future
--   sector-context encounter" on both tables. A nebula/asteroid field is
--   frequently far larger than a single sector's cube (an emission nebula
--   can span up to 200 ly; the default sector edge is 4 pc, ~13 ly), so unlike
--   `star_systems.position_x/y/z_mpc` (relative to one owning sector's own
--   center) these are placed directly in the same galaxy-frame Cartesian
--   space `sectors.center_x/y/z_pc` already uses (`docs/design/galaxy/
--   galaxy-coordinate-system.md`) -- a sphere (`center_x/y/z_pc`,
--   `radius_ly`) that may overlap zero, one, or several sectors' cubes,
--   not a single sector-relative offset. `nebulae`/`asteroid_fields` each
--   gain `center_x_pc`/`center_y_pc`/`center_z_pc`/`galactic_radius_pc` --
--   NULL together, mirroring `sectors`' own v4 null-together convention --
--   plus `sector_id` is now actually populated (`phenomenonGen.py
--   --sector-id`) with the *nearest* already-generated sector to the
--   phenomenon's own center, a convenience "home" link for browsing, not
--   the authoritative geometry (that's the sphere itself -- see
--   `queryDb.phenomena_near_sector`, which finds every phenomenon whose
--   sphere overlaps a given sector's cube by real distance, regardless of
--   which sector it's linked to).
--
-- v19: star-bound comets (`cometData.Comet` -- see
--   docs/design/systems/comet-orbital-realism.md), propagated via real two-body
--   Kepler/Barker orbital mechanics (`keplerMotion.py`) rather than a
--   fixed hyperbolic excess speed. New tables `comets` (plus child
--   `comet_composition`, mirroring `interstellar_comet_composition`'s
--   plain-component-list shape, not `asteroid_belt_composition`'s
--   component+concentration one) and a `comet` branch on `sector_objects`
--   below. Deliberately NOT part of `planets`/`asteroid_belts` (no
--   `orbital_index`): a comet's `distance_km` is its current,
--   continuously-varying position along an eccentric/parabolic orbit, not
--   a fixed slot in the orbital-spacing sequence those tables share --
--   see `systemData.StarSystem._generate_comets`'s docstring. `star_id`
--   follows `planets.star_id`'s own real semantics (a real `stars.id` for
--   a single star or a 'wide' binary's comet; NULL only for a 'close'
--   binary's, which orbits the merged pair, not one stored star row).
--   New tables, so (like v16/v17's new tables) no `ALTER TABLE` is needed
--   in `_migrate_v18_to_v19`.
--
-- v20: proper two-body (barycentric) trajectories for binary stars,
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
--   nominal fixed point instead -- see `orbits.calculate_reflex_offset`'s
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
--   `stars`, `planets`), so `_migrate_v19_to_v20` needs real `ALTER TABLE`
--   statements, the same as v17's/v18's column additions.
--
-- v21: sector-level exotic phenomena (`sectorGen.generate_sector_phenomena`
--   -- every generated sector now also seeds a realistically sparse
--   population of phenomenonGen.py's own seven types, not just
--   phenomenonGen.py's separate on-demand CLI). `black_holes`/
--   `neutron_stars` gain the same nullable `sector_id`/`center_x/y/z_pc`/
--   `galactic_radius_pc` galaxy-frame placement shape v18 already gave
--   `nebulae`/`asteroid_fields` -- a standalone compact remnant is a real,
--   stellar-mass gravitating body, so (unlike every other phenomenon type)
--   its own in-sector position is the exact one
--   `SpaceSector.add_phenomenon`'s Hill-sphere-aware placement computed at
--   generation time (never within a neighboring star system's or another
--   compact remnant's own Hill sphere), converted to an absolute
--   galaxy-frame center via `_db._galaxy_placement_from_sector_offset`
--   rather than `compute_phenomenon_placement`'s independent random
--   jitter (which still applies as before for `phenomenonGen.py --type
--   black-hole/neutron-star --sector-id`'s own standalone use, where no
--   such specific in-sector position exists). `supernova_remnants`/
--   `rogue_planets`/`interstellar_comets` needed no schema change at all
--   -- their `sector_id` column already existed (v16, "reserved for a
--   future sector-context encounter"); `_db.save_phenomenon` now actually
--   populates it instead of leaving it perpetually NULL.
--
-- v22: search-facing indexes on `sectors`/`star_systems`/`stars`/`planets`/
--   `moons`.`name` (the `GET /api/search` autocomplete lists' own `SELECT
--   DISTINCT name ... ORDER BY name LIMIT n` -- an index lets this walk the
--   index in name order and stop at the limit instead of sorting the whole
--   table) and on `stars.yerkes_class`/`planets.planet_class`/
--   `planets.body_type`/`planets.life_class`/`moons.planet_class`/
--   `moons.body_type`/`moons.life_class`/`asteroid_belts.density` (the
--   search page's own facet counts -- `GROUP BY` on one of these -- and
--   result-panel `IN (...)` filters). None of these speed up a name *term*
--   search itself (`name LIKE '%text%'` -- a leading wildcard, which no
--   plain B-tree index can accelerate under any engine; that would need a
--   FULLTEXT index instead, a bigger change than this migration's own
--   scope), but every one of them was a genuine full-table scan/sort on
--   every single search-page visit regardless of whether any filter was
--   even active (`search()`'s facet+autocomplete computation always runs
--   up front) -- confirmed in production as the API timing out
--   (`urllib.error.URLError`/`TimeoutError` from `planetgen/web/lib/apiclient.py`)
--   on `/api/search` once the database grew past a trivial size.
--
-- v23: wiki publishing wired up (docs/TODO.md's "Wiki publishing isn't
--   wired up yet" item; see `src/wikiClient/`, `planetgen/api/routes.py`'s
--   `POST /api/systems/<id>/wiki`/`POST /api/sectors/<id>/wiki`, and
--   `html/system.py`/`html/sector.py`/`html/admin.py`). `sectors` gains
--   `wiki_url`, a single manually-set-or-uploaded link (a sector has no
--   generated page content of its own the way a system does -- see
--   `_sector_wiki_content` in `planetgen/api/routes.py` -- so one plain URL is
--   enough; there is no per-backend pair the way `star_systems` has,
--   since nothing here ever needs to know which wiki software served it,
--   only where to link). `star_systems.wikijs_url`/`mediawiki_url`
--   (present in this file's `CREATE TABLE` since before v23, but never
--   populated or migrated for an existing database -- see `store.py`'s
--   `_migrate_v22_to_v23`) are the per-system equivalent, one column per
--   backend since a system's own pre-rendered `wikitext_content`/
--   `markdown_content` can genuinely be uploaded to both a MediaWiki and
--   a Wiki.js instance independently. Existence on a given wiki is never
--   a separate stored flag -- it's exactly "this column is not NULL" --
--   `html/system.py`'s description section swaps to a same-tab-opens-a-
--   new-tab link whenever either is set, per the same convention
--   `html/sector.py`'s link swap uses for `wiki_url`.
--
-- v24: name-uniqueness registries (`planetgen/names/uniqueness.py`) --
--   `sector_name_registry`/`system_name_registry`/`body_name_registry`,
--   one per level of the sector > system > planet/moon naming hierarchy,
--   tracking each base name's progress through that level's own
--   collision-decoration scheme (Greek/Roman letters for sectors and
--   systems against their own kind, a diminutive prefix for a system
--   against a sector, a companion suffix for a planet/moon against
--   anything). `store.py`'s `insert_sector`/`insert_star_system`/
--   `insert_planet`/`insert_moon` all consult and update these now, so no
--   two rows anywhere in this database ever end up sharing a display
--   name. See each table's own comment below for the exact shape.
--
-- v25: a spatial index on `sectors.center_x_pc`/`center_y_pc`/
--   `center_z_pc` -- `queryDb.galaxy_sectors_in_view` (the interactive 3D
--   Galaxy Map's live-viewport `GET /api/galaxy/view`, called repeatedly
--   as its camera moves, plus once server-side on every `galaxy.py`
--   page load for the zoomed-all-the-way-out starting view) runs a
--   bounding-box `WHERE center_x_pc BETWEEN ... AND center_y_pc BETWEEN
--   ... AND center_z_pc BETWEEN ...` on every call; without an index
--   covering these columns that's a genuine full-table scan every time,
--   the exact same failure mode v22's note above already confirmed in
--   production for `/api/search` before that migration -- confirmed
--   again here (`TimeoutError`/"Truncated or oversized response headers"
--   from the WSGI daemon) once this feature shipped against a
--   non-trivial `sectors` table. A composite `(center_x_pc, center_y_pc,
--   center_z_pc)` index lets MySQL range-scan on the first column instead
--   of examining every row, same reasoning as v22's indexes.
--
-- v26: the same v25 spatial-index fix, applied to `nebulae`/
--   `asteroid_fields`/`black_holes`/`neutron_stars`. `queryDb.
--   phenomena_near_sector` -- called on every `GET /api/sectors/<id>`
--   (`html/sector.py`'s own page, and the Sector Map it renders) -- used
--   to read every galaxy-placed row from all four tables unconditionally
--   (`WHERE center_x_pc IS NOT NULL`, no bounding box at all) and filter
--   by distance in Python; confirmed in production as the same failure
--   mode v25's note documents (`TimeoutError`/"Truncated or oversized
--   response headers", severe enough that other, unrelated pages sharing
--   the single `planetgen-api` WSGIDaemonProcess -- see `docs/apache-
--   deployment.md`'s "MySQL accounts" note -- stalled behind it too,
--   until Apache was restarted). Now scoped to a SQL bounding box first
--   (`queryDb._placed_phenomenon_rows`'s own `bbox` argument), padded by
--   the widest `radius_ly` actually present in `nebulae`/`asteroid_fields`
--   at call time (queried directly, rather than assuming a fixed ceiling,
--   so this stays exactly as correct as the old unconditional scan for an
--   arbitrarily large placed phenomenon -- see that function's own
--   docstring) -- these four composite indexes are what let MySQL
--   range-scan that bounding box instead of still examining every row.
--
-- v27: row timestamps on the top-level tables -- `sectors`,
--   `star_systems` and the seven exotic-phenomenon tables (`store.py`'s
--   `TIMESTAMPED_TABLES`). Each gets `created_at` (`star_systems` already
--   had one) and `modified_at`, plus an index on `modified_at` so "what
--   changed since T" is a range scan rather than a full one -- the shape
--   a tile/page cache needs to decide what to rebuild.
--   - `modified_at` is `TIMESTAMP(3) ... ON UPDATE CURRENT_TIMESTAMP(3)`,
--     so MySQL itself bumps it on any `UPDATE` that changes the row --
--     no call site has to remember to. Millisecond precision so two
--     changes a fraction of a second apart still compare as different.
--   - Child rows (`stars`, `planets`, `moons`, `asteroid_belts`,
--     `comets` and their own detail tables) get no timestamps. A change
--     to one bumps its parent system's `modified_at` instead, through
--     `_db.touch_star_system` -- cheaper than a timestamp (and index) on
--     the largest tables, and a cache only ever needs to know that the
--     system changed. Deleting a system bumps its sector the same way
--     (`_db.touch_sector`). Generation writes a system and its children
--     together, so the system's own insert already covers them.
--   - `planetgen.cli.orbits`'s orbit ticks (`_db.advance_orbital_phases`) do
--     NOT count as a modification: each of its `UPDATE`s sets
--     `modified_at = modified_at`, which stops `ON UPDATE` from firing.
--     The simulation clock moving would otherwise mark every row changed
--     on every run; `orbit_simulation_state.last_updated_at` is the
--     timestamp for that instead.
--   - Done in application code rather than triggers, which on a
--     binary-logged server need privileges a web-app account usually
--     doesn't have, and couldn't tell an orbit tick from an edit anyway.
--   - `_migrate_v26_to_v27` adds the columns with `ALGORITHM=INSTANT`
--     where supported (a metadata-only change, no table rebuild) and the
--     indexes online (`ALGORITHM=INPLACE, LOCK=NONE`), so it's safe to
--     run against a large, live database. It then backfills in
--     primary-key batches (committing after each): a system's
--     `modified_at` starts at its own `created_at`, and a sector takes
--     its oldest system's `created_at` for both columns. Phenomena (and
--     sectors with no systems) have nothing to recover from, so they
--     keep the migration's own time.
--
-- v28: galaxy-frame placement for the last three phenomenon types --
--   `supernova_remnants`, `rogue_planets` and `interstellar_comets` gain
--   the same `center_x_pc`/`center_y_pc`/`center_z_pc`/
--   `galactic_radius_pc` columns (plus the null-together CHECK and the
--   v26-style spatial index) `nebulae`/`asteroid_fields` (v18) and
--   `black_holes`/`neutron_stars` (v21) already have. Before this, those
--   three were only ever linked to a sector by `sector_id`, so they never
--   appeared on the Sector Map or in `queryDb.phenomena_near_sector` even
--   though `sectorGen` placed each one at a real spot inside its sector.
--   `_db.insert_sector` now stores that spot, and a supernova remnant's
--   embedded black hole/neutron star gets the remnant's own sector and
--   center too.
--   - `_migrate_v27_to_v28` adds the columns with `ALGORITHM=INSTANT`
--     where supported and the indexes online, like v27. It then backfills
--     every existing row that has a galaxy-placed `sector_id`: the real
--     in-sector position was never saved, so each gets a random point
--     inside its own sector's cube (seeded by the row's id, so a re-run
--     picks the same point). Rows with no sector, or whose sector was
--     never placed, stay NULL.
--
-- v29: dropped `star_systems.wikitext_content`/`markdown_content`, the
--   full wiki page text rendered once at generation time. Every value that
--   text was built from already has its own column (the same reasoning v5
--   applied to the `table_*` snapshot columns), so both formats are now
--   rendered on demand: `_db.load_star_system` rebuilds the generation
--   object graph and `planetgen/db/render.py` renders it
--   (`StarSystem.__str__`, toggling `SystemConfig.MARKDOWN`). Rendered
--   fresh, a page follows every later change the stored copy never saw --
--   a rename, a name made unique after the text was rendered (planets and
--   moons get their suffix in `insert_planet`/`insert_moon`, after
--   `insert_star_system` used to render), and each orbit tick's current
--   wobble and comet positions. Checked before the drop by rendering every
--   system and comparing against the stored copy; the only fixes needed
--   were binary-pair orientation on load (`BinaryStarProxy.from_dict`/
--   `WideBinaryPair.from_dict` re-derive which star is the heavier,
--   matching generation). `_migrate_v28_to_v29` drops the two columns.
--
-- v30: no shape change -- a data cleanup of two generator bugs in
--   `planets`/`moons`. A body `planetPhysics.reconcile_zone_and_class`
--   moved into an airless class kept its old class's `atm_density`/
--   `atm_molar_density`/`scale_height_km` (now cleared on reclassification,
--   and NULLed on `atmosphere = 'None'` rows), and bodies around very dim
--   stars could come out below the cosmic microwave background (now
--   floored at `physical_constants.COSMIC_BACKGROUND_TEMPERATURE_K`,
--   2.725 K, and raised to it on existing rows). `_migrate_v29_to_v30` does
--   both updates.
-- v31: new `quasars` table -- a galaxy's active nucleus (`quasarData.Quasar`).
--   Only ever generated at the galactic center (0, 0, 0), by the one
--   core sector that hosts the nucleus (`generate.add_galactic_nucleus`;
--   ring 0, layer 0, slot 0 since v32), so it has the usual placement columns but no galactic-orbit columns:
--   it is the point everything else orbits. A brand-new table, so
--   `_ensure_schema`'s `CREATE TABLE IF NOT EXISTS` creates it;
--   `_migrate_v30_to_v31` only records the version.
--
-- v32: galaxy-placed sectors moved from spherical shells to a cylindrical
--   grid (see docs/design/galaxy/galaxy-coordinate-system.md, "Cylindrical
--   sector grid"). A sector's address is now `(ring_index, layer_index,
--   ring_slot_index)`: ring = cylindrical radius band, layer = height band
--   centered on the plane, slot = angular wedge, each about one sector
--   edge (11.5 ly) across. `sectors.shell_index`/`shell_slot_index` give
--   way to those three columns (UNIQUE together); `sector_vertices` is
--   dropped, since a cell's corners are closed-form
--   (`galaxyGeometry.sector_cell_vertices_pc`); `galaxy_shape.
--   outer_shell_index` becomes `outer_ring_index`; and `galaxy_shell_band`
--   becomes `galaxy_ring_band`, one row per ring holding the layer range
--   that can hold content. Star-system offsets
--   (`star_systems.position_x/y/z_mpc`) are now along the sector's
--   cylindrical axes (local +X radial, +Y along the ring, +Z galactic
--   north -- `galaxyGeometry.sector_orientation`). Shell-addressed sectors
--   can't be mapped onto the new cells, so `_migrate_v31_to_v32` deletes
--   every galaxy-placed sector together with its systems and phenomena
--   and rebuilds the skeleton from the stored shape; visiting the galaxy
--   regenerates them. Standalone (never placed) sectors are untouched.
--
-- v33: one sector standard (see docs/design/galaxy/galaxy-coordinate-system.md,
--   "Sector size" and "Layers"). The edge is a whole number of parsecs,
--   4 pc (~13.05 ly) by default, instead of 11.5 ly; ring `i` holds
--   `round(2*pi*(i + 1/2))` slots (3, 9, 16, 22, ...) instead of a
--   multiple of 4, identical on every layer so columns line up; and the
--   skeleton is stored per layer instead of per ring: `galaxy_ring_band`
--   becomes `galaxy_layer`, one row per layer from the highest to the
--   lowest holding the last ring that layer reaches, alongside
--   `galaxy_column`, one row per ring holding the lowest and highest
--   layer its column reaches. Every galaxy generation path checks an
--   address against this outline before generating anything. Almost every address
--   changes, so `_migrate_v32_to_v33` deletes every galaxy-placed sector
--   with its systems and phenomena (as v32 did) and rebuilds the skeleton
--   from the stored shape at the standard edge.
--
-- v34: planets and moons are named from their system (`bodyNames.py`):
--   `<system> I`, `<system> II`, moons `<system> IIa`, and a binary's stars
--   `<system> <word>`. The system name is already unique, so derived names
--   are too, and `body_name_registry` (v24's planet/moon registry) is
--   dropped. Existing rows keep their old names until regenerated.
--
-- v35: ring slot counts follow the hybrid master-wedge rule
--   (`galaxyGeometry.ring_sector_count`/`ring_master_count`, see
--   docs/design/galaxy/galaxy-coordinate-system.md): 3 master wedges at the
--   center, doubling once each would hold 8 slots, and each ring's count
--   the multiple of its master count nearest `2*pi*(i + 1/2)` (3, 9, 15,
--   21, 27, 36, ...), so slot boundaries line up on master lines from the
--   center out. No column changes, but stored `ring_slot_index` values
--   change meaning, so `_migrate_v34_to_v35` deletes every sector in a
--   ring whose count changed (all but 15 of 3,856) with its systems and
--   phenomena, as v32 and v33 did.
--   `galaxy_layer`/`galaxy_column` are keyed by ring and layer and stay.
--
-- v36: every galaxy gets a central supermassive black hole: when the
--   nucleus roll comes up quiescent (no quasar), generate.
--   add_galactic_nucleus places a `black_holes` row at the galactic
--   center. `black_holes.mass_class` ('stellar', 'intermediate',
--   'supermassive') records which kind a row is; `_migrate_v35_to_v36`
--   fills it from each existing row's mass. Intermediate-mass black holes
--   now span 1e2-1e5 Msun (log-uniform).
--
-- v37: research-based interstellar rates (docs/design/systems/
--   interstellar-object-rates.md). `rogue_planets.mass_bin` records the
--   mass bin a rogue was drawn from ('terrestrial', 'sub-neptune',
--   'saturn', 'jupiter') or 'brown-dwarf' (free-floating brown dwarfs
--   share the table); `star_systems.runaway_class`/`runaway_speed_kms`
--   flag runaway and hypervelocity stars. `_migrate_v36_to_v37` fills
--   `mass_bin` from each existing row's mass.
--
-- v38: letter classes (docs/design/systems/nebula-and-asteroid-field-classes.md).
--   `nebulae.nebula_class` (A-Q) and `supernova_remnants.remnant_class`
--   (R-W) come from program_constants.NEBULA_CLASSES, each with its
--   contents: `dominant_species`, `density_cm3` (nH), `temperature_k` and
--   `extinction_av` (magnitudes). `nebula_type` is now the class's family
--   and gains 'diffuse' (its CHECK is named chk_nebulae_type). Asteroid
--   fields get `field_class`, a letter from composition and density plus
--   a size digit (e.g. 'C3'), and `composition_family`.
--   `_migrate_v37_to_v38` infers classes for existing rows and fills
--   each class's typical contents.
--
-- v39: containment. star_systems, rogue_planets, interstellar_comets,
--   black_holes, neutron_stars, asteroid_fields and nebulae get
--   `inside_nebula_id`/`inside_remnant_id`: the innermost (smallest)
--   nebula or supernova remnant whose sphere holds the object, at most
--   one of the two set. `_db.refresh_containment` sets them by a 3D
--   distance test when a sector is generated, when a nebula or remnant is
--   placed (for every sector it reaches), and in `_migrate_v38_to_v39`.
--   A nebula only nests inside a larger cloud. ON DELETE SET NULL, so
--   deleting a cloud frees what was inside it.
--
-- v40: names that follow one standard. Black holes, neutron stars,
--   nebulae, supernova remnants, rogue planets and quasars draw their
--   names through system_name_registry like star systems do (a supernova
--   remnant's core is "<remnant> Core" and isn't registered separately),
--   so `first_star_system_id` is nullable and `first_object_table`/
--   `first_object_id` name a phenomenon holder. Comets and asteroid fields
--   get designations instead: `P/<host>-<n>` (periodic, under 200 years)
--   or `C/<host>-<n>` for a star-bound comet, `I/<sector>-<n>` for an
--   interstellar comet and `AF <field_class>-<sector>-<nn>` for an
--   asteroid field, where <sector> is the sector's grid designation
--   (`galaxyGeometry.provisional_sector_designation`) or its name. Each
--   of those phenomenon tables gains an index on `name`.
--   `_migrate_v39_to_v40` renames and registers existing rows.
--
-- v41: where phenomena sit, and what's near everything. Every placeable
--   phenomenon table (black_holes, neutron_stars, nebulae,
--   supernova_remnants, rogue_planets, interstellar_comets,
--   asteroid_fields, quasars) gains `quadrant`, the sector octant its
--   center sits in (the label `star_systems.quadrant` uses, in the same
--   `galaxyGeometry.sector_orientation` frame). New `nearest_systems`
--   table: the 3 nearest star systems to every placed system and
--   phenomenon, across sector boundaries. `_migrate_v40_to_v41` adds
--   both and fills them.
--
-- v42: facilities (starbases, colonies, outposts) -- the `facilities`
--   table below, one host each. `_migrate_v41_to_v42` creates it.
--
-- v43: bright-star pre-placement -- the `bright_stars` table below, and
--   `galaxy_shape.bright_star_min_luminosity_sol`/`bright_star_seed`.
--   `_migrate_v42_to_v43` adds them (empty: the next plan scatters).
--
-- v44: population and politics (POP.1 to POP.4,
--   docs/design/habitability/population-and-politics.md) -- the `species`, `polities`,
--   `system_owners` and `population_state` tables below, all filled by
--   `planetgen population` (planetgen/population/model.py).
--   `_migrate_v43_to_v44` creates them empty.
--
-- v45: id blocks for batched writes (PERF.13) -- the `id_blocks` table
--   below hands each writer a block of ids per table, so a sector's
--   rows can be written many at a time with their parents' ids already
--   known. `_migrate_v44_to_v45` creates it; the rows fill on first use.
--
-- v46: whole-word name search (PERF.16) -- a FULLTEXT index on `name` in
--   `sectors`, `star_systems`, `stars`, `planets` and `moons`, which
--   `queryDb`'s search and the Galaxy Map's address bar match words
--   against instead of scanning with `LIKE '%term%'`.
--   `_migrate_v45_to_v46` adds them.
--
-- v47: rogue planet classes (GEN.8) -- `rogue_planets.planet_class`, the
--   rogue's `PLANET_CLASSES` letter, drawn from the classes flagged "r"
--   (`roguePlanetData.choose_rogue_planet_class`). NULL for a brown dwarf
--   (a failed star, not a planet). `_migrate_v46_to_v47` adds the column
--   and gives every stored rogue its most probable fitting class
--   (`roguePlanetData.default_rogue_planet_class`).
--
-- v48: rogue planet surface conditions -- `rogue_planets` gains `age_gy`,
--   `internal_heat_flux_w_m2`, `effective_temperature_k`, `surface_regime`,
--   `surface_temperature_k` (a giant's at 1 bar), `surface_pressure_pa`,
--   `ice_shell_thickness_km`, `ocean_depth_km` and `has_liquid_water`
--   (`rogueSurface.rogue_surface_conditions`: with no star, a rogue's only
--   heat is its own). `_migrate_v47_to_v48` adds them and fills every
--   stored rogue, seeded by its name, and resets `has_internal_heat` to
--   the computed answer.
--
-- v49: bright-star backfill per sector block (GEN.23) -- the
--   `bright_star_blocks` table (dropped in v53 for `sector_stats`): for
--   each 3x3x3 sector block a backfill has reached, the dimmest
--   luminosity its stars have been
--   drawn down to. `_migrate_v48_to_v49` creates it empty.
--
-- v50: `system_configs` gains `comets` and `wide_binary`, so a stored
--   recipe keeps `--comets`/`--wide-binary` (both were dropped on save).
--   And a database migrated from an old version is brought
--   to exactly this file's (TEST.8 migrates every released schema and
--   compares). `_migrate_v49_to_v50` drops the placeholder DEFAULTs the
--   v9-v13 and v38 steps left on NOT NULL columns, and remakes
--   `nebulae`/`asteroid_fields`' `sector_id` foreign keys, still ON DELETE
--   CASCADE where v16/v17 made them, as ON DELETE SET NULL.
--
-- v51: `galaxy_shape.galaxy_seed`, the galaxy's 128-bit seed (GEN.39,
--   docs/design/generation/reproducible-galaxies.md): every sector, bright-star
--   layer and backfill block draws from a seed derived from it and its
--   address (`galaxySeed.unit_seed`). `BINARY(16)` rather than a
--   32-character hex string: half the bytes and no case to normalize; it
--   is shown and typed as 32 hex digits. Written once, when the galaxy is
--   first planned. `_migrate_v50_to_v51` adds it empty: a galaxy made
--   before has no seed, and none can be worked out for it.
--
-- v52: what made the galaxy and every run since (DB.6,
--   docs/design/generation/reproducible-galaxies.md): `galaxy_shape` gains the
--   22-hex-digit version key (`versionKey.version_key`), the PlanetGen,
--   Python and platform versions, written with the galaxy seed; and the
--   `generation_runs` table below holds one row per `planetgen` run that
--   changes the galaxy. A galaxy is built by a series of commands, not by
--   the seed alone, so these rows are what a rebuild replays.
--   `_migrate_v51_to_v52` adds both empty.
--
-- v53: one row of stats per sector (GEN.44, PERF.11, `sector_stats`
--   below), keyed by grid address so a sector nobody has filled can have
--   one: how deep its bright stars go (replacing v49's per-block
--   `bright_star_blocks`, which it drops), its expected density, and,
--   once filled, what it actually holds. `galaxy_shape` gains the
--   decaying average of actual against expected systems.
--   `_migrate_v52_to_v53` moves each block's level onto its unfilled
--   sectors and gives every filled sector a row at level 0.
--
-- v54: `system_configs` gains a `prevalence_<feature>` column per
--   forcing option (GEN.52): the sector or galaxy run's percentage
--   deviation from the feature's normal chance, 0 for none.
--   `_migrate_v53_to_v54` adds them at 0.
--
-- v55: `sector_stats` keeps raw statistics instead of a baked Galaxy Map
--   color (DB.14): `fill_share` and `color_r`/`color_g`/`color_b` go,
--   `mean_age_gy` and `total_luminosity_sol` come. `_migrate_v54_to_v55`
--   works the two new ones out for each filled sector from its stars.
--
-- v56: a nebula's shape is stored with it (GEN.75, `galaxy/nebula_shape.py`):
--   the single values as `nebulae.shape_*` columns (the ellipsoid's
--   stretch and turn, the warp, the isovalue, the scale that puts the
--   surface's farthest point at the radius) and its 4-8 metaballs as
--   rows of the new `nebula_shape_balls`. All NULL / absent for a nebula
--   saved before v56: `queryDb.nebula_shape` then draws the same shape
--   from the nebula's own properties, as a fresh save would have stored.
--
-- v57: `bright_stars.off_plane`, a virtual flag (1 when the star lies 250 pc
--   or more above or below the plane, `query.BRIGHT_STAR_PLANE_HALF_THICKNESS_PC`)
--   with `idx_bright_stars_off_plane (off_plane, luminosity_w)`, so
--   `queryDb.galaxy_brightest_stars` can list the most luminous stars off the
--   plane without reading every star (GEN.117). Computed from
--   `position_z_mpc`; nothing writes it.
--
-- v58: every object has a unique ID (GEN.69, GEN.68's plan, `galaxy/uid.py`):
--   a `uid` column, NULL for a row saved before v58. A sector's is its
--   designation as an integer (BIGINT UNSIGNED, unique); a star system's
--   and a phenomenon's is 96 bits (BINARY(12), unique in its table):
--   SHA-256(galaxy seed || "uid-<kind>:" parent "/" index) with the top bit
--   set, or, for an interstellar object or bright-sweep system that has one,
--   its GEN.64 position ID (76 bits, so never equal to a derived one).
--   A star's, planet's, moon's, belt's and comet's is 64 bits (BIGINT
--   UNSIGNED), unique under its system (`UNIQUE (star_system_id, uid)`). The
--   index is the object's slot in its parent as generated, so regenerating
--   in place from the same galaxy seed gives the same IDs. `store.assign_uids`
--   writes them after a sector or system is saved.
--
-- v59: `idx_bright_stars_population (population, off_plane, luminosity_w)`
--   replaces v57's `idx_bright_stars_off_plane`: the Galaxy Map picks a
--   tile's bright stars from each population in turn, so the old giants
--   of the bulge and the thick disk are listed beside the young blue
--   stars that outshine them (`queryDb.galaxy_brightest_stars`).
-- v60: velocity (GEN.121, GEN.122): `planets`, `moons` and `comets` each gain
--   `velocity_x_kms`/`_y_kms`/`_z_kms`, the body's velocity in km/s relative
--   to its orbital anchor (the same one `position_x/y/z_km` is measured
--   from), on the same axes. `advance_orbital_phases` and
--   `advance_comet_orbits` move it with the position, so it is always the
--   velocity at the stored phase or anomaly: a planet's or moon's is the
--   tangent of its circular orbit, `orbital_speed_kms` long, a comet's the
--   Kepler or Barker velocity. `_migrate_v59_to_v60` works it out for rows
--   saved before.
-- v61: `star_systems` gains `velocity_x_kms`/`_y_kms`/`_z_kms` (GEN.121), the
--   system's velocity in km/s on the galactic axes: the rotation curve's
--   tangent at its place (`galactic_orbital_speed_kms`,
--   `galaxy.system_position.galactic_velocity_ms`) plus, for a runaway or
--   hypervelocity system, `runaway_speed_kms` along a random direction.
--   `advance_galactic_positions` turns it with the position, and a loaded
--   sector puts it on the system's entry, stars and bodies. The direction
--   is not a column of its own: `_migrate_v60_to_v61` gives a system saved
--   before one from its id. Placed phenomena store none: they only follow
--   the rotation curve, so it is worked out from `galactic_orbital_speed_kms`.
-- v62: `sector_paths` and `sector_path_knots` (GEN.123): the path of each star
--   system, rogue planet and interstellar comet through its sector as
--   cubic Hermite spline knots -- see the tables' own comment.
-- v63: `generation_run_arguments` replaces the JSON `generation_runs.arguments` (DB.13).
-- v64: `facilities` gains `velocity_x_kms`/`_y_kms`/`_z_kms` (GEN.125): a
--   stand-alone facility's velocity in km/s on the galactic axes, the rotation
--   curve's tangent at its place (the same as a star system's before any
--   runaway speed). `advance_galactic_positions` turns it with the position.
--   Facilities on a body, in orbit or in a belt keep 0: their motion is the
--   `orbit_*` columns, relative to their host.
-- v65: `phenomenon_scatter` (GEN.100): the black holes, neutron stars,
--   planetary nebulae, supernova remnants, hypervelocity stars and nucleus
--   `planetgen plan` places galaxy-wide, and `galaxy_shape.phenomenon_scatter_seed`.
-- v66: every object the orbit update moves gains `epoch_unix` and an indexed
--   `next_update_due` (GEN.106), both in Unix seconds by the database
--   server's clock: `planets`, `moons`, `comets`, `star_systems` (its
--   galactic orbit; `binary_epoch_unix`/`binary_next_update_due` for a
--   pair's mutual orbit), `facilities`, and the six orbiting phenomenon
--   tables. `epoch_unix` is when the stored position holds (NULL: at
--   `orbit_simulation_state.last_updated_at`, the last update run);
--   `next_update_due` is when the object will have moved its threshold
--   (0.01 mpc on a galactic orbit, 0.01 AU in a system, 100,000 km round a
--   planet), from its speed, capped at a billion years. NULL means not yet
--   worked out: the next update fills it in. The update moves and counts
--   only rows that are due, and anything that changes an orbit sets
--   `next_update_due` back to NULL. See docs/design/orbits/orbital-updates.md
--   section 3.
-- v67: each sector records the code that generated it (DB.7): `version_key`
--   (the 22-digit key, `planetgen/galaxy/version_key.py`), `planetgen_version`,
--   `python_version` and `platform`, the same four values `galaxy_shape` keeps
--   for the galaxy's plan. NULL on a sector generated before v67. A galaxy
--   extended by a different release has sectors with different keys: it
--   reproduces only sector by sector, and the CLI and the Generate page warn
--   before extending it (`galaxy/version_check.py`).
-- v68: every rotating object stores its spin (GEN.104): `spin_axis_x/y/z`, a
--   unit vector, and `axial_tilt_deg`, its angle from the orbit normal (from
--   the galactic pole for stars, remnants, rogue planets and interstellar
--   comets), on `stars`, `planets`, `moons`, `comets`, `rogue_planets`,
--   `interstellar_comets` and standalone `black_holes`/`neutron_stars` (an
--   anchored remnant's is on its `stars` row). `stars`, `comets`,
--   `rogue_planets` and `interstellar_comets` gain `rotation_period_hours`;
--   the spin vector is the axis times 2 pi over that period. NULL on a row
--   generated before v68. Draws: `planetgen/physics/spin.py`, and
--   docs/design/orbits/orbital-updates.md section 6.
-- v69: `phenomenon_scatter.epoch_unix` (GEN.137): the orbit epoch the scatter
--   was drawn at, Unix seconds. A hypervelocity star travels in a straight
--   line, `p0 + v (t - epoch_unix)` (`galaxy/straight_line.py`); NULL on a row
--   drawn before v69, which holds at the database's orbit epoch.
-- v70: planets and moons store their mantle redox and their air (GEN.85):
--   `mantle_redox` ('reduced' at or below the iron-wustite buffer, 'oxidized'
--   from IW+2.5, else 'intermediate'; NULL for a gas giant), `mantle_delta_iw`
--   (log units above IW) and `p_<gas>_kpa`, the partial pressure in kPa of
--   O2, CO2, CO, N2, Ar, H2, H2O, CH4, H2S and SO2 (0 for an airless body).
--   The `atmosphere` text is written from them; a gas outside the ten (helium,
--   ammonia, sodium) is named there but not stored. NULL on a row generated
--   before v70. Draws: `planetgen/physics/atmosphere.py`.
-- v71: `galaxy_shape.phenomenon_min_mass_solar` (GEN.167): the lowest mass
--   of the neutron stars and black holes the phenomenon scatter placed
--   (`planetgen plan --phenomenon-min-mass`, default 8 solar masses). A
--   sector fill draws those below it itself, from its own stream (GEN.168).
--   NULL on a scatter drawn before v71, which placed every mass.
-- v72: stellar activity and planetary magnetic fields (GEN.86). `stars`:
--   `log_lx_lbol` (log10 of the coronal X-ray share of the luminosity; NULL
--   for a white dwarf or compact remnant), `l_xuv_w` (X-ray plus EUV, W,
--   including a hot photosphere's own output above 13.6 eV),
--   `xuv_saturated` (the corona still saturated), `flare_n33_per_yr`
--   (flares above 1e33 erg a year), `flare_alpha` (the flare-energy slope)
--   and `xuv_fluence_j` (the XUV given off over the star's age, J).
--   `planets` and `moons`: `magnetic_moment_a_m2`, `dipole_class` ('none',
--   'weak', 'earth-like', 'strong', 'multipolar'), `magnetopause_rp`
--   (standoff in body radii; NULL without a field), `xuv_flux_earth` (XUV
--   flux in units of Earth's today), `xuv_exposure_index` (lifetime XUV
--   against a 1 solar mass star's at 1 AU over 5 Gyr) and
--   `flare_irradiation_index` (flares above 1e33 erg a year over distance
--   squared, AU). NULL on a row generated before v72. Draws:
--   `planetgen/physics/activity.py` and `magnetism.py`; design
--   docs/design/habitability/activity-magnetism-radiation-hydrosphere.md sections 2 and 3.
-- v73: `galaxy_shape.bright_star_mass_limit_sol` (GEN.185): the star scatter
--   places in advance every star born with at least this mass (its mass
--   pass, the same mass limit as the phenomenon scatter's), then the
--   lighter stars at least `bright_star_min_luminosity_sol` bright (its
--   luminosity pass), skipping any sector the mass pass gave a star that
--   bright. A sector's own draw is lighter than the limit and dimmer than
--   the floor. NULL on a scatter drawn before v73, which had no mass pass.
--
-- v74: hydrosphere and ocean chemistry (GEN.88). `planets` and `moons`:
--   `water_mass_fraction`, `hydrosphere` ('dry', 'vapour', 'ice',
--   'ice-covered ocean', 'surface ocean', 'hycean'), `ocean_fraction` and
--   `land_fraction` (shares of the surface; water counts liquid or
--   frozen), `ocean_depth_km` (mean liquid depth where there is water),
--   `ice_shell_km` (the ice lid), `hp_ice_km` (high-pressure ice under the
--   deepest possible liquid), `ocean_class` ('ice-sealed', 'chloride
--   brine', 'acid sulfate', 'soda', 'neutral'), `ocean_ph`,
--   `water_activity` and `phosphorus` ('high', 'limited', 'starved'). All
--   NULL for a gas giant; the ocean fields NULL without liquid water.
--   `rogue_planets` gains `hp_ice_km`, and its `surface_regime` may be
--   'hycean'. NULL on a row generated before v74. Draws:
--   `planetgen/physics/hydrosphere.py`; design
--   docs/design/habitability/activity-magnetism-radiation-hydrosphere.md section 5.
--
-- v75: surface radiation dose, UV and the galactic hazard (GEN.87).
--   `planets` and `moons` (NULL for a gas giant): `surface_dose_msv_yr`
--   (the yearly dose at the surface, mSv, = `dose_gcr_msv_yr` x
--   `dose_helio_mult` + `dose_sep_msv_yr` + `dose_ground_msv_yr`: cosmic
--   rays after the air and any dipole, the compressed-heliosphere
--   multiplier 1 to 2.5 that `refresh_containment` updates when a nebula or
--   remnant presses on the star, stellar particles, and crust plus radon),
--   `uv_surface_index` (DNA-weighted UV against Earth's 1, after the ozone
--   layer) and `ozone_loss_flag` (ozone present but a galactic event rate
--   over 1 per 100 Myr, or a stellar-particle dose of 100 mSv/yr at the
--   ozone layer). `stars` gains `lethal_event_rate_per_gyr` (supernovae
--   within 8 pc per Gyr, from the distance to the galactic centre). NULL
--   on a row generated before v75. Draws: the crust abundance
--   (`planetgen/physics/radiation.py`); design
--   docs/design/habitability/activity-magnetism-radiation-hydrosphere.md section 4.
--
-- v76: the habitability score of every planet and moon (GEN.89).
--   `planets` and `moons` (NULL for a gas giant): `phi4` (the PHI-4 display
--   score, the geometric mean of the four domains), `phi4_pressure`,
--   `phi4_temperature`, `phi4_chemistry` and `phi4_radiation` (each 0 to 1) with
--   `tier_<domain>` (0 Blue, 1 Green, 2 Yellow, 3 Red) and `equipment_tier`
--   (0 ideal, 1 mask, 2 mask and scrubber, 3 environment suit, 4 full life
--   support), then `phi_bio` (microbial life), `phi_cpx` (complex life) and
--   `phi_tech` (human operability with equipment) and the four likelihoods
--   behind them (`l_solv`, `l_chem`, `l_ener`, `l_rad`). `hab_note` names the
--   domains short of Blue, or a compact host (pulsar, neutron star, black
--   hole) whose planet is rated at a lethal dose; `energy_flux_w_m2` is the
--   light or chemical power life can use. Re-scored by `refresh_containment`
--   when the dose changes. NULL on a row generated before v76.
--   `planetgen/physics/habitability_world.py`; design
--   docs/design/habitability/habitability-index.md section 8.
--
-- v78: object IDs (DB.20, GEN.170, `galaxy/object_uid.py`) replace v58's
--   hashed `uid`s. Every object table (`star_systems`, `stars`, `planets`,
--   `moons`, `asteroid_belts`, `comets`, the eight phenomenon tables and
--   `facilities`) has `uid BINARY(10)` (BINARY(12) since v93), UNIQUE on its own: 80
--   bits, the birth sector's address (40), a serial in that sector (28; the top two
--   bits say generated, run-time or field-drawn) and a body number (12; 0
--   for the top-level object, 1 and up for the stars, planets, moons,
--   belts and comets born in a system). Given once at birth and never
--   recomputed. `id_counters` holds the counters of run-time births: per
--   sector for serials, per system for body numbers. The row `id` stays the
--   foreign key. Existing rows are numbered by row order.
--
-- v79: the mass backfill (GEN.187). `sector_stats.bright_mass_sol`: the
--   initial mass, in solar masses, down to which a backfill placed every
--   living star of an unfilled sector (the sectors a face away from the
--   generated ones: 1, the next ring 2, then 5, then 8); NULL where none
--   did, so the sector holds the galaxy scatter's mass pass only (stars
--   born at or above `galaxy_shape.bright_star_mass_limit_sol`). Alongside
--   it `bright_level_sol` is the luminosity pass's floor, as before. A new
--   plan scatter clears it. `planetgen/generation/run_galaxy.py`
--   (`backfill_bright_stars_around`); design
--   docs/design/generation/mass-backfill.md.
--
-- v80: `phenomenon_scatter.mass_solar` (MAP.165): the solar mass each scattered
--   black hole (the nucleus too) and neutron star was drawn with, so the Galaxy
--   Map sizes it by mass and the sector build gives the object that same mass;
--   NULL for a quasar, hypervelocity star, nebula or remnant, and on a row
--   scattered before v80 (the map then sizes it by class, the build draws).
--
-- v81: galaxy shape density settings (ADM.49). `galaxy_shape.arm_level` (the
--   thin disk's level, the middle of the arm crest and the inter-arm trough:
--   the arm factor is `arm_level * (1 + arm_amplitude * cos)`, so the crest is
--   `arm_level * (1 + arm_amplitude)` and the trough `arm_level * (1 - arm_amplitude)`)
--   and `core_amplitude` (the galactic core's centre density, a Gaussian on top
--   of the bulge; 0 is none; replaced by `nuclear_amplitude` in v99). Defaults 1 and 0 are the model before v81.
--
-- v82: indexes and a summary table for the pages that timed out on a big galaxy
--   (PERF.68, PERF.70, PERF.74). `phenomenon_scatter` gains
--   `idx_phenomenon_scatter_class (kind, subtype, mass_solar)`, which serves the
--   Galaxy Map's coarse tiles (the few biggest classes) without reading the
--   address range; `star_systems.idx_star_systems_quadrant` becomes
--   `(quadrant, name, id)` and `idx_star_systems_binary_name (is_binary, name,
--   id)` is new, so the Systems list sorted by Octant or Binary reads one page
--   of rows per group; and `sector_system_counts` stores each sector's system
--   count, refreshed in the background by `query.list_sectors` (never more often
--   than `countcache` allows), so the Sectors list does not count every system
--   on each request.
--
-- v83: `nebula_field` (DB.24): one row per molecular cloud of classes M and N
--   the galaxy's seed places (`galaxy/nebula_field.py`), written at plan time,
--   so the Galaxy Map can find a nebula no generated sector holds yet. Seed
--   fixed, so nothing in a fill changes it except `is_stored` and `uid`, set when a
--   sector's save stores the cloud as a `nebulae` row.
--
-- v84: `phenomenon_scatter.map_luminosity_sol` (DB.23, GEN.200): the luminosity, in
--   solar luminosities, of a main-sequence star of the row's `mass_solar` by the
--   map's own mass-luminosity relation (`stellar_evolution.map_luminosity_sol`), so
--   the Galaxy Map can show a faint but massive neutron star or black hole where a
--   star of its mass would be visible. FLOAT, NULL for a row without a mass
--   (planetary nebulae, supernova remnants, hypervelocity stars, the quasar).
--   Written when the scatter inserts a row; `planetgen fill-map-luminosity` fills it
--   in a galaxy scattered before v84.
--
-- v85: `object_timing_stats` and `object_timing_samples` (PERF.92): how long each
--   kind of generated object takes, per run and per sector (count, total, minimum,
--   maximum and a 24-bucket histogram), plus a repeatable 1-in-1,000 sample of
--   the objects with their full timing; a row for every object only in a
--   benchmark run (`PLANETGEN_OBJECT_TIMING_FULL=1`). Written by
--   `generation/object_timing.py` with each sector a run saves.
--
-- v86: `galaxy_block_stats` (DB.25): the sums a drill-down stage shows (sectors,
--   systems, expected systems, stars, star age sum, luminosity) stored per block of
--   each ladder level above the sector, so a Galaxy Map window reads its cells in one
--   statement (`query.galaxy_block_stats`). Kept by `db/block_stats.py`.
--
-- v87: the Generation lane's columns on `planets` and `moons` (DB.26). GEN.178:
--   `induced_moment_a_m2` (the induced magnetic moment, after `dipole_class`; NULL
--   when not computed or none). GEN.148: `flux_variation` (the relative swing of
--   the stellar flux over the orbit or pair cycle, after `energy_flux_w_m2`; 0 is
--   steady). GEN.207 (PHI-6): `phi4_chemistry` and `tier_chemistry` are renamed
--   `phi4_air` and `tier_air`, and `phi4_water`, `tier_water`, `phi4_suitability`
--   and `tier_suitability` follow them. GEN.150: `nebulae` gains `age_myr` (Myr),
--   `q_total` (total ionizing photons per second) and `hosts` (host stars as counts
--   by spectral type, "O5V x1, O9V x4"), after `extinction_av`. All NULL on a row
--   generated before v87.
--
-- v88: the tech-level columns of `species` (POP.9). `energy_index`, `materials_index`,
--   `information_index`, `medical_index`, `propulsion_index` and `defense_index` (each
--   0 to 7, NULL until drawn), `tech_level` (DECIMAL(3,2), the weighted sum, indexed
--   for range search) and `tech_source` ('generated' or 'admin': a pass never
--   redraws a hand-edited species). NULL on a species saved before v88.
--
-- v89: facility types and the facility columns that go with them (DB.28, POP.10). A
--   *type refines a kind*: `facility_types` (name, kind, icon, description, order,
--   active) holds the admin-defined types; `facility_type_fields` the typed fields
--   a type carries (key, label, datatype int/decimal/text/bool/choice, unit,
--   required, bounds, choices, default; at most 24 per type and the key's pattern
--   `^[a-z][a-z0-9_]{0,31}$` are checked by the application); `facility_field_values`
--   one typed value per facility and field (retire a field, do not delete it: the
--   field FK is RESTRICT). `facilities` gains `type_id` (RESTRICT), the affiliation
--   (`affiliation_polity_id`, SET NULL, or free text `affiliation_name`), the five
--   Green/Yellow/Red ratings `rating_crime`, `rating_housing`, `rating_resources`,
--   `rating_maintenance` and `rating_health` (1 Green to 3 Red, NULL until rated),
--   and how they were set (`ratings_source`, `ratings_note`, `ratings_set_at`).
--   All NULL on a facility saved before v89. Design:
--   docs/design/habitability/population-and-politics.md, "Data model" of the facility section.
--
-- v90: `galaxy_block_stats.ladder_step` (MAP.151b): the step between the block ladder's
--   levels (`galaxy.drill.DRILL_STEP`, 3) the row's totals were computed for. The ladder
--   went from 243/27/3 to 243/81/27/9/3 and the same `m` means another block on each, so
--   the map reads only rows of its own step. The migration empties the table and fills it
--   again for the 3-ary ladder.
--
-- v91 (Spherical Sector Coordinates, DB.32): the sector grid is spherical shells
--   (`galaxy/spherical.py`, docs/design/galaxy/spherical-sector-coordinates.md). A sector's
--   address `(ring_index, layer_index, ring_slot_index)` becomes `(shell_index,
--   band_index, wedge_index)` in `sectors`, `bright_stars`, `sector_stats` and
--   `phenomenon_scatter` (same types); `galaxy_shape.outer_ring_index` becomes
--   `outer_shell_index`; `galaxy_layer` and `galaxy_column` are replaced by
--   `galaxy_shell` (one row per shell: the highest band its outline reaches and the
--   cell count of that run); `galaxy_block_stats` is dropped (the Galaxy Map's block
--   ladder is rebuilt on shells in v92). Every address means something else, so the migration CONVERTS:
--   every system, phenomenon and facility keeps its place, velocity and epoch and is
--   linked to the new sector that holds it, the sector rows are rebuilt on shells
--   (a cell that holds an object, or whose centre lay in a generated sector), object
--   IDs are renumbered for the new birth sector (body numbers kept), the scatter and
--   bright stars take the address of their cell, `sector_stats` is rebuilt for the
--   generated sectors, and the outline is built again from the stored shape.
--
-- v92: `galaxy_block_stats` again, for the block ladder on shells (MAP.205): keyed
--   `(m, shell_index, band_index, wedge_index)`, no `ladder_step` (the ladder is the
--   code's `galaxy/drill.py`; a ladder change is a migration that refills the table).
--
-- v93 (GEN.230, GEN.233, GEN.223): every `uid` is `BINARY(12)`, the 96-bit object ID
--   (19, 18 and 19 bits of shell, biased band and wedge, a 28-bit serial and a 12-bit
--   body number): the 80-bit ID's 12-bit band field could not name a polar cell past
--   shell 1,303. The migration renumbers by moving each field to its new place (the
--   serial and body are untouched), in the object tables, `nebula_field`,
--   `id_counters` (the scope of a sector's counter is 7 bytes now) and `sectors.uid`
--   (the designation; a sector with no grid address counts up from `1 << 56`).
--   `galaxy_shape.reach_shell_index`: the outermost shell the grid has an address for,
--   by default 65,535 (262 kpc); the outline (`outer_shell_index`) still says where the
--   disk's content ends. `bright_stars` and `phenomenon_scatter` store their point in
--   space as every other object does: position in milliparsecs (now DOUBLE, so a slow
--   drift is not rounded away), velocity vector in km/s, `epoch_unix` and an indexed
--   `next_update_due`; the orbit update moves them while they are unbuilt. Both also gain
--   `member_count` and `system_mass_solar` (NULL for a single object) and the new tables
--   `scatter_members` and `scatter_nodes` hold the companions of a scattered system (multi-star
--   note, sections 5, 13, 14); nothing writes them yet. Every table that holds a point in
--   the galaxy frame also has `cell_key`, a generated Morton key of its 0.25 pc cube, and the
--   covering index (cell_key, x, y, z) (region queries; `galaxy/cellkey.py`); `star_systems`
--   gets `galaxy_x/y/z_pc` for it, since its own position is relative to its sector.
--
-- v94 (GEN.234): `bright_stars.population` also takes 'halo' and 'cluster' (the stellar halo
--   replaces the flat density floor; the globular clusters get their own population).
--
-- v95 (GEN.158): `stars.metallicity_fe_h` and `bright_stars.metallicity_fe_h` (DOUBLE, nullable), a star's [Fe/H]
--   in dex against the Sun (`planetgen/physics/metallicity.py`); NULL on rows made before v95.
--
-- v96 (GEN.153): `neutron_stars.pulsar_type` also takes 'magnetar' (its CHECK, `neutron_stars_chk_1`, is
--   replaced; the quadrant CHECK stays `neutron_stars_chk_2`).
-- v97 (GEN.159): `galaxy_components`, the discrete stellar components of a galaxy (globular clusters now, satellites
--   and streams later) as density profiles with a position, velocity, mass, age and [Fe/H], and a `cell_key` like the
--   other tables with a point in the galaxy frame.
-- v98 (GEN.281): `galaxy_shape.cell_star_ceiling` (default 50,000 systems) and `galaxy_shape.fill_luminosity_cap_sol`
--   (default 100 solar luminosities): a cell that expects more systems than the ceiling is filled only with the stars
--   at or above the cap; the rest is computed from the address and never stored.
-- v99 (GEN.280): `galaxy_shape.core_amplitude` (ADM.49's Gaussian core) becomes `nuclear_amplitude`, how much of the
--   observed nuclear star cluster and nuclear stellar disc the centre holds (1 the observed pair scaled to the bulge,
--   0 none); `bright_stars.population` (now VARCHAR(16)) also takes `nuclear_cluster` and
--   `nuclear_disc`.
-- v101 (DB.54): `galaxy_shape.core_min_luminosity_sol` and `core_mass_limit_sol` (NULL until set: the core follows the
--   baseline scatter limits times 1.1) and `halo_min_luminosity_sol` (default 500, at least 50), the scatter limits by
--   place that the tiers of GEN.333 read.
-- v100 (GEN.281): `galaxy_shape.cell_star_ceiling` defaults to 5,000 systems (a cell of 50,000 took 14 minutes) and has a
--   CHECK, `chk_galaxy_shape_cell_star_ceiling`, of at least 1.
--
-- v102 (DB.55): the point of every kind of object in one additive revision (`planetgen/db/point_columns.py`, GEN.294).
--   A movement vector (`velocity_x/y/z_kms`, galactic axes, backfilled from the old scalar speed as the tangent of the
--   galaxy's rotation) on black holes, neutron stars, rogue planets, interstellar comets, asteroid fields, nebulae,
--   remnants, belts and quasars; `mass_kg` on interstellar comets, asteroid fields, nebulae, remnants, belts, comets and
--   facilities (a comet nucleus at 67P's density, a cloud's bounding sphere at its density, else a stated zero);
--   `epoch_unix` and `next_update_due` on belts and quasars, `epoch_unix` on `scatter_members`; a system-frame
--   `position_x/y/z_km` on belts; `sectors.max_mass_kg` and `max_mass_id`; `asteroid_fields.formed_at` and
--   `origin_event_id`; and `parent_type`, `parent_id` and `dynamical_state` on planets and moons (GEN.289). All
--   nullable except the planet and moon state, which default to a bound body of its star (planet) or planet (moon).
--   (A star's own columns came with v104.)
--
-- v103 (GEN.336): the scalar `galactic_orbital_speed_kms` is gone from `stars`, `black_holes`, `neutron_stars`, `nebulae`,
--   `supernova_remnants`, `rogue_planets`, `interstellar_comets` and `asteroid_fields`, and `binary_galactic_orbital_speed_kms`
--   from `star_systems`. The movement vector of v101 (`velocity_x/y/z_kms`) is the one record; a speed is its length. A star
--   takes its system's vector.
--
-- v104 (DB.53): the built system's orbit tree. `orbit_nodes` (key `star_system_id`, `node_no`; the columns of
--   `scatter_nodes` and a disk descriptor); on `stars` `parent_node_id`, the member columns of `scatter_members`
--   (`state`, `flags`, `offset_x/y/z_au`, `velocity_x/y/z_kms`, `decay_age_yr`, `epoch_unix`) and the disk descriptor
--   (`disk_kind`, `disk_mass_solar`, `disk_inner_au`, `disk_outer_au`); on `planets` `host_node_id`. A close or wide pair
--   becomes one node from its stored mutual orbit, and each star's offset is its stored binary position.
--
-- v105 (GEN.324): `bright_stars.population` also takes `satellite`, the stars of a satellite galaxy placed by the component
--   pass (a cluster's are `cluster`, v94).
--
-- v106 (DB.65): every generated `cell_key` holds a cube number of at most 2^21 - 1 (`LEAST(..., 2097151)` on each axis), so the
--   last double below the 262,144 pc edge, which rounds up to 2^21, no longer wraps to the opposite edge's key.
--
-- v107 (GEN.526): `bright_stars.population` (now VARCHAR(24)) also takes `satellite_young` and `satellite_intermediate`, the stars a
--   satellite galaxy's history puts in the last 0.1 Gy and the 3 Gy before that (`satellite` is the rest).
--
-- v108 (DB.72): `band_index` is INT in `bright_stars`, `sector_stats` and `phenomenon_scatter` (it was SMALLINT, which stops at
--   32,767; a shell's bands reach +-32,768 from shell 20,861 and +-102,942 at shell 65,535).
--
-- v109 (GEN.289): `orbital_phase_deg` is gone from `planets` and `moons`. A planet or moon is its state vector (position,
--   velocity, mass, epoch, parent and dynamical state, the columns of v102) and moves by `kepler.universal_step`; the
--   phase was a second record of the same place on a circle. A circumbinary planet's parent is its system.
--
-- v110 (DB.50): `planets` and `moons` gain `life_stage`, `evolutionary_pace`, `life_milestone_age_gyr`, `life_system_age_gyr` and
--   `flavor_key`, the structure that replaces the evolutionary paragraphs and `flavor_text`. NULL until generation writes them
--   (DB.51); nothing is dropped.
--
-- v111 (DB.52): the prose is gone. `planet_evolutionary_paragraphs`, `moon_evolutionary_paragraphs`, `flavor_text` and `hab_note` (on
--   `planets` and `moons`) and `system_configs.markdown` are dropped: the life structure of v110 replaced the first three, the
--   tier columns the note, and nothing read the flag. `flavor_text_count` stays (the generator still counts the line).
--
-- v112 (GEN.244): `bright_stars.population` also takes `open_cluster_00` to `open_cluster_11`, the age bins (a quarter dex wide,
--   1 Myr to 1 Gyr) of a generated galaxy's young open clusters.
--
-- v113 (DB.79): `object_index`, one row per object of every kind the site lists (kind and 96-bit object ID, name, sector, system,
--   parent, class, mass, radius, placed), with keyset keys per kind and sort and a FULLTEXT index on the name. Empty until the
--   index builder (PERF.140) fills it; no other table changes.
--
-- v115 (DB.44): `sectors.storage_checksum` (BINARY(32), NULL) and `sectors.parity_dirty` (TINYINT NOT NULL DEFAULT 1, so a sector that
--   existed before counts as not yet protected), the two columns DB.9's storage protection reads. Nothing reads them yet.
--
-- v116 (MAP.172): `nebula_block_stats`, `nebula_block_classes` and `nebula_block_state`, what the nebulae add up to in each block of the Galaxy Map's block
--   ladder (`galaxy/drill.py`, levels 243 to 3): count, volume, centroid, bounding radius, the volume of each family, the dominant
--   class and the cover, kept by `db/nebula_blocks.py`. Made from `nebula_field` and the stored nebulae; no sector is needed.
--
-- v117 (DB.66): `galaxy_shape.family` (`spiral`, `lenticular`, `elliptical`, `irregular`, `dwarf_spheroid`) and the columns of
--   the components a family may have that a spiral does not draw today: the bar's axis ratios and angle, the thick disk's
--   amplitude, the Sersic spheroid (amplitude, index, effective radius, axis ratios, three orientation angles) and the
--   irregular family's clump share. Existing galaxies become spirals with no spheroid; nothing reads the columns yet.
--
-- v118 (DB.56): `orbital_events` (one row per orbital event: a collision outcome, merger, split, ejection, capture, unbind or warning, with the
--   numbers it was decided from, the conservation check, the products and the flags of docs/design/orbits/collisions-and-mergers.md section 7) and
--   `orbital_event_bodies` (its participants and products, one row each). Nothing writes them yet (GEN.197's `record_event` does).
--
-- v119 (DB.30): `species_limit_values`, the environment limits a species does not share with the human visitor: one row per limit it
--   replaces (species, the limit's name in `protection_limits.LIMIT_ROWS`, the value in its unit). A limit with no row is the human
--   value. Empty until the species editor (POP.13) writes it; no other table changes.
--
-- v120 (DB.49): `sectors.leaf_digest` (BINARY(32), NULL), the sector's fingerprint stored at save time (NULL for a sector saved
--   before it), and `sector_tree_nodes`, the digests above the leaves: a band node (level 1: shell and band, over the wedges'
--   leaves) and a shell node (level 2: shell, band_index 0, over its band nodes), each with the count of what is under it and a
--   dirty flag; the galaxy root is worked out on demand. Nothing reads either until GEN.305.
--
-- v121 (DB.16): a sector records `generator_epoch` (SMALLINT UNSIGNED, NULL: unknown before this) and `run_id` (the `generation_runs` row of the run that
--   saved it) instead of the four version text columns of v67, which are gone: the version key, release, Python and platform are looked up from the run.
--   A sector with no recorded version (made before v67) has a NULL `run_id`.
--
-- v122 (DB.36): `change_plans` (who asked and in which role, the request, the state, the role that must approve, who decided and
--   when, the job that applied it), `change_plan_ops` (a plan's operations in order: kind, target object ID, before and after as
--   JSON, status) and `edit_journal` (the append-only list of every applied operation: plan, object ID, sector, aspect, before,
--   after, generator epoch), for the 8.6 edit-everything plans (docs/design/interface/edit-everything.md). Nothing reads them yet.
--
-- v123 (DB.37): `object_locks`, a standing keep mark on an object (`*`) or one aspect of it, one row per object and aspect, with who locked it,
--   when and why, for the 8.6 edit-everything plans (docs/design/interface/edit-everything.md). Nothing reads it yet.
--
-- v124 (DB.38): `region_density_overrides`, the append-only list of an admin's density multipliers for a region (its shape and parameters as JSON, the
--   multiplier, who set it and when, and when it was undone), applied when sectors there are filled (MAP.119, ADM.57). Nothing reads it yet.
--
-- v125 (DB.60): `sectors.density_multiplier_used` (DECIMAL(10,4), NULL) and `sectors.override_id` (BIGINT UNSIGNED, NULL, the `id` of the
--   `region_density_overrides` row it came from), so a rebuilt galaxy can tell which density override a sector was filled under; both NULL for a
--   sector filled with none. Nothing writes or reads them yet.
--
-- v126 (DB.39): `request_comments` (the conversation on an edit request: plan, author and role, body) and `request_notifications` (a notice owed to an
--   account about a request, and when it was sent), for the 10.0 roles; a request is a `change_plans` row (DB.36). Nothing reads them yet.
--
-- v127 (DB.61): `galaxy_components.name` is nullable (a cluster nobody has named, the unique key stays over the named rows), and `stars` and
--   `bright_stars` gain `cluster_id`, a nullable indexed reference to `galaxy_components.id` (cleared if the cluster is removed), so a star can
--   point to the cluster that is its parent. Every existing cluster keeps its name; nothing writes or reads the link yet.
--
-- v128 (DB.64): `constellations`, the stored constellations of a system's sky (VIEW.3): the system's uid and the constellation's number in it, the object
--   ID of its name, the member stars and the edges between them (JSON), and the seed and version key it was built under. Nothing writes or reads it yet.
--
-- v129 (DB.63): `vessels`, the ships and stations (9.3): object ID, kind, name, whether generated or user-made and its owner, the position, point mass and
--   movement vector of GEN.294 at an epoch, the course (JSON) and the seed and version key it was made under. Nothing writes or reads it yet.
--
-- v130 (DB.85): `galaxy_containers`, each globular cluster and each superbubble stored as a container: its object ID (`uid`), kind, name, the cluster component or the bubble's field address that
--   made it, its centre and extent in the galaxy frame, and the sector its centre is born in. The plan writes them from the galaxy seed; `db/containers.py` asks which hold a point and what one holds.
--
-- v131 (OPS.58): `galaxy_shape` also keeps the generator epoch, the fingerprint format and the generator battery digest (with its cases) of the code that planned the galaxy, and `generation_runs`
--   the epoch each run generated under, so a galaxy says which epoch each pass of it was made at. `planetgen.galaxy.version_check` compares them with the running code.
-- v132 (GEN.458): `stars.snow_line_km`, the star's own snow line (drawn once, 0.2 dex around the 170 K blackbody position, then stored), and `planets.born_beyond_snow_line`
--   and `moons.born_beyond_snow_line`, whether the body formed beyond its star's line (a moon takes its planet's). The fixed 2.7 AU constant is gone from the generator.
-- v144 (DB.69, GEN.428): `planets` and `moons` gain the bulk composition: `metal_mass_fraction`, `rock_mass_fraction`, `ice_mass_fraction` and `envelope_mass_fraction`
--   (with the existing `water_mass_fraction`, the five shares of the whole body's mass add up to 1), the rock's `rock_feo_fraction`, `rock_al2o3_fraction` and
--   `rock_mg_si_ratio`, the envelope's `envelope_helium_fraction` and `envelope_z_fraction`, a small body's `porosity`, `carbon_rock` and `composition_origin`
--   ('drawn' or 'inverse'). The mass, radius and density of a rocky body follow from them (`physics/bulk_structure.py`). NULL on a row saved before v145.
-- v133 (GEN.630): `stars.role` also takes 'member' (a star of a system of more than two), and `stars.slot` says which child (0 or 1) of its parent node a star is.
--   A pair already stored gets its slots from its roles. Nothing else changes.
-- v134 (GEN.267): `system_configs` also keeps the star count a system was forced to (`stars`), the cap on it (`max_stars`, 12 unless raised) and the prevalence of five or more stars
--   (`prevalence_multi_star`), so a stored config regenerates with the same count.
-- v135 (GEN.413): `planets.ozone_column_du` and `moons.ozone_column_du`, the ozone column in Dobson units (NULL for a gas giant and for a row made before v135; 0 with no oxygen),
--   drawn from the oxygen and the star's temperature (`physics/ozone.py`); the UV index reads it.
-- v137 (GEN.412): `planets` and `moons` keep the Bond albedo (`bond_albedo`) and the greenhouse quantity (`greenhouse_factor`, G = (T_s / T_e)^4 - 1) their surface temperature is made of,
--   and `albedo_basis` (0 recorded at generation, 1 the exact inverse of a row saved without them, 2 a seeded class draw where the inverse had no solution).
-- v150 (GEN.649): `planets` and `moons` keep four life scores, `life_microbial`, `life_flora`, `life_fauna` and `life_sentient` (0 none to 9 abundant, NULL when the body has no score),
--   in place of `phi_bio`, `phi_cpx` and `phi_tech`. A row saved before gets them from its two likelihoods and its evolutionary path.
-- v152 (GEN.700): `planets` and `moons` keep the tide heat flow (NULL for a planet, which has none) (`tidal_flux_w_m2`) of an icy moon under a frozen surface, which the chemical energy of its ocean's life counts. A row saved before has NULL.
-- v153 (GEN.699): `planets` and `moons` also keep the two draws of the mantle's rheology (`mantle_eta0_pa_s`, `mantle_activation_kj_mol`) and the heat that comes out, `internal_heat_flux_w_m2`
--   (radiogenic plus primordial, `physics/interior_heat.py`). `urey_ratio` is now an output of that thermal history, not a draw. NULL for a gas giant and a row made before v153.
-- v156 (DB.100): `planets` and `moons` keep the world's whole milestone schedule, the age its pace reaches each of the five milestones (`life_age_abiogenesis_gyr` to `life_age_technological_civilization_gyr`, Gy; NULL with no life chemistry),
--   so the life paragraph and the timeline paragraph are written from the row alone. A row saved before gets the ages of its stored pace.
-- v157 (GEN.654, DB.92): `planets` and `moons` keep the aquifer a world with no liquid at its surface holds underground, `aquifer_depth_km` (to the top of the water-bearing layer) and `aquifer_thickness_km` (NULL for a world with none, a gas giant and a row made before v157).
-- v160 (GEN.732): `planets` and `moons` keep the ring system of a body (`ring_kind` dusty, narrow or broad; `ring_origin`; `ring_inner_km` and `ring_outer_km`; `ring_mass_kg`;
--   `ring_optical_depth`, a mean; `ring_ice_fraction`). NULL for a body with no ring, a moon, and a row made before v160, which gets its rings when the galaxy is regenerated.
-- v162 (DB.96): `planets` and `moons` keep the colour of the body's air as one integer (`atmosphere_colour`, sRGB as red*65536 + green*256 + blue; NULL for a body with no air, no star
--   temperature or no mixture, and a row made before v162, which gets it when the galaxy is regenerated). The cloud cover it was worked out with is `cloud_cover` (v138).
-- v171 (GEN.836): `rogue_planets` keeps the same seven ring columns as `planets` (v160); a moon may hold a ring system too.
-- v168 (DB.108): the mass of each generated sector and of the blocks and regions above it: `sector_stats.total_mass_kg`, `galaxy_block_stats.mass_sum_kg` and `region_stats.mass_sum_kg`, kilograms. See `db/sector_mass.py`.
-- v169 (GEN.721): `planets` and `moons` also keep `lid_regime`, 'mobile' (plate tectonics) or 'stagnant' (one rigid plate), drawn once after the
--   hydrosphere (`physics/interior_heat.py`); the heat flux law, the dynamo lifetime and the ground dose's felsic crust read it. NULL for a gas giant
--   and a row made before v169, whose heat flux was a stagnant lid's.
-- v192 (GEN.1057): `planets`, `moons` and `rogue_planets` also keep `ring_age_yr`, the years since the event that made the ring (`physics/ring_origin.py`); a ring is now an outcome of the body's history
--   (a moon that fell inside the Roche limit, a comet that passed inside it), so a body with neither has none. NULL for a body with no ring, a resonance ring and a row made before v191.
-- v191 (GEN.1034): `planets` and `moons` also keep the open-ocean equilibrium tide of a world with a sea (`physics/ocean_tides.py`): `tide_range_spring_m` and `tide_range_neap_m` (DOUBLE, metres,
--   0 for an ice-covered ocean, NULL for a body without a free sea) and `ocean_liquid` (VARCHAR(40), the liquid of that sea, water until GEN.1035). NULL for a row made before v190.
-- v170 (GEN.767): `planets` and `moons` also keep the day and night shares and the hours of daylight (`physics/daylight.py`): `day_share_equator`, `day_share_mean`,
--   `parent_shadow_share` (a moon in its planet's shadow), `solar_day_hours`, `daylight_hours` and `polar_circle_deg`. NULL for a row made before v170.
--
-- v136 (GEN.655, DB.93): the survival equipment line replaces `planets.equipment_tier` and `moons.equipment_tier`: `eq_body`, `eq_body_kind`, `eq_body_rec`, `eq_body_rec_kind`, `eq_rad_suit`,
--   `eq_rad_suit_rec`, `eq_mobility`, `eq_dust`, `eq_gas_filter`, `eq_co2_scrub`, `eq_supply`, `eq_supply_rec`, `eq_shelter` and `eq_shelter_rec` (TINYINT, the levels of
--   `physics/equipment.py`, required and recommended) and `shelter_column_g_cm2` and `outside_hours` (DOUBLE); all NULL for a gas giant and an unscored body.
-- v137 (GEN.656): `planets` and `moons` also keep where the mean surface temperature is higher or lower (`physics/regional_temperature.py`, kelvin, NULL for a gas giant and an unscored body):
--   a world that turns on its own has the annual means of its latitude bands `temp_polar_k`, `temp_north_k`, `temp_equator_k` and `temp_south_k` and the day and night swing
--   `temp_swing_k`; a tidally locked world has `temp_day_side_k`, `temp_terminator_k` and `temp_night_side_k` instead. A pure function of stored columns, so no epoch.
-- v138 (GEN.414, DB.68): `planets` and `moons` also keep the haze, dust and cloud optical depths at 550 nm (`tau_haze`, `tau_dust`, `tau_cloud`), the cloud cover and the number that picks
--   the aerosol optics inside their brackets (`aerosol_draw`). A row saved before v138 has NULL: the regenerated galaxy draws them.
-- v139 (MAP.143): `sector_stats.habitable_worlds` (the planets and moons of a generated sector's systems that are Blue or Green overall, NULL before a fill) and
--   `galaxy_block_stats.habitable_worlds` (their sum over a block's generated sectors), so the Galaxy Map can colour a sector or block by how many it holds.
-- v140 (DB.73): one parent link for every kind. `parent_type` (`star`, `planet`, `system`, `node`, `rogue_planet` or `free`) and
--   `parent_id` (the row of that table; for `node` the `orbit_nodes.node_no` of the row's own system; NULL for `system`'s
--   head when unknown and for `free`) are on planets and moons since v102 and now also on stars, asteroid_belts, comets, black_holes,
--   neutron_stars, rogue_planets and interstellar_comets, backfilled from `parent_node_id` and `star_id`, with an index on both
--   columns. `parent_type` is VARCHAR(16). `moons.planet_id` and `moons.star_system_id` are nullable: a moon of a free planet has
--   `parent_type` 'rogue_planet', the rogue planet's row in `parent_id` and in `rogue_planet_id` (a foreign key that cascades), and both NULL. A planet that leaves its star keeps its
--   row and ID and becomes `free`. There is no foreign key for a polymorphic `parent_id`. The `sector_objects` view lists the bodies of systems only, so a moon of a free planet is not in it.
-- v141 (GEN.459): `planets` and `moons` keep the radiogenic content of the body (`c_k0_ppm`, `c_th0_ppm`, `c_u0_ppm`: potassium, thorium and uranium in the bulk silicate at
--   formation, ppm by mass; `al26_ratio`: 26Al/27Al at formation; `urey_ratio`: the radiogenic share of the heat flowing out), the one set its interior heat and its ground dose read.
--   Rows from before v141 have NULL there and keep the heat and the dose they were generated with.
-- v145 (GEN.443): `stars.uv_c_share`, the share of its blackbody's UV-C a star gives off (the Sun's measured 0.256 is one value of the open bracket 0.1 to 1), drawn once from the star's own stream.
--   A row saved before v145 has NULL, which the radiation model reads as the blackbody (1).
-- v146 (MAP.258): `galaxy_block_stats.m` is INT UNSIGNED (the block ladder now runs from 177,147 down to 3, `galaxy/drill.py`). Every level below 243 changed shape with the levels put above it, so the table is filled again with a row for every level and the nebula block tables are made again.
-- v147 (MAP.259): `galaxy_block_lights` is new: the summed light of the pre-placed bright stars and scattered phenomena in each block of the ladder, filled by the migration from the stored scatter.
-- v148 (MAP.220): `habitable_systems`, one row per star system with a habitable world (its planets and moons that are Blue or Green overall),
--   with its count and position, so the Galaxy Map's marker layer finds the systems with one in a box. (`star_systems` has a FULLTEXT key,
--   so a column could only be added by copying the whole table.)
-- v149 (OPS.69): `sectors.wiki_url`, `star_systems.wikijs_url` and `star_systems.mediawiki_url` are gone with the wiki upload.
-- v151 (GEN.461): `planets` and `moons` also keep `water_base_pressure_pa`, the pressure under the water over the covered area, Pa (liquid, ice or both; the compressed column of
--   `physics/water_layer.py`). NULL for a gas giant, a body with no water layer and a row made before v151, whose `ocean_depth_km`, `ice_shell_km` and `hp_ice_km` are water-equivalent km.
-- v154 (DB.89, GEN.642): `planets` and `moons` gain the summary of a giant's interior (`interior_kind`, `core_mass_earth`, `dilute_core_mass_earth`,
--   `metallization_pressure_gpa`, `helium_depletion`, `has_helium_shell`), and `planet_interior_layers` and `moon_interior_layers` hold its layers, outside in: the
--   kind, the state, the pressure and radius fraction at the top and the bottom, and the mass of each material (Earth masses). NULL and no rows on a body without one.
-- v155 (DB.71, GEN.424): the air moves out of `planets` and `moons` into the child tables `planet_gases` and `moon_gases` (one row per gas above zero: `gas`, its formula in the gas table, and `partial_pa`,
--   its partial pressure in Pa); the ten `p_<gas>_kpa` columns are gone. A class's balance gas (helium, sodium, ammonia) is a row like the others, so the pressure is the sum of the rows. The migration
--   copies each row's ten partial pressures (kPa times 1000) into the child tables.
-- v159 (GEN.441): `planets` and `moons` gain the surface composition: `ground_type` and `ground_grain` and the areal share of each ground type (`ground_ocean_fraction`, `ground_ice_fraction`, `ground_ultramafic_fraction`, ..., twelve shares that add up to 1), NULL for a giant and for a row saved before v152. The minerals, ices and ocean the shares stand for are rows of the materials table (`physics/surface_composition.py`).
-- v180 (DB.70): `planets` and `moons` gain `moon_origin` (how a moon came to orbit its parent: `disk`, `impact`, `captured`, `ring` or `binary`; NULL for a planet and for a moon saved before) and a giant's circumplanetary disk (`moon_disk_lifetime_myr`, `moon_disk_alpha_over_f`, `moon_disk_f_inf`, `moon_disk_ice_line_km`; NULL for a body with no disk). Nothing is backfilled.
-- v164 (GEN.764): `scatter_epochs` holds one row per kind of scatter drawn (`mass`, `luminosity`, `band`, `phenomena`, `nebula-field`), the generator epoch and the run that drew it, and the row goes when the
--   scatter is cleared; `species` gains `generator_epoch`, the epoch the population pass made it under. A sector has its own epoch (v121) and the plan its own (v131), so every stage's results say what made them.
-- v165 (MAP.305, MAP.303): `system_interests` replaces `habitable_systems`: one row for each reason a star system is worth a mark on the Galaxy Map (a planet or moon rated
--   Yellow or better overall, a gas giant), with how many bodies give it, the largest score among them and the system's position.
-- v166 (MAP.321): `system_interests` also takes the reasons `uninhabitable` (a scored terrestrial planet or moon under the Yellow floor) and `asteroid-field` (a system with a belt); no
--   column changes, the migration fills the two reasons from the planets, moons and belts already stored.
-- v167 (MAP.321): `system_interests` also takes the reason `inhabited` (a world whose timeline reached a technological civilization, or with a terrestrial colony on it); no column
--   changes, the migration fills it from the planets, moons and colonies already stored.
-- v172 (MAP.220): `system_interests` also takes one reason for each planet class, `class-A` to `class-V` (a system with a planet or moon of it, `n` how many, no score); no column
--   changes, the migration fills them from the planets and moons already stored.
-- v173 (MAP.311): `system_interests` also takes the reasons `civilization` (a world whose timeline reached a technological civilization), `sentient` (a sentient life score of 1 or more) and `life` (any of the
--   four life scores of 1 or more), no score kept; no column changes, the migration fills them from the planets and moons already stored.
-- v174 (MAP.344): the `civilization` rows of `system_interests` keep, as `best_score`, the technology level of the system's most advanced species (NULL until its species is drawn); no column
--   changes, the migration sets it from the species already stored.
-- v175 (GEN.889, GEN.890): the life class switch. `planets` and `moons` lose `life_chemical` and gain `life_class` (the class number that lives on the body) and `life_class_variant` (its pigment or
--   habitat letter), `life_class_best` (the class id the body supports best) and `life_class_best_score` (0 to 9), `life_animal_grade` (0 none to 3 large and active) and `life_complexity` (the
--   complexity index of the class, `physics/life_complexity.py`); the classes a body supports are the child tables `planet_life_class_scores` and `moon_life_class_scores` (one row for each class
--   with a score of 1 or more); `species.life_chemical` becomes `species.life_class` (the id). A row made before v175 keeps the class its chemistry stood for and nothing else
--   (see `migrations/versions/0175_life_classes.py`); the rest comes when the galaxy is regenerated.
-- v176 (MAP.351): `system_interests` also takes the reasons `lifeclass-1a` to `lifeclass-9`, one for each life class that lives on a planet or moon of the system (`life_class` and
--   `life_class_variant`), no score kept; no column changes, the migration fills them from the planets and moons already stored.
-- v182 (GEN.848): `planets` and `moons` keep a day and night swing for every place that has a temperature: `temp_polar_swing_k`, `temp_north_swing_k`, `temp_south_swing_k` and
--   `temp_mean_swing_k` for a world that turns on its own (the equator's is `temp_swing_k`), `temp_day_side_swing_k`, `temp_terminator_swing_k` and `temp_night_side_swing_k` (0, a locked
--   world has no day) for a tidally locked one; kelvin of difference, peak to peak. NULL for a body with no regional temperatures and a row made before v182, which gets them when the galaxy
--   is regenerated.
-- v183 (GEN.834): `planets` and `moons` split `ice_mass_fraction` (every ice but water) by species: `ice_co2_mass_fraction`, `ice_ch4_mass_fraction`, `ice_nh3_mass_fraction`, `ice_n2_mass_fraction`,
--   `ice_so2_mass_fraction`, `ice_co_mass_fraction`, `ice_hcn_mass_fraction` and `ice_other_mass_fraction` (ethane, methanol, formaldehyde, hydrogen sulfide), shares of the whole body's mass that add up to
--   `ice_mass_fraction`. NULL for a row made before v183 and for a body with no composition; 0 for a body with no ice. Nothing is backfilled.
-- v184 (GEN.952): `scatter_members.radius_au` (DOUBLE, NULL), a member's extent in AU, so a page can work out the X-ray luminosity of an accreting pair; NULL on a row made before v184, nothing is backfilled.
-- v185 (GEN.846): `bright_stars.rotation_period_hours` (DOUBLE, NULL), the spin period in hours of a star the scatter placed, drawn from the star's own seed; NULL on a row made before v184, nothing is backfilled.
-- v187 (GEN.768): `planets` and `moons` keep the hours of one local day a person can be outside (`outside_hours_day`, DOUBLE) and the limit that sets them (`outside_hours_limit`, VARCHAR(16): daylight, dose,
--   heat, cold, air, uv or vehicle; `physics/outside_hours.py`). NULL for a gas giant, an unscored body, a world whose daylight is not stored and a row made before v187, which gets them when the galaxy is regenerated.
-- v189 (GEN.1009): `scatter_members.radius_max_au` and `mass_max_solar` (DOUBLE, NULL), the largest radius and mass a member has had (a white dwarf's red-giant tip and birth mass), so the history rule and the build read what the tree drew from; NULL on a row made before v187, nothing is backfilled.
-- v190 (GEN.474): `planets` and `moons` gain `ground_intermediate_fraction`, the areal share of intermediate (andesitic) ground, a ground type the volcanic crust of a body takes when the melt of its mantle
--   is 52 to 63 percent silica; NULL on a row made before v190 (the regenerated galaxy replaces it; see `migrations/versions/0190_ground_intermediate.py`). (v190, GEN.475) `asteroid_belts` and `asteroid_fields`
--   also gain `chondrite_group`, the nearest meteorite group of the body's mean composition (`none` far from every group), NULL before v190.
-- v193 (GEN.1097): `planets` and `moons` keep the state of their day and night (`day_night_state`, VARCHAR(24): day and night, tidally locked, spin-orbit resonance, non-rotating, polar day and night, no solid surface;
--   `physics/day_night.py`), the hours of night in a solar day at the equator (`night_hours`, DOUBLE) and the solar days in one orbit (`days_per_year`, DOUBLE). NULL for a body whose spin is not stored, `night_hours` and
--   `days_per_year` NULL for a world that turns once an orbit, and all NULL on a row made before v193, which gets them when the galaxy is regenerated.
-- v194 (GEN.1013): `dark_members`, one row for each member of a built system's tree that is not one of its stars (a brown dwarf, neutron star, black hole or stripped helium star): `parent_node_no`, `slot`, `kind`, `mass_solar`, `radius_au`. A system stored before has none; nothing is backfilled.
-- v196 (GEN.1108): `planets` and `moons` keep the climate of the body in seven bands of latitude (`physics/latitude_bands.py`, `physics/regional_temperature.py`), `band1` (the north pole) to `band7`
--   (the south pole), each with `bandN_mean_k` (the annual-mean temperature), `bandN_swing_k` (the day and night swing, peak to peak, kelvin of difference), `bandN_day_hours` and `bandN_night_hours`
--   (hours of day and of night in a solar day), `bandN_longest_day_days` and `bandN_longest_night_days` (the longest polar day and night of the year at the poleward edge) and `bandN_state`
--   (VARCHAR(24): day and night, polar day and night, non-rotating; for a world locked to its star the seven are rings from the day centre to the night centre, with no hours, and the state is
--   permanent day side, terminator band or permanent night side). `temp_mean_swing_k` stays (the swing of the whole surface). The columns `temp_polar_k`, `temp_north_k`, `temp_equator_k`, `temp_south_k`,
--   `temp_swing_k`, `temp_day_side_k`, `temp_terminator_k`, `temp_night_side_k` and the swings `temp_polar_swing_k`, `temp_north_swing_k`, `temp_south_swing_k`, `temp_day_side_swing_k`,
--   `temp_terminator_swing_k` and `temp_night_side_swing_k` are gone. NULL for a body with no temperature or rotation period and a row made before v194, which gets them when the galaxy is regenerated.
-- v197 (GEN.1174): `stars` and `bright_stars` gain `nuclear_arrival_age_gy` (Gyr, NULL), the age at which MIST's tracks put a pre-main-sequence star on the main sequence (hydrogen burning settled), kept beside `phase_end_age_gy` of `bright_stars`, which for such a star is the earlier age its radius settles (BHAC15). NULL for every star that is not pre-main-sequence and for a row made before v195; nothing is backfilled.
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
-- v197 (GEN.1215): `planets` and `moons` also keep what their surface liquids hold: `surface_liquid_acid_wt` and `surface_liquid_acid_layer` (the strength of the liquid sulfuric acid the air makes,
--   percent of its mass, and where it is, pool or cloud (`physics/acid_medium.py`)) and `surface_liquid_ammonia_wt` (the ammonia of a liquid water ocean, percent of the liquid's mass). NULL for a body without such a liquid and a row made before v197.


-- v39: star_systems and several phenomenon tables point at nebulae and
-- supernova_remnants, which are created further down, so foreign key
-- checks are off while this file runs (MySQL and MariaDB then accept a
-- reference to a table that doesn't exist yet) and back on at the end.
SET FOREIGN_KEY_CHECKS = 0;

-- ---------------------------------------------------------------------
-- schema_migrations -- DDL-level structure version tracking. Replaces
-- SQLite's `PRAGMA user_version` (see the MySQL port note above);
-- `store.py`'s `migrate_database` reads `MAX(version)` from this table
-- and inserts a row per migration step it applies.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS schema_migrations (
    version      INT NOT NULL PRIMARY KEY,
    applied_at   TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- alembic_version -- Alembic's own record of the revision a database is
-- at (DB.11), declared here as Alembic would create it so that a new
-- database and a migrated one have the same shape. Its single row is the
-- schema version as four digits ("0061"); `alembic_runner.upgrade` keeps
-- it in step with `schema_migrations`.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS alembic_version (
    version_num  VARCHAR(32) NOT NULL,
    CONSTRAINT alembic_version_pkc PRIMARY KEY (version_num)
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

    -- v30, re-keyed to the spherical grid in v91: the sector's cell (shell =
    -- radius band 4 pc deep, band = latitude ring, signed, 0 on the galactic
    -- plane, wedge = longitude slice of the band) -- see the header comment's
    -- "v91" note. Deliberately not part of the CHECK below.
    shell_index         INT,
    band_index          INT,
    wedge_index         INT,

    -- v27: row timestamps -- see the header comment's "v27" note.
    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BIGINT UNSIGNED,
    max_mass_kg                 DOUBLE NULL,  -- v102 (DB.55)
    max_mass_id                 BINARY(12) NULL,  -- v102 (DB.55)
    storage_checksum            BINARY(32) NULL,  -- v115 (DB.44)
    parity_dirty                TINYINT NOT NULL DEFAULT 1,  -- v115 (DB.44)
    leaf_digest                 BINARY(32) NULL,  -- v120 (DB.49)
    generator_epoch             SMALLINT UNSIGNED NULL,  -- v121 (DB.16)
    run_id                      BIGINT UNSIGNED NULL,  -- v121 (DB.16)
    density_multiplier_used     DECIMAL(10,4) NULL,  -- v125 (DB.60)
    override_id                 BIGINT UNSIGNED NULL,  -- v125 (DB.60)

    UNIQUE KEY uq_sectors_uid (uid),
    CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    -- At most one sector per address (v8, re-keyed in v30) --
    -- `generate.ensure_sector_generated` relies on this to turn a
    -- lazy-generation race into a clear IntegrityError rather than a
    -- duplicate sector. Never-placed sectors (NULL address) don't collide
    -- -- ordinary SQL NULL semantics for UNIQUE.
    UNIQUE KEY uq_sectors_address (shell_index, band_index, wedge_index),

    KEY idx_sectors_galactic_radius_pc (galactic_radius_pc),
    KEY idx_sectors_name (name),
    FULLTEXT KEY ft_sectors_name (name),
    -- v25: lets `queryDb.galaxy_sectors_in_view`'s bounding-box query
    -- range-scan on center_x_pc instead of a full table scan -- see the
    -- header comment's "v25" note.
    KEY idx_sectors_center (center_x_pc, center_y_pc, center_z_pc),
    -- v27: see the header comment's "v27" note.
    KEY idx_sectors_modified_at (modified_at),
    -- v117 (DB.16): lets the mixed-version check group sectors by run without a table scan.
    KEY idx_sectors_run (run_id),
    CONSTRAINT fk_sectors_run FOREIGN KEY (run_id) REFERENCES generation_runs(id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- galaxy_shape / galaxy_shell -- the galaxy-wide density "skeleton"
-- (v8, re-keyed to rings in v32, to layers in v33 and to shells in v91 --
-- see the header comment's notes).
-- `galaxy_shape` is a singleton (`id` pinned to 1, enforced by the CHECK)
-- -- there is exactly one galaxy. Building or rebuilding the skeleton
-- (`planetgen plan`) replaces this row and every `galaxy_shell` row
-- wholesale; neither
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
    arm_level                   DOUBLE NOT NULL DEFAULT 1,  -- v81 (ADM.49): see header comment
    nuclear_amplitude           DOUBLE NOT NULL DEFAULT 0,  -- v99 (GEN.280): see header comment
    spiral_reference_radius_pc  DOUBLE NOT NULL,
    spiral_reference_angle_rad  DOUBLE NOT NULL,
    k_norm                      DOUBLE NOT NULL,

    -- The sector edge length this skeleton was built at, in parsecs
    -- (galaxyGeometry's own scope -- every grid formula takes this as a
    -- parameter, as the shell depth, rather than assuming a fixed constant).
    edge_pc                     DOUBLE NOT NULL,

    -- SpaceSector(edge_ly=...).expected_system_count() at relative_density
    -- = 1 -- cached because every qualification check
    -- (predicted_star_count >= 1.0) needs it, and it's a fixed galaxy-wide
    -- fact, not worth re-deriving from edge_pc on every call.
    expected_system_count_at_density_1  DOUBLE NOT NULL,

    -- The last shell with any qualifying content (its equatorial band
    -- reaches farthest) -- this galaxy's real edge, found by `planetgen
    -- plan`, not picked as an arbitrary radius. (v32; `outer_ring_index`
    -- until v91.)
    outer_shell_index            INT NOT NULL,
    reach_shell_index            INT NOT NULL DEFAULT 65535,  -- v93 (GEN.233), see header comment
    -- v98 (GEN.281): a cell that expects more systems than `cell_star_ceiling` is filled only with the stars at or
    -- above `fill_luminosity_cap_sol` (`galaxy/resolution.py`); the rest is unresolved, computed and never stored.
    cell_star_ceiling            DOUBLE NOT NULL DEFAULT 5000,
    fill_luminosity_cap_sol      DOUBLE NOT NULL DEFAULT 100,
    -- v101 (DB.54): the scatter limits by place (docs/design/galaxy/realistic-galactic-core.md 5.1). The core's luminosity
    -- floor and mass limit are NULL until set, and then follow the baseline's times 1.1; the halo has a luminosity
    -- floor and no mass limit of its own.
    core_min_luminosity_sol      DOUBLE,
    core_mass_limit_sol          DOUBLE,
    halo_min_luminosity_sol      DOUBLE NOT NULL DEFAULT 500,
    -- v114 (DB.81): the stellar mass, solar masses, the plan was told to build the galaxy at (a pinned mass, or the one
    -- two pinned Scale levers give); NULL when the seed drew it. A loaded shape draws its mass chain and structure
    -- again at this mass; `galaxy_shape_pins` lists what the plan was told.
    pinned_stellar_mass_solar    DOUBLE,

    -- v117 (DB.66, docs/design/galaxy/galaxy-types.md section 10): the galaxy family and the numbers of the components a
    -- family may have that a spiral does not draw today. Nothing reads them yet (GEN.379 does). A galaxy planned
    -- before v117 is a spiral with no spheroid and no clump share. The bar and the thick disk were never stored, so
    -- NULL means "as the structure draw gives them from the seed"; an amplitude of 0 removes a component exactly.
    family                          VARCHAR(16) NOT NULL DEFAULT 'spiral',
    bar_present                     TINYINT,
    bar_axis_ratio_b                DOUBLE,
    bar_axis_ratio_c                DOUBLE,
    bar_angle_rad                   DOUBLE,
    thick_disk_amplitude            DOUBLE,
    spheroid_amplitude              DOUBLE NOT NULL DEFAULT 0,
    spheroid_sersic_n               DOUBLE NOT NULL DEFAULT 0,
    spheroid_effective_radius_pc    DOUBLE NOT NULL DEFAULT 0,
    spheroid_axis_ratio_b           DOUBLE NOT NULL DEFAULT 1,
    spheroid_axis_ratio_c           DOUBLE NOT NULL DEFAULT 1,
    spheroid_orientation_1_rad      DOUBLE NOT NULL DEFAULT 0,
    spheroid_orientation_2_rad      DOUBLE NOT NULL DEFAULT 0,
    spheroid_orientation_3_rad      DOUBLE NOT NULL DEFAULT 0,
    body_clump_share                DOUBLE NOT NULL DEFAULT 0,

    -- v43: the luminosity threshold (L_sun) and seed the bright-star
    -- scatter used (`bright_stars`); NULL when none has run. A fill reads
    -- this, not `program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL`.
    bright_star_min_luminosity_sol  DOUBLE,
    bright_star_seed                BIGINT UNSIGNED,
    -- v73 (GEN.185): the mass limit (solar masses) of the star scatter's mass
    -- pass: every star born at or above it is placed, so a fill and the
    -- luminosity pass draw only lighter ones. NULL when the scatter had no
    -- mass pass (a scatter drawn before v73).
    bright_star_mass_limit_sol      DOUBLE,

    -- v64 (GEN.100): the seed the phenomenon scatter (`phenomenon_scatter`)
    -- used; NULL when none has run (a sector's fill then rolls its own
    -- black holes, neutron stars, planetary nebulae and supernova remnants).
    phenomenon_scatter_seed         BIGINT UNSIGNED,
    -- v71 (GEN.167): the scatter's mass cut, see header comment.
    phenomenon_min_mass_solar       DOUBLE,

    -- v51: the galaxy's 128-bit seed (GEN.39), set by `planetgen plan`
    -- (`--seed`, or drawn at random) and kept by every later plan.
    -- NULL only on a galaxy planned before v51.
    galaxy_seed                     BINARY(16),

    -- v52 (DB.6): the code that wrote `galaxy_seed` -- the version key
    -- (MAJOR 4, REVISION 4, BUILD 6, Python 2+2+2, OS 1, architecture 1
    -- hex digits), the full release, Python version and platform.
    version_key                     CHAR(22),
    planetgen_version               VARCHAR(32),
    python_version                  VARCHAR(16),
    platform                        VARCHAR(64),

    -- v53 (PERF.11): a decaying average, over every sector fill, of the
    -- systems a sector got against the systems expected there
    -- (`sector_stats`), and how many fills went into it.
    density_ratio_avg               DOUBLE,
    density_ratio_samples           BIGINT UNSIGNED NOT NULL DEFAULT 0,

    -- v131 (OPS.58): what else the code that wrote `galaxy_seed` was -- the
    -- generator epoch (`util/generator_epoch.py`), the fingerprint format
    -- (`util/fingerprint_encoding.FP_SPEC`) and the generator battery's overall
    -- digest with its cases (`generation/generator_battery.py`, OPS.57), as
    -- `{case: full digest}`. NULL on a galaxy planned before v131.
    generator_epoch                 SMALLINT UNSIGNED NULL,
    fp_spec                         SMALLINT UNSIGNED NULL,
    battery_digest                  CHAR(64) NULL,
    battery_cases                   JSON NULL,

    CONSTRAINT chk_galaxy_shape_cell_star_ceiling CHECK (cell_star_ceiling >= 1),  -- v100 (GEN.281)
    CONSTRAINT chk_galaxy_shape_scatter_limits CHECK (  -- v101 (DB.54)
        halo_min_luminosity_sol >= 50 AND (core_min_luminosity_sol IS NULL OR core_min_luminosity_sol > 0)
        AND (core_mass_limit_sol IS NULL OR core_mass_limit_sol > 0)),
    CONSTRAINT chk_galaxy_shape_family CHECK (  -- v117 (DB.66)
        family IN ('spiral', 'lenticular', 'elliptical', 'irregular', 'dwarf_spheroid')),
    CONSTRAINT chk_galaxy_shape_bar_axis_ratio_b CHECK (bar_axis_ratio_b > 0 AND bar_axis_ratio_b <= 1),  -- v117 (DB.66), NULL passes
    CONSTRAINT chk_galaxy_shape_bar_axis_ratio_c CHECK (bar_axis_ratio_c > 0 AND bar_axis_ratio_c <= 1),
    CONSTRAINT chk_galaxy_shape_thick_disk_amplitude CHECK (thick_disk_amplitude >= 0),
    CONSTRAINT chk_galaxy_shape_spheroid_amplitude CHECK (spheroid_amplitude >= 0),
    CONSTRAINT chk_galaxy_shape_spheroid_sersic_n CHECK (spheroid_sersic_n >= 0),
    CONSTRAINT chk_galaxy_shape_spheroid_effective_radius CHECK (spheroid_effective_radius_pc >= 0),
    CONSTRAINT chk_galaxy_shape_spheroid_axis_ratio_b CHECK (spheroid_axis_ratio_b > 0 AND spheroid_axis_ratio_b <= 1),
    CONSTRAINT chk_galaxy_shape_spheroid_axis_ratio_c CHECK (spheroid_axis_ratio_c > 0 AND spheroid_axis_ratio_c <= 1),
    CONSTRAINT chk_galaxy_shape_body_clump_share CHECK (body_clump_share BETWEEN 0 AND 1)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- galaxy_shape_pins -- the shape levers a plan was told (v114, DB.81)
-- One row per lever the plan pinned (`planetgen.galaxy.dials`, the lever's name), with the value in the table's units
-- (lengths in parsecs). A lever with no row was left to the seed (Auto). Rebuilt with `galaxy_shape` by every plan;
-- the galaxy page and the API say "mass pinned at X, the rest drawn" from it.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS galaxy_shape_pins (
    name                        VARCHAR(40) NOT NULL PRIMARY KEY,
    value                       DOUBLE NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- One row per shell that holds content (v91; replaces `galaxy_layer` and
-- `galaxy_column`): the galaxy as a stack of nested spheres, from the core
-- (shell 0) out to `outer_shell_index`. Each shell holds bands
-- -`high_band_index` through `high_band_index` (both polar caps too when
-- `high_band_index` is the shell's top band), the widest run around the plane
-- that holds a sector center clearing the qualification threshold for some
-- longitude (`galaxySkeleton.build_shell_extents`). Always symmetric about the
-- plane. `cell_count` is the number of sectors that run holds, worked out once
-- when the outline is built. A shell with no row has no content at all.
CREATE TABLE IF NOT EXISTS galaxy_shell (
    shell_index       INT NOT NULL PRIMARY KEY,
    high_band_index   INT NOT NULL,
    cell_count        BIGINT UNSIGNED NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- galaxy_block_stats (v92, MAP.205): one row per block `(m, shell_index,
-- band_index, wedge_index)` of a level of the Galaxy Map's block ladder above the
-- sector (`galaxy.drill`, docs/design/galaxy/spherical-block-ladder.md: blocks `m` shells
-- deep, rows of blocks by colatitude, wedges by longitude) that holds at least one
-- generated sector, with the sums of what those sectors hold: `sector_count`,
-- `systems`, `expected_systems`, `stars`, `age_sum` (each sector's mean age times its
-- stars) and `luminosity_sol`. A stage's child stats are these rows (a level-3
-- block's children are its sectors, read from `sector_stats`). Written by
-- `db/block_stats.py` when sectors are saved or deleted, and by `planetgen plan
-- --fill-block-stats`; no foreign key, as the blocks are not rows of their own.
-- v146: `m` is INT UNSIGNED: the ladder's root is 177,147 shells deep.
CREATE TABLE IF NOT EXISTS galaxy_block_stats (
    m                 INT UNSIGNED NOT NULL,
    shell_index       INT NOT NULL,
    band_index        INT NOT NULL,
    wedge_index       INT NOT NULL,
    sector_count      INT UNSIGNED NOT NULL,
    systems           BIGINT UNSIGNED NOT NULL,
    expected_systems  DOUBLE NOT NULL,
    stars             BIGINT UNSIGNED NOT NULL,
    age_sum           DOUBLE NOT NULL,
    luminosity_sol    DOUBLE NOT NULL,
    habitable_worlds  BIGINT UNSIGNED NOT NULL DEFAULT 0,   -- v139 (MAP.143)
    mass_sum_kg       DOUBLE NOT NULL DEFAULT 0,   -- v168 (DB.108): the summed `sector_stats.total_mass_kg` of its generated sectors
    PRIMARY KEY (m, shell_index, band_index, wedge_index)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- galaxy_block_lights (v147, MAP.259): the light of the objects pre-placed in a block, for the map's
-- sector-block nodes. One row per block `(m, shell_index, band_index, wedge_index)` of a level above a
-- sector that holds at least one bright star no system was built around (`bright_stars.star_system_id IS
-- NULL`) or scattered phenomenon not yet built (`phenomenon_scatter.built_at IS NULL`): their count and
-- summed luminosity (Lsun), luminosity times position (`*_lpx`, `*_lpy`, `*_lpz`, Lsun times parsecs, so
-- the centre of light is each over the luminosity) and, for stars, luminosity times temperature
-- (`star_lt`). Built objects are in `galaxy_block_stats`, so no light is counted twice. Written by
-- `db/block_lights.py` when sectors are saved or deleted and when a scatter is drawn or cleared, and by
-- `planetgen fill-block-totals`; no foreign key, as the blocks are not rows of their own.
CREATE TABLE IF NOT EXISTS galaxy_block_lights (
    m                    INT UNSIGNED NOT NULL,
    shell_index          INT NOT NULL,
    band_index           INT NOT NULL,
    wedge_index          INT NOT NULL,
    star_count           BIGINT UNSIGNED NOT NULL,
    star_luminosity_sol  DOUBLE NOT NULL,
    star_lpx              DOUBLE NOT NULL,
    star_lpy              DOUBLE NOT NULL,
    star_lpz              DOUBLE NOT NULL,
    star_lt              DOUBLE NOT NULL,
    point_count          BIGINT UNSIGNED NOT NULL,
    point_luminosity_sol DOUBLE NOT NULL,
    point_lpx             DOUBLE NOT NULL,
    point_lpy             DOUBLE NOT NULL,
    point_lpz             DOUBLE NOT NULL,
    PRIMARY KEY (m, shell_index, band_index, wedge_index)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- region_stats (v161, DB.107): one row per named region of the Galaxy Map (`galaxy/named_regions.py`): a 9-block, `(shell_index,
-- band_index, wedge_index)` (the level, 9, is implied), that holds at least one generated sector. `sector_count` is how many of
-- its sectors are generated and `capacity` how many its outline holds (the share charted is `sector_count / capacity`),
-- `clump_count` how many face-connected clumps those sectors make, and `extent_x_pc`, `extent_y_pc`, `extent_z_pc` and
-- `extent_radius_pc` the bounding sphere of the largest clump (never smaller than the region's framing sphere), galaxy
-- frame, parsecs: what the map frames and draws the region with. A cell with no generated sector has no row, and the
-- region does not exist. Written by `db/region_stats.py`, which `db/block_stats.py` calls whenever it puts right the
-- blocks above a sector (a save, a deletion, a rebuild); no foreign key, as the cell is not a row of its own.
CREATE TABLE IF NOT EXISTS region_stats (
    shell_index       INT NOT NULL,
    band_index        INT NOT NULL,
    wedge_index       INT NOT NULL,
    sector_count      INT UNSIGNED NOT NULL,
    capacity          INT UNSIGNED NOT NULL,
    clump_count       INT UNSIGNED NOT NULL,
    extent_x_pc       DOUBLE NOT NULL,
    extent_y_pc       DOUBLE NOT NULL,
    extent_z_pc       DOUBLE NOT NULL,
    extent_radius_pc  DOUBLE NOT NULL,
    mass_sum_kg       DOUBLE NOT NULL DEFAULT 0,   -- v168 (DB.108): the summed `sector_stats.total_mass_kg` of its generated sectors
    PRIMARY KEY (shell_index, band_index, wedge_index)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Singleton row (same pattern as galaxy_shape above) tracking when
-- planetgen.cli.orbits last advanced every planet's/moon's orbital_phase_deg in
-- this database, so the next run knows how much real time has actually
-- elapsed since then. Absent entirely until planetgen.cli.orbits's first run
-- against a given database (it creates this row itself).
CREATE TABLE IF NOT EXISTS orbit_simulation_state (
    id               BIGINT UNSIGNED PRIMARY KEY CHECK (id = 1),
    last_updated_at  TIMESTAMP NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- system_configs -- one row per SystemConfig "recipe"
-- (planetgen/generation/config.py:29-33 SERIALIZABLE_FIELDS)
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS system_configs (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    habitable_world   TINYINT(1) CHECK (habitable_world IN (0, 1)),
    asteroid_belt     TINYINT(1) CHECK (asteroid_belt IN (0, 1)),
    comets            TINYINT(1) CHECK (comets IN (0, 1)),          -- v50
    large_star        TINYINT(1) CHECK (large_star IN (0, 1)),
    moons             TINYINT(1) CHECK (moons IN (0, 1)),
    max_planets       TINYINT(1) CHECK (max_planets IN (0, 1)),
    planets           TINYINT(1) CHECK (planets IN (0, 1)),
    star_type         VARCHAR(64),
    name              VARCHAR(255),
    age               VARCHAR(16) CHECK (age IN ('young', 'old')),
    intelligent_life  TINYINT(1) CHECK (intelligent_life IN (0, 1)),
    binary_system     TINYINT(1) CHECK (binary_system IN (0, 1)),
    wide_binary       TINYINT(1) CHECK (wide_binary IN (0, 1)),     -- v50
    num_orbits        INT,
    prevalence_habitable_world DOUBLE NOT NULL DEFAULT 0,  -- v54
    prevalence_asteroid_belt DOUBLE NOT NULL DEFAULT 0,
    prevalence_comets DOUBLE NOT NULL DEFAULT 0,
    prevalence_large_star DOUBLE NOT NULL DEFAULT 0,
    prevalence_moons DOUBLE NOT NULL DEFAULT 0,
    prevalence_max_planets DOUBLE NOT NULL DEFAULT 0,
    prevalence_intelligent_life DOUBLE NOT NULL DEFAULT 0,
    prevalence_binary_system DOUBLE NOT NULL DEFAULT 0,
    prevalence_wide_binary DOUBLE NOT NULL DEFAULT 0,
    prevalence_planets DOUBLE NOT NULL DEFAULT 0,
    stars             TINYINT UNSIGNED,  -- v134 (GEN.267): the forced number of stars, NULL when left to chance
    max_stars         TINYINT UNSIGNED NOT NULL DEFAULT 12,  -- v134: the cap on the star count
    prevalence_multi_star DOUBLE NOT NULL DEFAULT 0  -- v134: the chance of five or more stars, percent from usual
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
    -- v93: the system's point in the galaxy frame, parsecs: its sector's center plus the position above, turned to the
    -- galaxy's axes (`frame.local_to_galaxy_pc`). NULL iff the position or the sector's center is NULL.
    galaxy_x_pc          DOUBLE,
    galaxy_y_pc          DOUBLE,
    galaxy_z_pc          DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(galaxy_x_pc) < 262144 AND ABS(galaxy_y_pc) < 262144 AND ABS(galaxy_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((galaxy_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((galaxy_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((galaxy_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((galaxy_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((galaxy_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((galaxy_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((galaxy_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((galaxy_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((galaxy_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    KEY idx_star_systems_cell (cell_key, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc),

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
    binary_epoch_unix                          DOUBLE,  -- v66 (GEN.106)
    binary_next_update_due                     DOUBLE,
    binary_mutual_min_update_interval_years    DOUBLE,
    -- v14 (see header comment): the secondary's position relative to the
    -- primary, kept in lockstep with binary_mutual_orbital_phase_deg --
    -- same "position relative to whatever this orbit is around" convention
    -- as planets/moons' own position_x/y/z_km.
    binary_mutual_position_x_km   DOUBLE,
    binary_mutual_position_y_km   DOUBLE,
    binary_mutual_position_z_km   DOUBLE,
    -- v20 (see header comment): each binary member's own offset from the
    -- pair's barycenter -- secondary_position = primary_position +
    -- binary_mutual_position above always holds, but both are stored
    -- explicitly rather than derived on read.
    binary_primary_position_x_km    DOUBLE,
    binary_primary_position_y_km    DOUBLE,
    binary_primary_position_z_km    DOUBLE,
    binary_secondary_position_x_km  DOUBLE,
    binary_secondary_position_y_km  DOUBLE,
    binary_secondary_position_z_km  DOUBLE,
    -- v20: secondary_mass / (primary_mass + secondary_mass), constant --
    -- stored so advance_orbital_phases never needs to join back to `stars`.
    binary_secondary_mass_fraction  DOUBLE,
    -- v20: additional pair-wide wobble from circumbinary (P-type) planets
    -- (star_id IS NULL) -- NULL unless binary_configuration = 'close'.
    binary_planetary_wobble_x_km    DOUBLE,
    binary_planetary_wobble_y_km    DOUBLE,
    binary_planetary_wobble_z_km    DOUBLE,

    system_flavor_text   TEXT,
    -- v37: 'runaway' (> 30 km/s) or 'hypervelocity' (> 500 km/s) when the
    -- system moves unusually fast (generate.flag_fast_stars); NULL
    -- together otherwise.
    runaway_class        VARCHAR(16),
    runaway_speed_kms    DOUBLE,
    -- v61: velocity on the galactic axes, km/s (see the header comment's "v61" note).
    velocity_x_kms       DOUBLE NOT NULL DEFAULT 0,
    velocity_y_kms       DOUBLE NOT NULL DEFAULT 0,
    velocity_z_kms       DOUBLE NOT NULL DEFAULT 0,
    epoch_unix           DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due      DOUBLE,
    schema_version       INT NOT NULL DEFAULT 1,

    -- No stored page text since v29: Markdown and plain text are rendered on
    -- demand from the rows above -- see the header comment's "v29" note.

    -- v39: the innermost nebula or supernova remnant this sits inside
    -- (at most one set; _db.refresh_containment).
    inside_nebula_id    BIGINT UNSIGNED,
    inside_remnant_id   BIGINT UNSIGNED,

    created_at           TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at          TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),

    UNIQUE KEY uq_star_systems_uid (uid),
    CONSTRAINT fk_star_systems_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE SET NULL,
    CONSTRAINT fk_star_systems_config
        FOREIGN KEY (system_config_id) REFERENCES system_configs(id),
    KEY idx_star_systems_sector_id (sector_id),
    KEY idx_star_systems_next_update_due (next_update_due),
    KEY idx_star_systems_binary_next_update_due (binary_next_update_due),
    KEY idx_star_systems_system_config_id (system_config_id),
    KEY idx_star_systems_name (name),
    FULLTEXT KEY ft_star_systems_name (name),
    -- v27: see the header comment's "v27" note.
    KEY idx_star_systems_modified_at (modified_at),
    -- v63 (DB.13): the Systems table's octant filter and sort.
    -- v82 (PERF.70): (quadrant, name, id) serves the sort by Octant a group at a time.
    KEY idx_star_systems_quadrant (quadrant, name, id),
    -- v82 (PERF.70): the same for the sort by Binary.
    KEY idx_star_systems_binary_name (is_binary, name, id),
    -- v39: containment -- see the header comment's "v39" note.
    CONSTRAINT fk_star_systems_inside_nebula
        FOREIGN KEY (inside_nebula_id) REFERENCES nebulae(id) ON DELETE SET NULL,
    CONSTRAINT fk_star_systems_inside_remnant
        FOREIGN KEY (inside_remnant_id) REFERENCES supernova_remnants(id) ON DELETE SET NULL
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
    role                      VARCHAR(16) NOT NULL CHECK (role IN ('primary', 'secondary', 'single', 'member')),  -- v133 (GEN.630): 'member'
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
    galactic_orbital_period_gy    DOUBLE NOT NULL,
    galactic_orbital_phase_deg          DOUBLE NOT NULL,  -- v13, see header comment
    galactic_min_update_interval_years  DOUBLE NOT NULL,  -- v13, see header comment
    -- v15: this star's own Holman & Wiegert (1999) critical semi-major axis
    -- -- NULL except for a constituent of a 'wide' binary. See the header
    -- comment's "v15" note.
    wide_binary_a_crit_km               DOUBLE,
    -- v20 (see header comment): this star's own reflex-offset "wobble"
    -- from the combined pull of every planet orbiting it directly --
    -- NULL/0 with no planets. Class-blind: applies identically to an
    -- anchored black_holes/neutron_stars row.
    reflex_offset_x_km       DOUBLE,
    reflex_offset_y_km       DOUBLE,
    reflex_offset_z_km       DOUBLE,
    -- v68 (GEN.104): the spin, see header comment.
    rotation_period_hours       DOUBLE,
    spin_axis_x                 DOUBLE,
    spin_axis_y                 DOUBLE,
    spin_axis_z                 DOUBLE,
    axial_tilt_deg              DOUBLE,
    -- v72 (GEN.86): stellar activity, see header comment.
    log_lx_lbol                 DOUBLE,
    l_xuv_w                     DOUBLE,
    xuv_saturated               BOOLEAN,
    flare_n33_per_yr            DOUBLE,
    flare_alpha                 DOUBLE,
    xuv_fluence_j               DOUBLE,
    -- v75 (GEN.87): supernovae within 8 pc per Gyr here, see header comment.
    lethal_event_rate_per_gyr   DOUBLE,
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),
    -- v95 (GEN.158): the star's [Fe/H] in dex against the Sun; NULL for a row saved before v95.
    metallicity_fe_h      DOUBLE,
    parent_node_id              SMALLINT NULL,  -- v104 (DB.53)
    slot                        TINYINT UNSIGNED NULL,  -- v133 (GEN.630): 0 or 1, the child of its parent node this star is
    state                       VARCHAR(16) NULL,  -- v104 (DB.53)
    flags                       SMALLINT UNSIGNED NULL,  -- v104 (DB.53)
    offset_x_au                 DOUBLE NULL,  -- v104 (DB.53)
    offset_y_au                 DOUBLE NULL,  -- v104 (DB.53)
    offset_z_au                 DOUBLE NULL,  -- v104 (DB.53)
    velocity_x_kms              DOUBLE NULL,  -- v104 (DB.53)
    velocity_y_kms              DOUBLE NULL,  -- v104 (DB.53)
    velocity_z_kms              DOUBLE NULL,  -- v104 (DB.53)
    decay_age_yr                DOUBLE NULL,  -- v104 (DB.53)
    nuclear_arrival_age_gy      DOUBLE NULL,  -- v197 (GEN.1174): the age MIST's tracks reach the main sequence, see the header comment
    epoch_unix                  DOUBLE NULL,  -- v104 (DB.53)
    disk_kind                   VARCHAR(16) NULL,  -- v104 (DB.53)
    disk_mass_solar             DOUBLE NULL,  -- v104 (DB.53)
    disk_inner_au               DOUBLE NULL,  -- v104 (DB.53)
    disk_outer_au               DOUBLE NULL,  -- v104 (DB.53)
    cluster_id                  BIGINT UNSIGNED NULL,  -- v127 (DB.61)
    snow_line_km                DOUBLE NOT NULL,  -- v132 (GEN.458): this star's own snow line, drawn once around the 170 K blackbody position
    uv_c_share                  DOUBLE NULL,  -- v145 (GEN.443): the UV-C share of the star's blackbody, drawn once (log-uniform 0.1 to 1), NULL for a remnant and for a row saved before v145

    parent_type                 VARCHAR(16) NOT NULL DEFAULT 'system',  -- v140 (DB.73): what it orbits, see the v140 note
    parent_id                   BIGINT UNSIGNED NULL,  -- v140 (DB.73)

    UNIQUE KEY uq_stars_uid (uid),
    KEY idx_stars_parent (parent_type, parent_id),
    CONSTRAINT chk_stars_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),
    CONSTRAINT fk_stars_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    KEY idx_stars_star_system_id (star_system_id),
    KEY idx_stars_name (name),
    FULLTEXT KEY ft_stars_name (name),
    KEY idx_stars_yerkes_class (yerkes_class),
    KEY idx_stars_radius_km (radius_km),
    KEY idx_stars_star_type (star_type),
    KEY idx_stars_cluster (cluster_id),
    CONSTRAINT fk_stars_cluster
        FOREIGN KEY (cluster_id) REFERENCES galaxy_components(id) ON DELETE SET NULL
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
    flavor_text_count         INT NOT NULL DEFAULT 0,
    life_stage                VARCHAR(32),      -- v110 (DB.50): a MILESTONE_KEYS key, NULL before any milestone
    evolutionary_pace         VARCHAR(16),      -- v110: the pace life evolves at (fast, normal, slow), drawn once, set with any life chemistry (v143)
    life_milestone_age_gyr    DOUBLE,           -- v110: the age of that milestone, Gy
    life_system_age_gyr       DOUBLE,           -- v110: the age of the system the timeline was written at, Gy
    life_age_abiogenesis_gyr  DOUBLE,           -- v156 (DB.100): the age the pace reaches abiogenesis, Gy
    life_age_photosynthesis_gyr DOUBLE,           -- v156 (DB.100): the age the pace reaches photosynthesis, Gy
    life_age_complex_cells_gyr DOUBLE,           -- v156 (DB.100): the age the pace reaches complex cells, Gy
    life_age_multicellularity_gyr DOUBLE,           -- v156 (DB.100): the age the pace reaches multicellularity, Gy
    life_age_technological_civilization_gyr DOUBLE,           -- v156 (DB.100): the age the pace reaches technological civilization, Gy
    flavor_key                VARCHAR(32),      -- v110: which flavour line (replaces flavor_text)
    -- Orbital motion (v9, see header comment) -- inclination/ascending
    -- node are fixed at generation time. rotation_period_hours is a
    -- separate, static "day length" stat. There is no stored phase (v109,
    -- GEN.289): the body is its state vector below.
    orbital_inclination_deg     DOUBLE NOT NULL,
    orbital_ascending_node_deg  DOUBLE NOT NULL,
    eccentricity                DOUBLE NOT NULL DEFAULT 0,   -- v159 (GEN.727): of the orbit the body was generated on (distance_km is its semi-major axis)
    -- Position/speed (v11, see header comment) -- Cartesian position
    -- relative to this planet's orbital anchor (the star, or the
    -- BinaryStarProxy's combined center for a binary system), made from
    -- distance_km and the orbital-motion columns above and moved by
    -- planetgen.cli.orbits with kepler.universal_step, with the velocity
    -- below and epoch_unix as its state.
    -- orbital_speed_kms is constant around a circular orbit -- it only
    -- changes if distance_km/period_years do.
    position_x_km               DOUBLE NOT NULL,
    position_y_km               DOUBLE NOT NULL,
    position_z_km               DOUBLE NOT NULL,
    orbital_speed_kms           DOUBLE NOT NULL,
    -- v60: velocity relative to the orbital anchor, km/s, in step with the position (see the header comment's "v60" note).
    velocity_x_kms              DOUBLE NOT NULL DEFAULT 0,
    velocity_y_kms              DOUBLE NOT NULL DEFAULT 0,
    velocity_z_kms              DOUBLE NOT NULL DEFAULT 0,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,
    min_update_interval_years   DOUBLE NOT NULL,  -- v12, see header comment
    rotation_period_hours       DOUBLE NOT NULL,
    -- v68 (GEN.104): the spin, see header comment.
    spin_axis_x                 DOUBLE,
    spin_axis_y                 DOUBLE,
    spin_axis_z                 DOUBLE,
    axial_tilt_deg              DOUBLE,
    -- v70 (GEN.85): mantle redox, see header comment (the air's partial pressures went to `planet_gases` and `moon_gases` in v155).
    mantle_redox                VARCHAR(16) CHECK (mantle_redox IN ('reduced', 'intermediate', 'oxidized')),
    mantle_delta_iw             DOUBLE,
    -- v72 (GEN.86): magnetic field and the star's XUV and flares here, see header comment.
    magnetic_moment_a_m2        DOUBLE,
    dipole_class                VARCHAR(16) CHECK (dipole_class IN ('none', 'weak', 'earth-like', 'strong', 'multipolar')),
    induced_moment_a_m2         DOUBLE,
    magnetopause_rp             DOUBLE,
    xuv_flux_earth              DOUBLE,
    xuv_exposure_index          DOUBLE,
    flare_irradiation_index     DOUBLE,
    -- v74 (GEN.88): water, land, ice and the ocean's chemistry, see header comment.
    water_mass_fraction         DOUBLE,
    hydrosphere                 VARCHAR(20) CHECK (hydrosphere IN ('dry', 'vapour', 'ice', 'ice-covered ocean', 'surface ocean', 'brine seep', 'hycean')),
    ocean_fraction              DOUBLE,
    land_fraction               DOUBLE,
    ocean_depth_km              DOUBLE,
    ice_shell_km                DOUBLE,
    hp_ice_km                   DOUBLE,
    ocean_class                 VARCHAR(16) CHECK (ocean_class IN ('ice-sealed', 'chloride brine', 'acid sulfate', 'soda', 'neutral')),
    ocean_ph                    DOUBLE,
    water_activity              DOUBLE,
    phosphorus                  VARCHAR(8) CHECK (phosphorus IN ('high', 'limited', 'starved')),
    -- v75 (GEN.87): surface dose, UV and ozone, see header comment.
    surface_dose_msv_yr         DOUBLE,
    dose_gcr_msv_yr             DOUBLE,
    dose_sep_msv_yr             DOUBLE,
    dose_ground_msv_yr          DOUBLE,
    dose_belt_msv_yr            DOUBLE,  -- v187 (GEN.983): the dose a moon takes from its parent's trapped-particle belts, mSv/yr (0 for a planet)
    dose_helio_mult             DOUBLE,
    uv_surface_index            DOUBLE,
    ozone_loss_flag             BOOLEAN,
    -- v76 (GEN.89): the habitability score, see header comment.
    phi4                        DOUBLE,
    phi4_pressure               DOUBLE,
    phi4_temperature            DOUBLE,
    phi4_air                    DOUBLE,
    phi4_water                  DOUBLE,
    phi4_suitability            DOUBLE,
    phi4_radiation              DOUBLE,
    tier_pressure               TINYINT,
    tier_temperature            TINYINT,
    tier_air                    TINYINT,
    tier_water                  TINYINT,
    tier_suitability            TINYINT,
    tier_radiation              TINYINT,
    life_microbial              TINYINT UNSIGNED NULL,  -- v150 (GEN.649): 0 none to 9 abundant, from the evolutionary path and the surface
    life_flora                  TINYINT UNSIGNED NULL,  -- v150
    life_fauna                  TINYINT UNSIGNED NULL,  -- v150
    life_sentient               TINYINT UNSIGNED NULL,  -- v150
    l_solv                      DOUBLE,
    l_chem                      DOUBLE,
    l_ener                      DOUBLE,
    l_rad                       DOUBLE,
    eq_body                     TINYINT,
    eq_body_kind                TINYINT,
    eq_body_rec                 TINYINT,
    eq_body_rec_kind            TINYINT,
    eq_rad_suit                 TINYINT,
    eq_rad_suit_rec             TINYINT,
    eq_mobility                 TINYINT,
    eq_dust                     TINYINT,
    eq_gas_filter               TINYINT,
    eq_co2_scrub                TINYINT,
    eq_supply                   TINYINT,
    eq_supply_rec               TINYINT,
    eq_shelter                  TINYINT,
    eq_shelter_rec              TINYINT,
    shelter_column_g_cm2        DOUBLE,
    outside_hours               DOUBLE,
    band1_mean_k                DOUBLE,  -- v196 (GEN.1108): the climate of the body in seven bands of latitude (`physics/latitude_bands.py`): `band1` is the north pole, `band7` the south pole (a world locked to its star: the day centre to the night centre)
    band1_swing_k               DOUBLE,
    band1_day_hours             DOUBLE,
    band1_night_hours           DOUBLE,
    band1_longest_day_days      DOUBLE,
    band1_longest_night_days    DOUBLE,
    band1_state                 VARCHAR(24),
    band2_mean_k                DOUBLE,
    band2_swing_k               DOUBLE,
    band2_day_hours             DOUBLE,
    band2_night_hours           DOUBLE,
    band2_longest_day_days      DOUBLE,
    band2_longest_night_days    DOUBLE,
    band2_state                 VARCHAR(24),
    band3_mean_k                DOUBLE,
    band3_swing_k               DOUBLE,
    band3_day_hours             DOUBLE,
    band3_night_hours           DOUBLE,
    band3_longest_day_days      DOUBLE,
    band3_longest_night_days    DOUBLE,
    band3_state                 VARCHAR(24),
    band4_mean_k                DOUBLE,
    band4_swing_k               DOUBLE,
    band4_day_hours             DOUBLE,
    band4_night_hours           DOUBLE,
    band4_longest_day_days      DOUBLE,
    band4_longest_night_days    DOUBLE,
    band4_state                 VARCHAR(24),
    band5_mean_k                DOUBLE,
    band5_swing_k               DOUBLE,
    band5_day_hours             DOUBLE,
    band5_night_hours           DOUBLE,
    band5_longest_day_days      DOUBLE,
    band5_longest_night_days    DOUBLE,
    band5_state                 VARCHAR(24),
    band6_mean_k                DOUBLE,
    band6_swing_k               DOUBLE,
    band6_day_hours             DOUBLE,
    band6_night_hours           DOUBLE,
    band6_longest_day_days      DOUBLE,
    band6_longest_night_days    DOUBLE,
    band6_state                 VARCHAR(24),
    band7_mean_k                DOUBLE,
    band7_swing_k               DOUBLE,
    band7_day_hours             DOUBLE,
    band7_night_hours           DOUBLE,
    band7_longest_day_days      DOUBLE,
    band7_longest_night_days    DOUBLE,
    band7_state                 VARCHAR(24),
    temp_mean_swing_k           DOUBLE,  -- v182 (GEN.848): the day and night swing of the whole surface, peak to peak, kelvin of difference
    outside_hours_day           DOUBLE,  -- v187 (GEN.768): the hours of one local day a person can be outside in the listed equipment
    outside_hours_limit         VARCHAR(16),  -- v187 (GEN.768): what sets them: daylight, dose, heat, cold, air, uv or vehicle
    day_night_state             VARCHAR(24),  -- v188 (GEN.1097): the state of its day and night: day and night, tidally locked, spin-orbit resonance, non-rotating, polar day and night or no solid surface
    night_hours                 DOUBLE,  -- v188 (GEN.1097): hours of night in one solar day at the equator, NULL for a world that turns once an orbit
    days_per_year               DOUBLE,  -- v188 (GEN.1097): solar days in one orbit round the star, NULL for a world that turns once an orbit
    energy_flux_w_m2            DOUBLE,
    flux_variation              DOUBLE,
    -- v20 (see header comment): this planet's own reflex-offset "wobble"
    -- from the combined pull of its own moons -- NULL/0 with no moons.
    reflex_offset_x_km       DOUBLE,
    reflex_offset_y_km       DOUBLE,
    reflex_offset_z_km       DOUBLE,
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),
    parent_type                 VARCHAR(16) NOT NULL DEFAULT 'star',  -- v102 (DB.55), v140 (DB.73): wider
    parent_id                   BIGINT UNSIGNED NULL,  -- v102 (DB.55)
    dynamical_state             VARCHAR(8) NOT NULL DEFAULT 'bound',  -- v102 (DB.55)
    host_node_id                SMALLINT NULL,  -- v104 (DB.53): the orbit node the planet orbits (NULL: its star)
    born_beyond_snow_line       BOOLEAN NULL,  -- v132 (GEN.458): formed beyond its star's snow line, NULL when the host has none to compare
    ozone_column_du             DOUBLE NULL,  -- v135 (GEN.413): the ozone column, Dobson units
    tau_haze                    DOUBLE NULL,  -- v138 (GEN.414, DB.68): haze optical depth at 550 nm, NULL for a row saved before v138 or a body with no air
    tau_dust                    DOUBLE NULL,  -- v138 (GEN.414, DB.68): dust optical depth at 550 nm
    tau_cloud                   DOUBLE NULL,  -- v138 (GEN.414, DB.68): in-cloud optical depth at 550 nm (grey)
    cloud_cover                 DOUBLE NULL,  -- v138 (GEN.414, DB.68): share of the sky under cloud, 0 to 1
    aerosol_draw                DOUBLE NULL,  -- v138 (GEN.414, DB.68): the number in [0, 1) that picks the optics inside the brackets
    c_k0_ppm                    DOUBLE,  -- v141 (GEN.459): potassium in the bulk silicate at formation, ppm by mass
    c_th0_ppm                   DOUBLE,  -- v141: thorium, the same
    c_u0_ppm                    DOUBLE,  -- v141: uranium, the same
    al26_ratio                  DOUBLE,  -- v141: 26Al/27Al at formation
    urey_ratio                  DOUBLE,  -- v141: the radiogenic share of the heat flowing out, today (an output of the thermal history since v146)
    bond_albedo                 DOUBLE NULL,  -- v142 (GEN.412): the Bond albedo drawn for the temperature
    greenhouse_factor           DOUBLE NULL,  -- v142: the greenhouse quantity G (T_s = T_e (1 + G)^(1/4)), 0 for no air
    albedo_basis                TINYINT UNSIGNED NULL,  -- v142: 0 recorded at generation, 1 exact inverse of a stored row, 2 seeded class draw
    metal_mass_fraction         DOUBLE,  -- v144 (DB.69, GEN.428): the iron-rich core's share of the whole body's mass, light elements included
    rock_mass_fraction          DOUBLE,  -- v144: the silicate (or silicon carbide and graphite) share
    ice_mass_fraction           DOUBLE,  -- v144: every ice other than water, as one share
    ice_co2_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is carbon dioxide ice
    ice_ch4_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is methane ice
    ice_nh3_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is ammonia ice
    ice_n2_mass_fraction        DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is nitrogen ice
    ice_so2_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is sulfur dioxide ice
    ice_co_mass_fraction        DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is carbon monoxide ice
    ice_hcn_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is hydrogen cyanide ice
    ice_other_mass_fraction     DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is every other ice (ethane, methanol, formaldehyde, hydrogen sulfide) ice
    envelope_mass_fraction      DOUBLE,  -- v144: the hydrogen and helium share
    rock_feo_fraction           DOUBLE,  -- v144: FeO, a mass fraction of the rock
    rock_al2o3_fraction         DOUBLE,  -- v144: Al2O3, a mass fraction of the rock
    rock_mg_si_ratio            DOUBLE,  -- v144: the rock's molar Mg/Si
    envelope_helium_fraction    DOUBLE,  -- v144: Y, a mass fraction of the envelope
    envelope_z_fraction         DOUBLE,  -- v144: heavy elements, a mass fraction of the envelope
    porosity                    DOUBLE,  -- v144: the void share of a small body's volume
    carbon_rock                 BOOLEAN,  -- v144: the rock is silicon carbide and graphite, not silicate
    composition_origin          VARCHAR(8) CHECK (composition_origin IN ('drawn', 'inverse')),  -- v144: drawn (composition first) or inverse (worked back from a mass and a radius)
    ground_type                   VARCHAR(16),  -- v159 (GEN.441): the ground type with the largest areal share (ocean, ice or a rock type)
    ground_grain                  VARCHAR(8),  -- v159: slab, crushed or powder
    ground_ocean_fraction         DOUBLE,  -- v159: the areal share of the surface that is ocean (the twelve shares add up to 1)
    ground_ice_fraction           DOUBLE,  -- v159: the areal share of the surface that is ice
    ground_ultramafic_fraction    DOUBLE,  -- v159: the areal share of the surface that is ultramafic
    ground_anorthositic_fraction  DOUBLE,  -- v159: the areal share of the surface that is anorthositic
    ground_basaltic_fraction      DOUBLE,  -- v159: the areal share of the surface that is basaltic
    ground_felsic_fraction        DOUBLE,  -- v159: the areal share of the surface that is felsic
    ground_clay_fraction          DOUBLE,  -- v159: the areal share of the surface that is clay
    ground_oxidised_fraction      DOUBLE,  -- v159: the areal share of the surface that is oxidised
    ground_metal_fraction         DOUBLE,  -- v159: the areal share of the surface that is metal
    ground_carbon_fraction        DOUBLE,  -- v159: the areal share of the surface that is carbon
    ground_evaporite_fraction     DOUBLE,  -- v159: the areal share of the surface that is evaporite
    ground_sulfur_fraction        DOUBLE,  -- v159: the areal share of the surface that is sulfur
    ground_intermediate_fraction  DOUBLE,  -- v190 (GEN.474): the areal share of the surface that is intermediate (andesitic) ground
    surface_liquid                VARCHAR(32),  -- v188 (GEN.1035): the liquid the surface holds (an id of physics/surface_liquids.py, NULL: none)
    surface_liquid_fraction       DOUBLE,  -- v188: the share of the surface that liquid covers
    surface_liquid_depth_km       DOUBLE,  -- v188: its mean depth over that share, km (NULL: not known)
    surface_liquid_acid_wt        DOUBLE,  -- v197 (GEN.1215): the strength of the liquid sulfuric acid the air makes, percent of its mass (NULL: none)
    surface_liquid_acid_layer     VARCHAR(8),  -- v197: where that acid is, pool (at the ground) or cloud
    surface_liquid_ammonia_wt     DOUBLE,  -- v197: the ammonia of a liquid water ocean, percent of the liquid's mass (NULL: no ocean or no ammonia)
    moon_origin                   VARCHAR(8),  -- v180 (DB.70): how a moon came to orbit its parent (NULL for a planet and for a moon saved before)
    moon_disk_lifetime_myr        DOUBLE,  -- v180: a giant planet's circumplanetary disk: its lifetime, Myr (NULL for a body with no disk)
    moon_disk_alpha_over_f        DOUBLE,  -- v180: the disk's alpha over f (the viscosity over the supply ratio of Canup and Ward)
    moon_disk_f_inf               DOUBLE,  -- v180: the asymptotic ice fraction of the disk's ramp
    moon_disk_ice_line_km         DOUBLE,  -- v180: the distance of the disk's ice line from the giant, km
    water_base_pressure_pa      DOUBLE NULL,  -- v151 (GEN.461): the pressure at the bottom of the water layer, Pa
    tidal_flux_w_m2             DOUBLE NULL,  -- v152 (GEN.700): the eccentricity-tide heat flow of an icy moon under ice, W/m2 (NULL for a planet)
    aquifer_depth_km            DOUBLE NULL,  -- v157 (GEN.654, DB.92): the depth to the top of the water-bearing layer underground, km (NULL: none)
    aquifer_thickness_km        DOUBLE NULL,  -- v156: the thickness of that layer, km
    mantle_eta0_pa_s            DOUBLE NULL,  -- v153 (GEN.699): the viscosity prefactor of the mantle, eta = eta0 exp(Ea / R T), Pa s
    mantle_activation_kj_mol    DOUBLE NULL,  -- v153: its viscosity activation energy Ea, kJ/mol
    internal_heat_flux_w_m2     DOUBLE NULL,  -- v153: the heat flowing out of the interior, radiogenic plus primordial, over the body's surface
    interior_kind               VARCHAR(16) CHECK (interior_kind IN ('gas_giant', 'ice_giant', 'brown_dwarf')),  -- v154 (DB.89, GEN.642): the kind of interior a giant was drawn with
    core_mass_earth             DOUBLE,  -- v154: the compact core's (an ice giant's rock core's) mass in Earth masses, 0 with none
    dilute_core_mass_earth      DOUBLE,  -- v154: the dilute core's mass in Earth masses, 0 with none
    metallization_pressure_gpa  DOUBLE,  -- v154: the drawn molecular-to-metallic pressure, GPa
    helium_depletion            DOUBLE,  -- v154: the molecular envelope's helium share against the natal share, 1 with no helium rain
    has_helium_shell            BOOLEAN,  -- v154: a helium-rich shell above the metallic envelope
    ring_kind                   VARCHAR(8) CHECK (ring_kind IN ('dusty', 'narrow', 'broad')),  -- v160 (GEN.732): the kind of the ring system, NULL with none
    ring_origin                 VARCHAR(16) CHECK (ring_origin IN ('moon_erosion', 'disrupted_moon', 'captured_body')),  -- v160: where its material came from
    ring_inner_km               DOUBLE,  -- v160: the inner edge of the ring system, km from the body's centre
    ring_outer_km               DOUBLE,  -- v160: its outer edge, km (inside the Roche limit of its own material and the innermost moon)
    ring_mass_kg                DOUBLE,  -- v160: its mass, kg
    ring_optical_depth          DOUBLE,  -- v160: its mean optical depth over the area it covers
    ring_ice_fraction           DOUBLE,  -- v160: its water-ice share by mass
    ring_age_yr                 DOUBLE,  -- v192 (GEN.1057): years since the event that made it (a fallen moon or a passing comet), NULL for a resonance ring
    atmosphere_colour           INT UNSIGNED,  -- v162 (DB.96): the colour of the body's air, sRGB as red*65536 + green*256 + blue, NULL with no air
    lid_regime                  VARCHAR(10) NULL,  -- v169 (GEN.721): mobile (plate tectonics) or stagnant (one rigid plate), drawn once, NULL for a gas giant
    day_share_equator           DOUBLE NULL,  -- v170 (GEN.767): share of the time the star is up at the equator over a year, 0 to 1, NULL for a world that turns once an orbit
    day_share_mean              DOUBLE NULL,  -- v170 (GEN.767): share of the time the star is up, averaged over the surface, 0 to 1
    parent_shadow_share         DOUBLE NULL,  -- v170 (GEN.767): a moon only, the share of its time in its planet's shadow, 0 to 1
    solar_day_hours             DOUBLE NULL,  -- v170 (GEN.767): hours from one noon to the next, NULL for a world that turns once an orbit
    daylight_hours              DOUBLE NULL,  -- v170 (GEN.767): hours of daylight in one solar day at the equator, NULL for a world that turns once an orbit
    polar_circle_deg            DOUBLE NULL,  -- v170 (GEN.767): latitude beyond which the star stays up or down for at least a day, degrees
    life_class                  TINYINT UNSIGNED NULL,  -- v175 (GEN.889): the number of the life class that lives on the body (physics/life_classes.py), NULL for none
    life_class_variant          CHAR(1) NULL,  -- v175: its pigment or habitat letter (1a, 1b, 1c, 3a ...), NULL for a class with none
    life_class_best             VARCHAR(2) NULL,  -- v175: the id of the best-supported class, NULL when no class is supported
    life_class_best_score       TINYINT UNSIGNED NULL,  -- v175: its score, 0 to 9
    life_animal_grade           TINYINT UNSIGNED NULL,  -- v175 (GEN.890): 0 no animals, 1 sessile, 2 small and mobile, 3 large and active
    life_complexity             DOUBLE NULL,  -- v175: the complexity index of its life class, 1 for an Earth-like biosphere
    tide_range_spring_m         DOUBLE NULL,  -- v191 (GEN.1034): open-ocean equilibrium spring tide range, metres, 0 under an ice cover, NULL without a sea
    tide_range_neap_m           DOUBLE NULL,  -- v191 (GEN.1034): the neap range, metres
    ocean_liquid                VARCHAR(40) NULL,  -- v191 (GEN.1034): the liquid of the sea the tide is of (water until GEN.1035), NULL without a sea

    UNIQUE KEY uq_planets_uid (uid),
    CONSTRAINT fk_planets_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_planets_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,
    KEY idx_planets_star_system_id (star_system_id),
    KEY idx_planets_next_update_due (next_update_due),
    KEY idx_planets_star_id (star_id),
    KEY idx_planets_name (name),
    FULLTEXT KEY ft_planets_name (name),
    KEY idx_planets_planet_class (planet_class),
    KEY idx_planets_body_type (body_type),
    KEY idx_planets_life_class (life_class),
    KEY idx_planets_radius_km (radius_km),
    CONSTRAINT chk_planets_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),
    KEY idx_planets_parent (parent_type, parent_id),
    CONSTRAINT chk_planets_dynamical_state CHECK (dynamical_state IN ('bound', 'unbound'))
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

-- A giant's layers, outside in (DB.89, GEN.642): one row per layer of the interior drawn for the planet, in Earth masses per material.
CREATE TABLE IF NOT EXISTS planet_interior_layers (
    id                      BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    planet_id               BIGINT UNSIGNED NOT NULL,
    layer_index             SMALLINT UNSIGNED NOT NULL,
    kind                    VARCHAR(24) NOT NULL CHECK (kind IN ('molecular_envelope', 'metallic_envelope', 'helium_gradient', 'helium_shell', 'dilute_core', 'compact_core', 'ice_mantle', 'rock_core')),
    state                   VARCHAR(16) NOT NULL CHECK (state IN ('molecular_fluid', 'metallic_fluid', 'ionic_fluid', 'superionic', 'molten')),
    top_pressure_gpa        DOUBLE NOT NULL,
    bottom_pressure_gpa     DOUBLE NOT NULL,
    mass_earth              DOUBLE NOT NULL,
    top_radius_fraction     DOUBLE NOT NULL,
    bottom_radius_fraction  DOUBLE NOT NULL,
    hydrogen                DOUBLE NOT NULL,
    helium                  DOUBLE NOT NULL,
    iron                    DOUBLE NOT NULL,
    rock                    DOUBLE NOT NULL,
    water_liquid            DOUBLE NOT NULL,
    ice_x                   DOUBLE NOT NULL,
    nh3_ice                 DOUBLE NOT NULL,
    ch4_ice                 DOUBLE NOT NULL,

    UNIQUE KEY uq_planet_interior_layers_index (planet_id, layer_index),
    CONSTRAINT fk_planet_interior_layers_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- The air of a planet (v155, DB.71, GEN.424): one row per gas above zero, its formula in the gas table (`physics/data/gases.json`) and its partial pressure.
-- A class's balance gas (helium, sodium, ammonia) is a row like the others. No row means no air (or a row made before GEN.85).
CREATE TABLE IF NOT EXISTS planet_gases (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    planet_id  BIGINT UNSIGNED NOT NULL,
    gas         VARCHAR(16) NOT NULL,
    partial_pa  DOUBLE NOT NULL CHECK (partial_pa > 0),

    CONSTRAINT fk_planet_gases_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    UNIQUE KEY uq_planet_gases (planet_id, gas)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- The life classes a planet supports (v175, GEN.889): one row for each class of `physics/life_classes.py` whose score for the body is 1 or more (0 to 9). Indexed on the class and the score so the
-- map filter and the search facet are cheap.
CREATE TABLE IF NOT EXISTS planet_life_class_scores (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    planet_id   BIGINT UNSIGNED NOT NULL,
    life_class  VARCHAR(2) NOT NULL,
    score       TINYINT UNSIGNED NOT NULL CHECK (score BETWEEN 1 AND 9),

    CONSTRAINT fk_planet_life_class_scores_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    UNIQUE KEY uq_planet_life_class_scores (planet_id, life_class),
    KEY idx_planet_life_class_scores_class (life_class, score)
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
    planet_id                 BIGINT UNSIGNED NULL,  -- v140 (DB.73): NULL for a moon of a free planet
    star_system_id            BIGINT UNSIGNED NULL,  -- v140 (DB.73): NULL for a moon of a free planet
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
    flavor_text_count         INT NOT NULL DEFAULT 0,
    life_stage                VARCHAR(32),      -- v110 (DB.50): a MILESTONE_KEYS key, NULL before any milestone
    evolutionary_pace         VARCHAR(16),      -- v110: the pace life evolves at (fast, normal, slow), drawn once, set with any life chemistry (v143)
    life_milestone_age_gyr    DOUBLE,           -- v110: the age of that milestone, Gy
    life_system_age_gyr       DOUBLE,           -- v110: the age of the system the timeline was written at, Gy
    life_age_abiogenesis_gyr  DOUBLE,           -- v156 (DB.100): the age the pace reaches abiogenesis, Gy
    life_age_photosynthesis_gyr DOUBLE,           -- v156 (DB.100): the age the pace reaches photosynthesis, Gy
    life_age_complex_cells_gyr DOUBLE,           -- v156 (DB.100): the age the pace reaches complex cells, Gy
    life_age_multicellularity_gyr DOUBLE,           -- v156 (DB.100): the age the pace reaches multicellularity, Gy
    life_age_technological_civilization_gyr DOUBLE,           -- v156 (DB.100): the age the pace reaches technological civilization, Gy
    flavor_key                VARCHAR(32),      -- v110: which flavour line (replaces flavor_text)
    -- Orbital motion (v9, see header comment) -- inclination/ascending
    -- node are fixed at generation time. rotation_period_hours is a
    -- separate, static "day length" stat. There is no stored phase (v109,
    -- GEN.289): the body is its state vector below.
    orbital_inclination_deg     DOUBLE NOT NULL,
    orbital_ascending_node_deg  DOUBLE NOT NULL,
    eccentricity                DOUBLE NOT NULL DEFAULT 0,   -- v159 (GEN.727): of the orbit the body was generated on (distance_km is its semi-major axis)
    -- Position/speed (v11, see header comment) -- Cartesian position
    -- relative to this moon's orbital anchor, its parent planet; made
    -- from distance_km and the orbital-motion columns above the same way
    -- planets.position_x/y/z_km are (see that table's comment).
    position_x_km               DOUBLE NOT NULL,
    position_y_km               DOUBLE NOT NULL,
    position_z_km               DOUBLE NOT NULL,
    orbital_speed_kms           DOUBLE NOT NULL,
    -- v60: velocity relative to the orbital anchor, km/s, in step with the position (see the header comment's "v60" note).
    velocity_x_kms              DOUBLE NOT NULL DEFAULT 0,
    velocity_y_kms              DOUBLE NOT NULL DEFAULT 0,
    velocity_z_kms              DOUBLE NOT NULL DEFAULT 0,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,
    min_update_interval_years   DOUBLE NOT NULL,  -- v12, see header comment
    rotation_period_hours       DOUBLE NOT NULL,
    -- v68 (GEN.104): the spin, see header comment.
    spin_axis_x                 DOUBLE,
    spin_axis_y                 DOUBLE,
    spin_axis_z                 DOUBLE,
    axial_tilt_deg              DOUBLE,
    -- v70 (GEN.85): mantle redox, see header comment (the air's partial pressures went to `planet_gases` and `moon_gases` in v155).
    mantle_redox                VARCHAR(16) CHECK (mantle_redox IN ('reduced', 'intermediate', 'oxidized')),
    mantle_delta_iw             DOUBLE,
    -- v72 (GEN.86): magnetic field and the star's XUV and flares here, see header comment.
    magnetic_moment_a_m2        DOUBLE,
    dipole_class                VARCHAR(16) CHECK (dipole_class IN ('none', 'weak', 'earth-like', 'strong', 'multipolar')),
    induced_moment_a_m2         DOUBLE,
    magnetopause_rp             DOUBLE,
    xuv_flux_earth              DOUBLE,
    xuv_exposure_index          DOUBLE,
    flare_irradiation_index     DOUBLE,
    -- v74 (GEN.88): water, land, ice and the ocean's chemistry, see header comment.
    water_mass_fraction         DOUBLE,
    hydrosphere                 VARCHAR(20) CHECK (hydrosphere IN ('dry', 'vapour', 'ice', 'ice-covered ocean', 'surface ocean', 'brine seep', 'hycean')),
    ocean_fraction              DOUBLE,
    land_fraction               DOUBLE,
    ocean_depth_km              DOUBLE,
    ice_shell_km                DOUBLE,
    hp_ice_km                   DOUBLE,
    ocean_class                 VARCHAR(16) CHECK (ocean_class IN ('ice-sealed', 'chloride brine', 'acid sulfate', 'soda', 'neutral')),
    ocean_ph                    DOUBLE,
    water_activity              DOUBLE,
    phosphorus                  VARCHAR(8) CHECK (phosphorus IN ('high', 'limited', 'starved')),
    -- v75 (GEN.87): surface dose, UV and ozone, see header comment.
    surface_dose_msv_yr         DOUBLE,
    dose_gcr_msv_yr             DOUBLE,
    dose_sep_msv_yr             DOUBLE,
    dose_ground_msv_yr          DOUBLE,
    dose_belt_msv_yr            DOUBLE,  -- v187 (GEN.983): the dose a moon takes from its parent's trapped-particle belts, mSv/yr (0 for a planet)
    dose_helio_mult             DOUBLE,
    uv_surface_index            DOUBLE,
    ozone_loss_flag             BOOLEAN,
    -- v76 (GEN.89): the habitability score, see header comment.
    phi4                        DOUBLE,
    phi4_pressure               DOUBLE,
    phi4_temperature            DOUBLE,
    phi4_air                    DOUBLE,
    phi4_water                  DOUBLE,
    phi4_suitability            DOUBLE,
    phi4_radiation              DOUBLE,
    tier_pressure               TINYINT,
    tier_temperature            TINYINT,
    tier_air                    TINYINT,
    tier_water                  TINYINT,
    tier_suitability            TINYINT,
    tier_radiation              TINYINT,
    life_microbial              TINYINT UNSIGNED NULL,  -- v150 (GEN.649): 0 none to 9 abundant, from the evolutionary path and the surface
    life_flora                  TINYINT UNSIGNED NULL,  -- v150
    life_fauna                  TINYINT UNSIGNED NULL,  -- v150
    life_sentient               TINYINT UNSIGNED NULL,  -- v150
    l_solv                      DOUBLE,
    l_chem                      DOUBLE,
    l_ener                      DOUBLE,
    l_rad                       DOUBLE,
    eq_body                     TINYINT,
    eq_body_kind                TINYINT,
    eq_body_rec                 TINYINT,
    eq_body_rec_kind            TINYINT,
    eq_rad_suit                 TINYINT,
    eq_rad_suit_rec             TINYINT,
    eq_mobility                 TINYINT,
    eq_dust                     TINYINT,
    eq_gas_filter               TINYINT,
    eq_co2_scrub                TINYINT,
    eq_supply                   TINYINT,
    eq_supply_rec               TINYINT,
    eq_shelter                  TINYINT,
    eq_shelter_rec              TINYINT,
    shelter_column_g_cm2        DOUBLE,
    outside_hours               DOUBLE,
    band1_mean_k                DOUBLE,  -- v196 (GEN.1108): the climate of the body in seven bands of latitude (`physics/latitude_bands.py`): `band1` is the north pole, `band7` the south pole (a world locked to its star: the day centre to the night centre)
    band1_swing_k               DOUBLE,
    band1_day_hours             DOUBLE,
    band1_night_hours           DOUBLE,
    band1_longest_day_days      DOUBLE,
    band1_longest_night_days    DOUBLE,
    band1_state                 VARCHAR(24),
    band2_mean_k                DOUBLE,
    band2_swing_k               DOUBLE,
    band2_day_hours             DOUBLE,
    band2_night_hours           DOUBLE,
    band2_longest_day_days      DOUBLE,
    band2_longest_night_days    DOUBLE,
    band2_state                 VARCHAR(24),
    band3_mean_k                DOUBLE,
    band3_swing_k               DOUBLE,
    band3_day_hours             DOUBLE,
    band3_night_hours           DOUBLE,
    band3_longest_day_days      DOUBLE,
    band3_longest_night_days    DOUBLE,
    band3_state                 VARCHAR(24),
    band4_mean_k                DOUBLE,
    band4_swing_k               DOUBLE,
    band4_day_hours             DOUBLE,
    band4_night_hours           DOUBLE,
    band4_longest_day_days      DOUBLE,
    band4_longest_night_days    DOUBLE,
    band4_state                 VARCHAR(24),
    band5_mean_k                DOUBLE,
    band5_swing_k               DOUBLE,
    band5_day_hours             DOUBLE,
    band5_night_hours           DOUBLE,
    band5_longest_day_days      DOUBLE,
    band5_longest_night_days    DOUBLE,
    band5_state                 VARCHAR(24),
    band6_mean_k                DOUBLE,
    band6_swing_k               DOUBLE,
    band6_day_hours             DOUBLE,
    band6_night_hours           DOUBLE,
    band6_longest_day_days      DOUBLE,
    band6_longest_night_days    DOUBLE,
    band6_state                 VARCHAR(24),
    band7_mean_k                DOUBLE,
    band7_swing_k               DOUBLE,
    band7_day_hours             DOUBLE,
    band7_night_hours           DOUBLE,
    band7_longest_day_days      DOUBLE,
    band7_longest_night_days    DOUBLE,
    band7_state                 VARCHAR(24),
    temp_mean_swing_k           DOUBLE,  -- v182 (GEN.848): the day and night swing of the whole surface, peak to peak, kelvin of difference
    outside_hours_day           DOUBLE,  -- v187 (GEN.768): the hours of one local day a person can be outside in the listed equipment
    outside_hours_limit         VARCHAR(16),  -- v187 (GEN.768): what sets them: daylight, dose, heat, cold, air, uv or vehicle
    day_night_state             VARCHAR(24),  -- v188 (GEN.1097): the state of its day and night: day and night, tidally locked, spin-orbit resonance, non-rotating, polar day and night or no solid surface
    night_hours                 DOUBLE,  -- v188 (GEN.1097): hours of night in one solar day at the equator, NULL for a world that turns once an orbit
    days_per_year               DOUBLE,  -- v188 (GEN.1097): solar days in one orbit round the star, NULL for a world that turns once an orbit
    energy_flux_w_m2            DOUBLE,
    flux_variation              DOUBLE,
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),
    parent_type                 VARCHAR(16) NOT NULL DEFAULT 'planet',  -- v102 (DB.55), v140 (DB.73): wider
    parent_id                   BIGINT UNSIGNED NULL,  -- v102 (DB.55)
    dynamical_state             VARCHAR(8) NOT NULL DEFAULT 'bound',  -- v102 (DB.55)
    born_beyond_snow_line       BOOLEAN NULL,  -- v132 (GEN.458): the parent planet's flag
    ozone_column_du             DOUBLE NULL,  -- v135 (GEN.413): the ozone column, Dobson units
    tau_haze                    DOUBLE NULL,  -- v138 (GEN.414, DB.68): haze optical depth at 550 nm, NULL for a row saved before v138 or a body with no air
    tau_dust                    DOUBLE NULL,  -- v138 (GEN.414, DB.68): dust optical depth at 550 nm
    tau_cloud                   DOUBLE NULL,  -- v138 (GEN.414, DB.68): in-cloud optical depth at 550 nm (grey)
    cloud_cover                 DOUBLE NULL,  -- v138 (GEN.414, DB.68): share of the sky under cloud, 0 to 1
    aerosol_draw                DOUBLE NULL,  -- v138 (GEN.414, DB.68): the number in [0, 1) that picks the optics inside the brackets
    rogue_planet_id             BIGINT UNSIGNED NULL,  -- v140 (DB.73): the rogue planet a moon of a free planet orbits
    c_k0_ppm                    DOUBLE,  -- v141 (GEN.459): potassium in the bulk silicate at formation, ppm by mass
    c_th0_ppm                   DOUBLE,  -- v141: thorium, the same
    c_u0_ppm                    DOUBLE,  -- v141: uranium, the same
    al26_ratio                  DOUBLE,  -- v141: 26Al/27Al at formation
    urey_ratio                  DOUBLE,  -- v141: the radiogenic share of the heat flowing out, today (an output of the thermal history since v146)
    bond_albedo                 DOUBLE NULL,  -- v142 (GEN.412): the Bond albedo drawn for the temperature
    greenhouse_factor           DOUBLE NULL,  -- v142: the greenhouse quantity G (T_s = T_e (1 + G)^(1/4)), 0 for no air
    albedo_basis                TINYINT UNSIGNED NULL,  -- v142: 0 recorded at generation, 1 exact inverse of a stored row, 2 seeded class draw
    metal_mass_fraction         DOUBLE,  -- v144 (DB.69, GEN.428): the iron-rich core's share of the whole body's mass, light elements included
    rock_mass_fraction          DOUBLE,  -- v144: the silicate (or silicon carbide and graphite) share
    ice_mass_fraction           DOUBLE,  -- v144: every ice other than water, as one share
    ice_co2_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is carbon dioxide ice
    ice_ch4_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is methane ice
    ice_nh3_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is ammonia ice
    ice_n2_mass_fraction        DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is nitrogen ice
    ice_so2_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is sulfur dioxide ice
    ice_co_mass_fraction        DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is carbon monoxide ice
    ice_hcn_mass_fraction       DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is hydrogen cyanide ice
    ice_other_mass_fraction     DOUBLE,  -- v183 (GEN.834): the share of the whole body's mass that is every other ice (ethane, methanol, formaldehyde, hydrogen sulfide) ice
    envelope_mass_fraction      DOUBLE,  -- v144: the hydrogen and helium share
    rock_feo_fraction           DOUBLE,  -- v144: FeO, a mass fraction of the rock
    rock_al2o3_fraction         DOUBLE,  -- v144: Al2O3, a mass fraction of the rock
    rock_mg_si_ratio            DOUBLE,  -- v144: the rock's molar Mg/Si
    envelope_helium_fraction    DOUBLE,  -- v144: Y, a mass fraction of the envelope
    envelope_z_fraction         DOUBLE,  -- v144: heavy elements, a mass fraction of the envelope
    porosity                    DOUBLE,  -- v144: the void share of a small body's volume
    carbon_rock                 BOOLEAN,  -- v144: the rock is silicon carbide and graphite, not silicate
    composition_origin          VARCHAR(8) CHECK (composition_origin IN ('drawn', 'inverse')),  -- v144: drawn (composition first) or inverse (worked back from a mass and a radius)
    ground_type                   VARCHAR(16),  -- v159 (GEN.441): the ground type with the largest areal share (ocean, ice or a rock type)
    ground_grain                  VARCHAR(8),  -- v159: slab, crushed or powder
    ground_ocean_fraction         DOUBLE,  -- v159: the areal share of the surface that is ocean (the twelve shares add up to 1)
    ground_ice_fraction           DOUBLE,  -- v159: the areal share of the surface that is ice
    ground_ultramafic_fraction    DOUBLE,  -- v159: the areal share of the surface that is ultramafic
    ground_anorthositic_fraction  DOUBLE,  -- v159: the areal share of the surface that is anorthositic
    ground_basaltic_fraction      DOUBLE,  -- v159: the areal share of the surface that is basaltic
    ground_felsic_fraction        DOUBLE,  -- v159: the areal share of the surface that is felsic
    ground_clay_fraction          DOUBLE,  -- v159: the areal share of the surface that is clay
    ground_oxidised_fraction      DOUBLE,  -- v159: the areal share of the surface that is oxidised
    ground_metal_fraction         DOUBLE,  -- v159: the areal share of the surface that is metal
    ground_carbon_fraction        DOUBLE,  -- v159: the areal share of the surface that is carbon
    ground_evaporite_fraction     DOUBLE,  -- v159: the areal share of the surface that is evaporite
    ground_sulfur_fraction        DOUBLE,  -- v159: the areal share of the surface that is sulfur
    ground_intermediate_fraction  DOUBLE,  -- v190 (GEN.474): the areal share of the surface that is intermediate (andesitic) ground
    surface_liquid                VARCHAR(32),  -- v188 (GEN.1035): the liquid the surface holds (an id of physics/surface_liquids.py, NULL: none)
    surface_liquid_fraction       DOUBLE,  -- v188: the share of the surface that liquid covers
    surface_liquid_depth_km       DOUBLE,  -- v188: its mean depth over that share, km (NULL: not known)
    surface_liquid_acid_wt        DOUBLE,  -- v197 (GEN.1215): the strength of the liquid sulfuric acid the air makes, percent of its mass (NULL: none)
    surface_liquid_acid_layer     VARCHAR(8),  -- v197: where that acid is, pool (at the ground) or cloud
    surface_liquid_ammonia_wt     DOUBLE,  -- v197: the ammonia of a liquid water ocean, percent of the liquid's mass (NULL: no ocean or no ammonia)
    moon_origin                   VARCHAR(8),  -- v180 (DB.70): how a moon came to orbit its parent (NULL for a planet and for a moon saved before)
    moon_disk_lifetime_myr        DOUBLE,  -- v180: a giant planet's circumplanetary disk: its lifetime, Myr (NULL for a body with no disk)
    moon_disk_alpha_over_f        DOUBLE,  -- v180: the disk's alpha over f (the viscosity over the supply ratio of Canup and Ward)
    moon_disk_f_inf               DOUBLE,  -- v180: the asymptotic ice fraction of the disk's ramp
    moon_disk_ice_line_km         DOUBLE,  -- v180: the distance of the disk's ice line from the giant, km
    water_base_pressure_pa      DOUBLE NULL,  -- v147 (GEN.461): the pressure at the bottom of the water layer, Pa
    tidal_flux_w_m2             DOUBLE NULL,  -- v152 (GEN.700): the eccentricity-tide heat flow of an icy moon under ice, W/m2
    aquifer_depth_km            DOUBLE NULL,  -- v157 (GEN.654, DB.92): the depth to the top of the water-bearing layer underground, km (NULL: none)
    aquifer_thickness_km        DOUBLE NULL,  -- v156: the thickness of that layer, km
    mantle_eta0_pa_s            DOUBLE NULL,  -- v153 (GEN.699): the viscosity prefactor of the mantle, eta = eta0 exp(Ea / R T), Pa s
    mantle_activation_kj_mol    DOUBLE NULL,  -- v153: its viscosity activation energy Ea, kJ/mol
    internal_heat_flux_w_m2     DOUBLE NULL,  -- v153: the heat flowing out of the interior, radiogenic plus primordial, over the body's surface
    interior_kind               VARCHAR(16) CHECK (interior_kind IN ('gas_giant', 'ice_giant', 'brown_dwarf')),  -- v154 (DB.89, GEN.642): the kind of interior a giant was drawn with
    core_mass_earth             DOUBLE,  -- v154: the compact core's (an ice giant's rock core's) mass in Earth masses, 0 with none
    dilute_core_mass_earth      DOUBLE,  -- v154: the dilute core's mass in Earth masses, 0 with none
    metallization_pressure_gpa  DOUBLE,  -- v154: the drawn molecular-to-metallic pressure, GPa
    helium_depletion            DOUBLE,  -- v154: the molecular envelope's helium share against the natal share, 1 with no helium rain
    has_helium_shell            BOOLEAN,  -- v154: a helium-rich shell above the metallic envelope
    ring_kind                   VARCHAR(8) CHECK (ring_kind IN ('dusty', 'narrow', 'broad')),  -- v160 (GEN.732): the kind of the ring system, NULL with none
    ring_origin                 VARCHAR(16) CHECK (ring_origin IN ('moon_erosion', 'disrupted_moon', 'captured_body')),  -- v160: where its material came from
    ring_inner_km               DOUBLE,  -- v160: the inner edge of the ring system, km from the body's centre
    ring_outer_km               DOUBLE,  -- v160: its outer edge, km (inside the Roche limit of its own material and the innermost moon)
    ring_mass_kg                DOUBLE,  -- v160: its mass, kg
    ring_optical_depth          DOUBLE,  -- v160: its mean optical depth over the area it covers
    ring_ice_fraction           DOUBLE,  -- v160: its water-ice share by mass
    ring_age_yr                 DOUBLE,  -- v192 (GEN.1057): years since the event that made it (a fallen moon or a passing comet), NULL for a resonance ring
    atmosphere_colour           INT UNSIGNED,  -- v161 (DB.96): the colour of the body's air, sRGB as red*65536 + green*256 + blue, NULL with no air
    lid_regime                  VARCHAR(10) NULL,  -- v169 (GEN.721): mobile (plate tectonics) or stagnant (one rigid plate), drawn once, NULL for a gas giant
    day_share_equator           DOUBLE NULL,  -- v170 (GEN.767): share of the time the star is up at the equator over a year, 0 to 1, NULL for a world that turns once an orbit
    day_share_mean              DOUBLE NULL,  -- v170 (GEN.767): share of the time the star is up, averaged over the surface, 0 to 1
    parent_shadow_share         DOUBLE NULL,  -- v170 (GEN.767): a moon only, the share of its time in its planet's shadow, 0 to 1
    solar_day_hours             DOUBLE NULL,  -- v170 (GEN.767): hours from one noon to the next, NULL for a world that turns once an orbit
    daylight_hours              DOUBLE NULL,  -- v170 (GEN.767): hours of daylight in one solar day at the equator, NULL for a world that turns once an orbit
    polar_circle_deg            DOUBLE NULL,  -- v170 (GEN.767): latitude beyond which the star stays up or down for at least a day, degrees
    life_class                  TINYINT UNSIGNED NULL,  -- v175 (GEN.889): the number of the life class that lives on the body (physics/life_classes.py), NULL for none
    life_class_variant          CHAR(1) NULL,  -- v175: its pigment or habitat letter (1a, 1b, 1c, 3a ...), NULL for a class with none
    life_class_best             VARCHAR(2) NULL,  -- v175: the id of the best-supported class, NULL when no class is supported
    life_class_best_score       TINYINT UNSIGNED NULL,  -- v175: its score, 0 to 9
    life_animal_grade           TINYINT UNSIGNED NULL,  -- v175 (GEN.890): 0 no animals, 1 sessile, 2 small and mobile, 3 large and active
    life_complexity             DOUBLE NULL,  -- v175: the complexity index of its life class, 1 for an Earth-like biosphere
    tide_range_spring_m         DOUBLE NULL,  -- v191 (GEN.1034): open-ocean equilibrium spring tide range, metres, 0 under an ice cover, NULL without a sea
    tide_range_neap_m           DOUBLE NULL,  -- v191 (GEN.1034): the neap range, metres
    ocean_liquid                VARCHAR(40) NULL,  -- v191 (GEN.1034): the liquid of the sea the tide is of (water until GEN.1035), NULL without a sea

    UNIQUE KEY uq_moons_uid (uid),
    CONSTRAINT fk_moons_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    CONSTRAINT fk_moons_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_moons_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,
    CONSTRAINT fk_moons_rogue_planet
        FOREIGN KEY (rogue_planet_id) REFERENCES rogue_planets(id) ON DELETE CASCADE,
    KEY idx_moons_rogue_planet_id (rogue_planet_id),
    KEY idx_moons_planet_id (planet_id),
    KEY idx_moons_next_update_due (next_update_due),
    KEY idx_moons_star_system_id (star_system_id),
    KEY idx_moons_star_id (star_id),
    KEY idx_moons_name (name),
    FULLTEXT KEY ft_moons_name (name),
    KEY idx_moons_planet_class (planet_class),
    KEY idx_moons_body_type (body_type),
    KEY idx_moons_life_class (life_class),
    KEY idx_moons_radius_km (radius_km),
    CONSTRAINT chk_moons_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),
    KEY idx_moons_parent (parent_type, parent_id),
    CONSTRAINT chk_moons_dynamical_state CHECK (dynamical_state IN ('bound', 'unbound'))
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

-- A giant's layers, outside in (DB.89, GEN.642): one row per layer of the interior drawn for the moon, in Earth masses per material.
CREATE TABLE IF NOT EXISTS moon_interior_layers (
    id                      BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    moon_id                 BIGINT UNSIGNED NOT NULL,
    layer_index             SMALLINT UNSIGNED NOT NULL,
    kind                    VARCHAR(24) NOT NULL CHECK (kind IN ('molecular_envelope', 'metallic_envelope', 'helium_gradient', 'helium_shell', 'dilute_core', 'compact_core', 'ice_mantle', 'rock_core')),
    state                   VARCHAR(16) NOT NULL CHECK (state IN ('molecular_fluid', 'metallic_fluid', 'ionic_fluid', 'superionic', 'molten')),
    top_pressure_gpa        DOUBLE NOT NULL,
    bottom_pressure_gpa     DOUBLE NOT NULL,
    mass_earth              DOUBLE NOT NULL,
    top_radius_fraction     DOUBLE NOT NULL,
    bottom_radius_fraction  DOUBLE NOT NULL,
    hydrogen                DOUBLE NOT NULL,
    helium                  DOUBLE NOT NULL,
    iron                    DOUBLE NOT NULL,
    rock                    DOUBLE NOT NULL,
    water_liquid            DOUBLE NOT NULL,
    ice_x                   DOUBLE NOT NULL,
    nh3_ice                 DOUBLE NOT NULL,
    ch4_ice                 DOUBLE NOT NULL,

    UNIQUE KEY uq_moon_interior_layers_index (moon_id, layer_index),
    CONSTRAINT fk_moon_interior_layers_moon
        FOREIGN KEY (moon_id) REFERENCES moons(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- The air of a moon (v155, DB.71, GEN.424): one row per gas above zero, its formula in the gas table (`physics/data/gases.json`) and its partial pressure.
-- A class's balance gas (helium, sodium, ammonia) is a row like the others. No row means no air (or a row made before GEN.85).
CREATE TABLE IF NOT EXISTS moon_gases (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    moon_id  BIGINT UNSIGNED NOT NULL,
    gas         VARCHAR(16) NOT NULL,
    partial_pa  DOUBLE NOT NULL CHECK (partial_pa > 0),

    CONSTRAINT fk_moon_gases_moon
        FOREIGN KEY (moon_id) REFERENCES moons(id) ON DELETE CASCADE,
    UNIQUE KEY uq_moon_gases (moon_id, gas)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- The life classes a moon supports (v175, GEN.889): one row for each class of `physics/life_classes.py` whose score for the body is 1 or more (0 to 9). Indexed on the class and the score so the
-- map filter and the search facet are cheap.
CREATE TABLE IF NOT EXISTS moon_life_class_scores (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    moon_id   BIGINT UNSIGNED NOT NULL,
    life_class  VARCHAR(2) NOT NULL,
    score       TINYINT UNSIGNED NOT NULL CHECK (score BETWEEN 1 AND 9),

    CONSTRAINT fk_moon_life_class_scores_moon
        FOREIGN KEY (moon_id) REFERENCES moons(id) ON DELETE CASCADE,
    UNIQUE KEY uq_moon_life_class_scores (moon_id, life_class),
    KEY idx_moon_life_class_scores_class (life_class, score)
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
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),
    position_x_km               DOUBLE NULL,  -- v102 (DB.55): system-frame offset from the system centre, km
    position_y_km               DOUBLE NULL,  -- v102 (DB.55): system-frame offset from the system centre, km
    position_z_km               DOUBLE NULL,  -- v102 (DB.55): system-frame offset from the system centre, km
    velocity_x_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_y_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_z_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    mass_kg                     DOUBLE NULL,  -- v102 (DB.55): point mass, a stated zero where the mass is not known
    epoch_unix                  DOUBLE NULL,  -- v102 (DB.55): when the position and velocity hold
    next_update_due             DOUBLE NULL,  -- v102 (DB.55): when the object next moves

    parent_type                 VARCHAR(16) NOT NULL DEFAULT 'system',  -- v140 (DB.73): what it orbits, see the v140 note
    parent_id                   BIGINT UNSIGNED NULL,  -- v140 (DB.73)
    chondrite_group             VARCHAR(4) NULL,  -- v190 (GEN.475): the nearest chondrite group of the belt's mean composition (CI, CM, H, L ...), `none` far from every group, NULL before v190

    UNIQUE KEY uq_asteroid_belts_uid (uid),
    KEY idx_asteroid_belts_parent (parent_type, parent_id),
    CONSTRAINT chk_asteroid_belts_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),
    CONSTRAINT fk_asteroid_belts_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_asteroid_belts_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,
    KEY idx_asteroid_belts_star_system_id (star_system_id),
    KEY idx_asteroid_belts_star_id (star_id),
    KEY idx_asteroid_belts_density (density)
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
-- comets -- v19, a comet gravitationally bound to a star
-- (planetgen/generation/comet.py:Comet), propagated via real two-body
-- Kepler/Barker orbital mechanics (planetgen/physics/kepler.py) --
-- contrast interstellar_comets below, an unbound object on a fixed
-- hyperbolic trajectory. Deliberately NOT part of the planets table (no
-- orbital_index/shared orbital-slot ordering): a comet's distance_km is
-- its current, continuously-varying position along its orbit, not a
-- fixed slot the way a planet's/belt's own distance is (see
-- systemData.StarSystem._generate_comets' docstring). star_system_id/
-- star_id follow planets.star_id's own real semantics (see that column's
-- comment above): a real stars.id for a single star or a 'wide' binary's
-- comet (each orbits one specific constituent star); NULL only for a
-- 'close' (P-type) binary's comet, which orbits the merged pair, not any
-- one individually-stored star row.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS comets (
    id                          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_system_id              BIGINT UNSIGNED NOT NULL,
    star_id                     BIGINT UNSIGNED,
    name                        VARCHAR(255) NOT NULL,
    orbit_type                  VARCHAR(16) NOT NULL CHECK (orbit_type IN ('elliptical', 'parabolic')),
    -- period_class: only set for orbit_type = 'elliptical' (flavor/
    -- plausibility metadata only -- see program_constants.COMET_PERIOD_CLASSES).
    period_class                VARCHAR(16) CHECK (period_class IN ('jupiter_family', 'halley_type', 'long_period')),
    nucleus_diameter_km         DOUBLE NOT NULL,
    -- Searchable composition summary, same role as interstellar_comets.composition_summary
    -- below -- structured per-component detail lives in comet_composition.
    composition_summary         TEXT NOT NULL,
    perihelion_distance_km      DOUBLE NOT NULL,
    eccentricity                DOUBLE NOT NULL,
    inclination_deg             DOUBLE NOT NULL,
    arg_periapsis_deg           DOUBLE NOT NULL,
    ascending_node_deg          DOUBLE NOT NULL,
    -- orbital_period_years/mean_anomaly_deg: only set for orbit_type = 'elliptical'.
    orbital_period_years        DOUBLE,
    mean_anomaly_deg            DOUBLE,
    -- parabolic_mean_anomaly/min_update_interval_years: parabolic_mean_anomaly
    -- only set for orbit_type = 'parabolic'; min_update_interval_years only
    -- for 'elliptical' (a parabolic comet's anomaly doesn't wrap, so it has
    -- no periodic floating-point-resolution floor to guard -- see Comet's
    -- own min_update_interval_years docstring).
    parabolic_mean_anomaly      DOUBLE,
    min_update_interval_years   DOUBLE,
    primary_mass_solar          DOUBLE NOT NULL,
    is_active                   TINYINT(1) NOT NULL CHECK (is_active IN (0, 1)),
    -- Current derived orbital state (Comet.update_orbital_state) --
    -- recomputed by _db.advance_comet_orbits as mean_anomaly_deg/
    -- parabolic_mean_anomaly advance over time, the Kepler/Barker-equation
    -- analog of planets.position_x/y/z_km's own periodic recomputation.
    distance_km                 DOUBLE NOT NULL,
    position_x_km               DOUBLE NOT NULL,
    position_y_km               DOUBLE NOT NULL,
    position_z_km               DOUBLE NOT NULL,
    orbital_speed_kms           DOUBLE NOT NULL,
    -- v60: velocity relative to the orbital anchor, km/s, in step with the position (see the header comment's "v60" note).
    velocity_x_kms              DOUBLE NOT NULL DEFAULT 0,
    velocity_y_kms              DOUBLE NOT NULL DEFAULT 0,
    velocity_z_kms              DOUBLE NOT NULL DEFAULT 0,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,
    -- v68 (GEN.104): the spin, see header comment.
    rotation_period_hours       DOUBLE,
    spin_axis_x                 DOUBLE,
    spin_axis_y                 DOUBLE,
    spin_axis_z                 DOUBLE,
    axial_tilt_deg              DOUBLE,
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),
    mass_kg                     DOUBLE NULL,  -- v102 (DB.55): point mass, a stated zero where the mass is not known

    parent_type                 VARCHAR(16) NOT NULL DEFAULT 'system',  -- v140 (DB.73): what it orbits, see the v140 note
    parent_id                   BIGINT UNSIGNED NULL,  -- v140 (DB.73)

    UNIQUE KEY uq_comets_uid (uid),
    KEY idx_comets_parent (parent_type, parent_id),
    CONSTRAINT chk_comets_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),
    CONSTRAINT fk_comets_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_comets_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE SET NULL,
    KEY idx_comets_star_system_id (star_system_id),
    KEY idx_comets_next_update_due (next_update_due),
    KEY idx_comets_star_id (star_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Structured per-component detail behind composition_summary above --
-- a plain component list (no concentration gradient), mirroring
-- interstellar_comet_composition's identical shape (Comet.composition and
-- InterstellarComet.composition are both plain lists, not the
-- (component, concentration) pairs asteroid_belt_composition uses).
CREATE TABLE IF NOT EXISTS comet_composition (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    comet_id    BIGINT UNSIGNED NOT NULL,
    position    INT NOT NULL,
    component   VARCHAR(64) NOT NULL,

    CONSTRAINT fk_comet_composition_comet
        FOREIGN KEY (comet_id) REFERENCES comets(id) ON DELETE CASCADE,
    KEY idx_comet_composition_comet_id (comet_id)
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
    -- v21: the sector a standalone black hole was generated as part of
    -- (sectorGen.generate_sector_phenomena), plus its own galaxy-frame
    -- center -- see this file's "v21" header note. Always NULL when
    -- star_id is set (an anchored remnant belongs to its own StarSystem,
    -- not directly to a sector). NULL together: never placed in the
    -- galaxy (a phenomenonGen.py --type black-hole run with no
    -- --sector-id, exactly like before v21).
    sector_id                 BIGINT UNSIGNED,
    name                      VARCHAR(255) NOT NULL,
    -- v36: 'stellar', 'intermediate' (IMBH, ~1e2-1e5 Msun) or
    -- 'supermassive' (a galaxy's quiescent central black hole, placed by
    -- generate.add_galactic_nucleus; no galactic orbit).
    mass_class                VARCHAR(16) NOT NULL DEFAULT 'stellar',
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
    galactic_orbital_period_gy           DOUBLE,
    galactic_orbital_phase_deg           DOUBLE,
    galactic_min_update_interval_years   DOUBLE,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,
    -- v68 (GEN.104): the spin, see header comment.
    spin_axis_x                 DOUBLE,
    spin_axis_y                 DOUBLE,
    spin_axis_z                 DOUBLE,
    axial_tilt_deg              DOUBLE,

    -- v21: this black hole's own galaxy-frame center -- see nebulae's
    -- identical "v18" column comment above (same shape, added here a
    -- version later).
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    velocity_x_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_y_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_z_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    KEY idx_black_holes_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    galactic_radius_pc  DOUBLE,
    -- v41: the sector octant the center sits in (star_systems.quadrant's labels).
    quadrant            VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- v27: row timestamps -- see the header comment's "v27" note.
    -- v39: the innermost nebula or supernova remnant this sits inside
    -- (at most one set; _db.refresh_containment).
    inside_nebula_id    BIGINT UNSIGNED,
    inside_remnant_id   BIGINT UNSIGNED,

    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),

    parent_type                 VARCHAR(16) NOT NULL DEFAULT 'free',  -- v140 (DB.73): what it orbits, see the v140 note
    parent_id                   BIGINT UNSIGNED NULL,  -- v140 (DB.73)

    UNIQUE KEY uq_black_holes_uid (uid),
    KEY idx_black_holes_parent (parent_type, parent_id),
    CONSTRAINT chk_black_holes_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),
    KEY idx_black_holes_name (name),  -- v40: name uniqueness lookups
    KEY idx_black_holes_placed_name (star_id, ((center_x_pc IS NOT NULL)), name),   -- v186 (PERF.211): the Phenomena table's sorts read in index order
    KEY idx_black_holes_descriptor_name (star_id, ((CASE WHEN has_accretion_disk THEN 'accreting' ELSE 'quiescent' END) COLLATE utf8mb4_unicode_ci), name),
    KEY idx_black_holes_next_update_due (next_update_due),
    CONSTRAINT chk_black_holes_mass_class CHECK (mass_class IN ('stellar', 'intermediate', 'supermassive')),
    CONSTRAINT chk_black_holes_placement CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    CONSTRAINT fk_black_holes_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE CASCADE,
    CONSTRAINT fk_black_holes_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE SET NULL,
    KEY idx_black_holes_star_id (star_id),
    KEY idx_black_holes_sector_id (sector_id),
    KEY idx_black_holes_galactic_radius_pc (galactic_radius_pc),
    -- v26: lets `queryDb.phenomena_near_sector`'s bounding-box query
    -- range-scan on center_x_pc instead of a full table scan -- see the
    -- header comment's "v26" note.
    KEY idx_black_holes_center (center_x_pc, center_y_pc, center_z_pc),
    -- v27: see the header comment's "v27" note.
    KEY idx_black_holes_modified_at (modified_at),
    -- v39: containment -- see the header comment's "v39" note.
    CONSTRAINT fk_black_holes_inside_nebula
        FOREIGN KEY (inside_nebula_id) REFERENCES nebulae(id) ON DELETE SET NULL,
    CONSTRAINT fk_black_holes_inside_remnant
        FOREIGN KEY (inside_remnant_id) REFERENCES supernova_remnants(id) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS neutron_stars (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    star_id                   BIGINT UNSIGNED,
    -- v21: same sector-link/placement convention as black_holes above.
    sector_id                 BIGINT UNSIGNED,
    name                      VARCHAR(255) NOT NULL,
    mass_solar                DOUBLE NOT NULL,
    radius_km                 DOUBLE NOT NULL,
    spin_period_ms            DOUBLE NOT NULL,
    magnetic_field_gauss      DOUBLE NOT NULL,
    -- v96 (GEN.153): 'magnetar' joins the pulsar types.
    pulsar_type               VARCHAR(16) NOT NULL CHECK (pulsar_type IN ('young', 'millisecond', 'magnetar', 'non-pulsing')),
    surface_temperature_k     DOUBLE NOT NULL,
    luminosity_w              DOUBLE NOT NULL,
    age_gy                    DOUBLE NOT NULL,
    -- v17: same nullable-when-anchored convention as black_holes above.
    galactic_orbital_period_gy           DOUBLE,
    galactic_orbital_phase_deg           DOUBLE,
    galactic_min_update_interval_years   DOUBLE,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,
    -- v68 (GEN.104): the spin, see header comment.
    spin_axis_x                 DOUBLE,
    spin_axis_y                 DOUBLE,
    spin_axis_z                 DOUBLE,
    axial_tilt_deg              DOUBLE,

    -- v21: same galaxy-frame center convention as black_holes above.
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    velocity_x_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_y_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_z_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    KEY idx_neutron_stars_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    galactic_radius_pc  DOUBLE,
    -- v41: the sector octant the center sits in (star_systems.quadrant's labels).
    quadrant            VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- v27: row timestamps -- see the header comment's "v27" note.
    -- v39: the innermost nebula or supernova remnant this sits inside
    -- (at most one set; _db.refresh_containment).
    inside_nebula_id    BIGINT UNSIGNED,
    inside_remnant_id   BIGINT UNSIGNED,

    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),

    parent_type                 VARCHAR(16) NOT NULL DEFAULT 'free',  -- v140 (DB.73): what it orbits, see the v140 note
    parent_id                   BIGINT UNSIGNED NULL,  -- v140 (DB.73)

    UNIQUE KEY uq_neutron_stars_uid (uid),
    KEY idx_neutron_stars_parent (parent_type, parent_id),
    CONSTRAINT chk_neutron_stars_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),
    KEY idx_neutron_stars_name (name),  -- v40: name uniqueness lookups
    KEY idx_neutron_stars_placed_name (star_id, ((center_x_pc IS NOT NULL)), name),   -- v186 (PERF.211): the Phenomena table's sorts read in index order
    KEY idx_neutron_stars_descriptor_name (star_id, pulsar_type, name),
    KEY idx_neutron_stars_next_update_due (next_update_due),
    CONSTRAINT chk_neutron_stars_placement CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    CONSTRAINT fk_neutron_stars_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE CASCADE,
    CONSTRAINT fk_neutron_stars_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE SET NULL,
    KEY idx_neutron_stars_star_id (star_id),
    KEY idx_neutron_stars_sector_id (sector_id),
    KEY idx_neutron_stars_galactic_radius_pc (galactic_radius_pc),
    -- v26: see black_holes' identical "v26" index comment above.
    KEY idx_neutron_stars_center (center_x_pc, center_y_pc, center_z_pc),
    -- v27: see the header comment's "v27" note.
    KEY idx_neutron_stars_modified_at (modified_at),
    -- v39: containment -- see the header comment's "v39" note.
    CONSTRAINT fk_neutron_stars_inside_nebula
        FOREIGN KEY (inside_nebula_id) REFERENCES nebulae(id) ON DELETE SET NULL,
    CONSTRAINT fk_neutron_stars_inside_remnant
        FOREIGN KEY (inside_remnant_id) REFERENCES supernova_remnants(id) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- nebulae -- v16 exotic phenomenon, always standalone. `sector_id` links
-- it to the nearest sector, populated by phenomenonGen.py --sector-id or
-- (v21 on) sectorGen.py's own per-sector generation -- see this file's
-- "v18"/"v21" header notes.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS nebulae (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    -- v18: the nearest already-generated sector to center_x/y/z_pc below --
    -- a convenience "home" link, not this phenomenon's real geometry (a
    -- nebula's sphere may overlap several sectors' cubes, or none). NULL
    -- iff center_x/y/z_pc are NULL (never placed in the galaxy at all).
    sector_id         BIGINT UNSIGNED,
    name              VARCHAR(255) NOT NULL,
    -- v38: letter class A-Q (program_constants.NEBULA_CLASSES); nebula_type
    -- is the class's family.
    nebula_class      CHAR(1) NOT NULL,
    nebula_type       VARCHAR(16) NOT NULL,
    radius_ly         DOUBLE NOT NULL,
    composition       TEXT NOT NULL,
    formation_cause   TEXT NOT NULL,
    -- v38: the class's contents -- see the header comment's "v38" note.
    dominant_species  VARCHAR(255) NOT NULL,
    density_cm3       DOUBLE NOT NULL,
    temperature_k     DOUBLE NOT NULL,
    extinction_av     DOUBLE NOT NULL,
    -- v87 (DB.26, GEN.150): age, ionizing photon rate and host stars, see header comment.
    age_myr           DOUBLE,
    q_total           DOUBLE,
    hosts             VARCHAR(512),
    -- v17: always populated (a nebula is always standalone). See this
    -- file's "v17" header note.
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,

    -- v18: this nebula's own galaxy-frame center -- see this file's "v18"
    -- header note. NULL together: never placed in the galaxy (the default,
    -- unless phenomenonGen.py --sector-id was given).
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    velocity_x_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_y_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_z_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    mass_kg                     DOUBLE NULL,  -- v102 (DB.55): point mass, a stated zero where the mass is not known
    KEY idx_nebulae_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    galactic_radius_pc  DOUBLE,
    -- v41: the sector octant the center sits in (star_systems.quadrant's labels).
    quadrant            VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- v56: the shape (see the header's "v56" note); NULL for a nebula
    -- saved before it.
    shape_axis_x          DOUBLE,
    shape_axis_y          DOUBLE,
    shape_axis_z          DOUBLE,
    shape_angle_x         DOUBLE,
    shape_angle_y         DOUBLE,
    shape_angle_z         DOUBLE,
    shape_warp_amplitude  DOUBLE,
    shape_warp_frequency  DOUBLE,
    shape_warp_octaves    INT,
    shape_noise_seed      BIGINT UNSIGNED,
    shape_iso             DOUBLE,
    shape_scale           DOUBLE,

    -- v27: row timestamps -- see the header comment's "v27" note.
    -- v39: the innermost nebula or supernova remnant this sits inside
    -- (at most one set; _db.refresh_containment).
    inside_nebula_id    BIGINT UNSIGNED,
    inside_remnant_id   BIGINT UNSIGNED,

    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),

    UNIQUE KEY uq_nebulae_uid (uid),
    KEY idx_nebulae_name (name),  -- v40: name uniqueness lookups
    KEY idx_nebulae_placed_name (((center_x_pc IS NOT NULL)), name),   -- v186 (PERF.211): the Phenomena table's sorts read in index order
    KEY idx_nebulae_descriptor_name (nebula_type, name),
    KEY idx_nebulae_radius_name (radius_ly, name),
    KEY idx_nebulae_next_update_due (next_update_due),

    -- Named explicitly (unlike sectors' own identical v4 CHECK above) so
    -- `_migrate_v17_to_v18` can add the exact same constraint by name to a
    -- migrated (not freshly created) database -- an anonymous CHECK gets
    -- an opaque, position-dependent auto-generated name (`nebulae_chk_2`,
    -- shifting if another CHECK is ever added/removed above it), which a
    -- migration step has no reliable way to reproduce or later reference.
    CONSTRAINT chk_nebulae_type CHECK (nebula_type IN ('diffuse', 'emission', 'reflection', 'planetary', 'dark')),
    CONSTRAINT chk_nebulae_placement CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    CONSTRAINT fk_nebulae_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE SET NULL,
    KEY idx_nebulae_sector_id (sector_id),
    KEY idx_nebulae_galactic_radius_pc (galactic_radius_pc),
    -- v26: see black_holes' identical "v26" index comment above.
    KEY idx_nebulae_center (center_x_pc, center_y_pc, center_z_pc),
    -- v27: see the header comment's "v27" note.
    KEY idx_nebulae_modified_at (modified_at),
    -- v39: containment -- see the header comment's "v39" note.
    CONSTRAINT fk_nebulae_inside_nebula
        FOREIGN KEY (inside_nebula_id) REFERENCES nebulae(id) ON DELETE SET NULL,
    CONSTRAINT fk_nebulae_inside_remnant
        FOREIGN KEY (inside_remnant_id) REFERENCES supernova_remnants(id) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- v56: a nebula's metaballs (GEN.75): 4 to 8 per nebula, in nebula-shape
-- units (`galaxy/nebula_shape.py`). `ball_index` is their order.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS nebula_shape_balls (
    nebula_id   BIGINT UNSIGNED NOT NULL,
    ball_index  INT NOT NULL,
    center_x    DOUBLE NOT NULL,
    center_y    DOUBLE NOT NULL,
    center_z    DOUBLE NOT NULL,
    radius      DOUBLE NOT NULL,

    PRIMARY KEY (nebula_id, ball_index),
    CONSTRAINT fk_nebula_shape_balls_nebula
        FOREIGN KEY (nebula_id) REFERENCES nebulae(id) ON DELETE CASCADE
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
    -- v38: letter class R-W (program_constants.NEBULA_CLASSES); the class
    -- fixes the morphology.
    remnant_class                     CHAR(1) NOT NULL,
    morphology                        VARCHAR(16) NOT NULL CHECK (morphology IN ('shell', 'plerion', 'composite')),
    age_years                         DOUBLE NOT NULL,
    radius_ly                         DOUBLE NOT NULL,
    progenitor_type                   VARCHAR(16) NOT NULL CHECK (progenitor_type IN ('Type Ia', 'core-collapse')),
    compact_remnant_kind              VARCHAR(16) CHECK (compact_remnant_kind IN ('black_hole', 'neutron_star')),
    compact_remnant_black_hole_id     BIGINT UNSIGNED,
    compact_remnant_neutron_star_id   BIGINT UNSIGNED,
    -- v38: the class's contents, as on nebulae.
    dominant_species                  VARCHAR(255) NOT NULL,
    density_cm3                       DOUBLE NOT NULL,
    temperature_k                     DOUBLE NOT NULL,
    extinction_av                     DOUBLE NOT NULL,
    -- v17: always populated (a supernova remnant is always standalone).
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,

    -- v28: this remnant's own galaxy-frame center -- see this file's "v28"
    -- header note. NULL together: never placed in the galaxy.
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    velocity_x_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_y_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_z_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    mass_kg                     DOUBLE NULL,  -- v102 (DB.55): point mass, a stated zero where the mass is not known
    KEY idx_supernova_remnants_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    galactic_radius_pc  DOUBLE,
    -- v41: the sector octant the center sits in (star_systems.quadrant's labels).
    quadrant            VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- v27: row timestamps -- see the header comment's "v27" note.
    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),

    UNIQUE KEY uq_supernova_remnants_uid (uid),
    KEY idx_supernova_remnants_name (name),  -- v40: name uniqueness lookups
    KEY idx_supernova_remnants_placed_name (((center_x_pc IS NOT NULL)), name),   -- v186 (PERF.211): the Phenomena table's sorts read in index order
    KEY idx_supernova_remnants_descriptor_name (morphology, name),
    KEY idx_supernova_remnants_radius_name (radius_ly, name),
    KEY idx_supernova_remnants_next_update_due (next_update_due),

    -- v28: named explicitly -- see `nebulae`'s identical "v18" CHECK comment.
    CONSTRAINT chk_supernova_remnants_placement CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    CONSTRAINT fk_supernova_remnants_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    CONSTRAINT fk_supernova_remnants_black_hole
        FOREIGN KEY (compact_remnant_black_hole_id) REFERENCES black_holes(id) ON DELETE SET NULL,
    CONSTRAINT fk_supernova_remnants_neutron_star
        FOREIGN KEY (compact_remnant_neutron_star_id) REFERENCES neutron_stars(id) ON DELETE SET NULL,
    KEY idx_supernova_remnants_sector_id (sector_id),
    -- v28: see the header comment's "v28" note.
    KEY idx_supernova_remnants_galactic_radius_pc (galactic_radius_pc),
    KEY idx_supernova_remnants_center (center_x_pc, center_y_pc, center_z_pc),
    -- v27: see the header comment's "v27" note.
    KEY idx_supernova_remnants_modified_at (modified_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- quasars -- v31 exotic phenomenon: a galaxy's active nucleus, always at
-- the galactic center. See this file's "v31" header note.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS quasars (
    id                              BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id                       BIGINT UNSIGNED,
    name                            VARCHAR(255) NOT NULL,
    black_hole_mass_solar           DOUBLE NOT NULL,
    event_horizon_radius_km         DOUBLE NOT NULL,
    eddington_ratio                 DOUBLE NOT NULL,
    luminosity_w                    DOUBLE NOT NULL,
    accretion_rate_solar_per_year   DOUBLE NOT NULL,
    broad_line_region_light_days    DOUBLE NOT NULL,
    is_radio_loud                   TINYINT(1) NOT NULL CHECK (is_radio_loud IN (0, 1)),
    -- NULL when radio-quiet (no jets).
    jet_length_ly                   DOUBLE,
    active_age_years                DOUBLE NOT NULL,

    -- The galactic center, when placed. NULL together: never placed.
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    velocity_x_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_y_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_z_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    epoch_unix                  DOUBLE NULL,  -- v102 (DB.55): when the position and velocity hold
    next_update_due             DOUBLE NULL,  -- v102 (DB.55): when the object next moves
    KEY idx_quasars_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    galactic_radius_pc  DOUBLE,
    -- v41: the sector octant the center sits in (star_systems.quadrant's labels).
    quadrant            VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- v27-style row timestamps -- see the header comment's "v27" note.
    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),

    UNIQUE KEY uq_quasars_uid (uid),
    KEY idx_quasars_name (name),  -- v40: name uniqueness lookups
    KEY idx_quasars_placed_name (((center_x_pc IS NOT NULL)), name),   -- v186 (PERF.211): the Phenomena table's sorts read in index order
    KEY idx_quasars_descriptor_name (((CASE WHEN is_radio_loud THEN 'radio-loud' ELSE 'radio-quiet' END) COLLATE utf8mb4_unicode_ci), name),
    CONSTRAINT chk_quasars_placement CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    CONSTRAINT fk_quasars_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    KEY idx_quasars_sector_id (sector_id),
    KEY idx_quasars_galactic_radius_pc (galactic_radius_pc),
    KEY idx_quasars_center (center_x_pc, center_y_pc, center_z_pc),
    KEY idx_quasars_modified_at (modified_at)
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
    -- v47: its PLANET_CLASSES letter (GEN.8); NULL for a brown dwarf.
    planet_class        VARCHAR(4),
    -- v37: the mass bin it was drawn from (program_constants.
    -- ROGUE_PLANET_MASS_BINS), or 'brown-dwarf' (13-80 Mjup).
    mass_bin            VARCHAR(16) NOT NULL DEFAULT 'terrestrial',
    mass_kg             DOUBLE NOT NULL,
    radius_km           DOUBLE NOT NULL,
    composition         TEXT NOT NULL,
    has_internal_heat   TINYINT(1) NOT NULL CHECK (has_internal_heat IN (0, 1)),
    has_moons           TINYINT(1) NOT NULL CHECK (has_moons IN (0, 1)),
    -- v48: surface conditions (rogueSurface.rogue_surface_conditions).
    age_gy                   DOUBLE,
    internal_heat_flux_w_m2  DOUBLE,
    effective_temperature_k  DOUBLE,
    surface_regime           VARCHAR(24),
    surface_temperature_k    DOUBLE,
    surface_pressure_pa      DOUBLE,
    ice_shell_thickness_km   DOUBLE,
    ocean_depth_km           DOUBLE,
    has_liquid_water         TINYINT(1),
    hp_ice_km                DOUBLE,  -- v74 (GEN.88)
    -- v171 (GEN.836): the ring system, as on `planets` (v160); NULL for a planet with none.
    ring_kind                VARCHAR(8) CHECK (ring_kind IN ('dusty', 'narrow', 'broad')),
    ring_origin              VARCHAR(16) CHECK (ring_origin IN ('moon_erosion', 'disrupted_moon', 'captured_body')),
    ring_inner_km            DOUBLE,
    ring_outer_km            DOUBLE,
    ring_mass_kg             DOUBLE,
    ring_optical_depth       DOUBLE,
    ring_ice_fraction        DOUBLE,
    ring_age_yr              DOUBLE,  -- v192 (GEN.1057): as on `planets`
    -- v17: always populated (a rogue planet is always standalone).
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,
    -- v68 (GEN.104): the spin, see header comment.
    rotation_period_hours       DOUBLE,
    spin_axis_x                 DOUBLE,
    spin_axis_y                 DOUBLE,
    spin_axis_z                 DOUBLE,
    axial_tilt_deg              DOUBLE,

    -- v28: this rogue planet's own galaxy-frame center -- see this file's "v28"
    -- header note. NULL together: never placed in the galaxy.
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    velocity_x_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_y_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_z_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    KEY idx_rogue_planets_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    galactic_radius_pc  DOUBLE,
    -- v41: the sector octant the center sits in (star_systems.quadrant's labels).
    quadrant            VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- v27: row timestamps -- see the header comment's "v27" note.
    -- v39: the innermost nebula or supernova remnant this sits inside
    -- (at most one set; _db.refresh_containment).
    inside_nebula_id    BIGINT UNSIGNED,
    inside_remnant_id   BIGINT UNSIGNED,

    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),

    parent_type                 VARCHAR(16) NOT NULL DEFAULT 'free',  -- v140 (DB.73): what it orbits, see the v140 note
    parent_id                   BIGINT UNSIGNED NULL,  -- v140 (DB.73)
    phi4                        DOUBLE NULL,  -- v158 (GEN.701): the habitability score of the surface a visitor would find (rogue_habitability.py)
    phi4_pressure               DOUBLE NULL,  -- v153
    phi4_temperature            DOUBLE NULL,  -- v153
    phi4_air                    DOUBLE NULL,  -- v153
    phi4_water                  DOUBLE NULL,  -- v153
    phi4_radiation              DOUBLE NULL,  -- v153
    phi4_suitability            DOUBLE NULL,  -- v153
    tier_pressure               TINYINT NULL,  -- v153
    tier_temperature            TINYINT NULL,  -- v153
    tier_air                    TINYINT NULL,  -- v153
    tier_water                  TINYINT NULL,  -- v153
    tier_radiation              TINYINT NULL,  -- v153
    tier_suitability            TINYINT NULL,  -- v153
    life_microbial              TINYINT UNSIGNED NULL,  -- v153
    life_flora                  TINYINT UNSIGNED NULL,  -- v153
    life_fauna                  TINYINT UNSIGNED NULL,  -- v153
    life_sentient               TINYINT UNSIGNED NULL,  -- v153
    l_solv                      DOUBLE NULL,  -- v153
    l_chem                      DOUBLE NULL,  -- v153
    l_ener                      DOUBLE NULL,  -- v153
    l_rad                       DOUBLE NULL,  -- v153
    eq_body                     TINYINT NULL,  -- v153
    eq_body_kind                TINYINT NULL,  -- v153
    eq_body_rec                 TINYINT NULL,  -- v153
    eq_body_rec_kind            TINYINT NULL,  -- v153
    eq_rad_suit                 TINYINT NULL,  -- v153
    eq_rad_suit_rec             TINYINT NULL,  -- v153
    eq_mobility                 TINYINT NULL,  -- v153
    eq_dust                     TINYINT NULL,  -- v153
    eq_gas_filter               TINYINT NULL,  -- v153
    eq_co2_scrub                TINYINT NULL,  -- v153
    eq_supply                   TINYINT NULL,  -- v153
    eq_supply_rec               TINYINT NULL,  -- v153
    eq_shelter                  TINYINT NULL,  -- v153
    eq_shelter_rec              TINYINT NULL,  -- v153
    shelter_column_g_cm2        DOUBLE NULL,  -- v153
    energy_flux_w_m2            DOUBLE NULL,  -- v153

    UNIQUE KEY uq_rogue_planets_uid (uid),
    KEY idx_rogue_planets_parent (parent_type, parent_id),
    CONSTRAINT chk_rogue_planets_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),
    KEY idx_rogue_planets_name (name),  -- v40: name uniqueness lookups
    KEY idx_rogue_planets_placed_name (((center_x_pc IS NOT NULL)), name),   -- v186 (PERF.211): the Phenomena table's sorts read in index order
    KEY idx_rogue_planets_descriptor_name (((CASE WHEN mass_bin = 'brown-dwarf' THEN 'brown dwarf' WHEN planet_type = 'g' THEN 'gas giant' ELSE 'terrestrial' END) COLLATE utf8mb4_unicode_ci), name),
    KEY idx_rogue_planets_next_update_due (next_update_due),
    -- v28: named explicitly -- see `nebulae`'s identical "v18" CHECK comment.
    CONSTRAINT chk_rogue_planets_placement CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    CONSTRAINT fk_rogue_planets_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    KEY idx_rogue_planets_sector_id (sector_id),
    -- v28: see the header comment's "v28" note.
    KEY idx_rogue_planets_galactic_radius_pc (galactic_radius_pc),
    KEY idx_rogue_planets_center (center_x_pc, center_y_pc, center_z_pc),
    -- v27: see the header comment's "v27" note.
    KEY idx_rogue_planets_modified_at (modified_at),
    -- v39: containment -- see the header comment's "v39" note.
    CONSTRAINT fk_rogue_planets_inside_nebula
        FOREIGN KEY (inside_nebula_id) REFERENCES nebulae(id) ON DELETE SET NULL,
    CONSTRAINT fk_rogue_planets_inside_remnant
        FOREIGN KEY (inside_remnant_id) REFERENCES supernova_remnants(id) ON DELETE SET NULL
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
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,
    -- v68 (GEN.104): the spin, see header comment.
    rotation_period_hours       DOUBLE,
    spin_axis_x                 DOUBLE,
    spin_axis_y                 DOUBLE,
    spin_axis_z                 DOUBLE,
    axial_tilt_deg              DOUBLE,

    -- v28: this comet's own galaxy-frame center -- see this file's "v28"
    -- header note. NULL together: never placed in the galaxy.
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    velocity_x_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_y_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_z_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    mass_kg                     DOUBLE NULL,  -- v102 (DB.55): point mass, a stated zero where the mass is not known
    KEY idx_interstellar_comets_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    galactic_radius_pc  DOUBLE,
    -- v41: the sector octant the center sits in (star_systems.quadrant's labels).
    quadrant            VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- v27: row timestamps -- see the header comment's "v27" note.
    -- v39: the innermost nebula or supernova remnant this sits inside
    -- (at most one set; _db.refresh_containment).
    inside_nebula_id    BIGINT UNSIGNED,
    inside_remnant_id   BIGINT UNSIGNED,

    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),

    -- v28: named explicitly -- see `nebulae`'s identical "v18" CHECK comment.
    parent_type                 VARCHAR(16) NOT NULL DEFAULT 'free',  -- v140 (DB.73): what it orbits, see the v140 note
    parent_id                   BIGINT UNSIGNED NULL,  -- v140 (DB.73)

    UNIQUE KEY uq_interstellar_comets_uid (uid),
    KEY idx_interstellar_comets_name (name),  -- v178 (PERF.205): the Phenomena table pages by name in index order
    KEY idx_interstellar_comets_placed_name (((center_x_pc IS NOT NULL)), name),   -- v186 (PERF.211): the Phenomena table's sorts read in index order
    KEY idx_interstellar_comets_descriptor_name (((CASE WHEN is_active THEN 'active' ELSE 'dormant' END) COLLATE utf8mb4_unicode_ci), name),
    KEY idx_interstellar_comets_parent (parent_type, parent_id),
    CONSTRAINT chk_interstellar_comets_parent_type CHECK (parent_type IN ('star', 'planet', 'system', 'node', 'rogue_planet', 'free')),
    CONSTRAINT chk_interstellar_comets_placement CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    CONSTRAINT fk_interstellar_comets_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    KEY idx_interstellar_comets_sector_id (sector_id),
    KEY idx_interstellar_comets_next_update_due (next_update_due),
    -- v28: see the header comment's "v28" note.
    KEY idx_interstellar_comets_galactic_radius_pc (galactic_radius_pc),
    KEY idx_interstellar_comets_center (center_x_pc, center_y_pc, center_z_pc),
    -- v27: see the header comment's "v27" note.
    KEY idx_interstellar_comets_modified_at (modified_at),
    -- v39: containment -- see the header comment's "v39" note.
    CONSTRAINT fk_interstellar_comets_inside_nebula
        FOREIGN KEY (inside_nebula_id) REFERENCES nebulae(id) ON DELETE SET NULL,
    CONSTRAINT fk_interstellar_comets_inside_remnant
        FOREIGN KEY (inside_remnant_id) REFERENCES supernova_remnants(id) ON DELETE SET NULL
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
-- asteroid_belts, which always orbits a star). `sector_id` links it to
-- the nearest sector, the same convention `nebulae` uses (see this file's
-- "v18"/"v21" header notes).
-- `asteroid_field_composition` mirrors `asteroid_belt_composition`
-- exactly (same per-component/concentration shape -- both are generated
-- via the same shared `asteroidData.generate_asteroid_composition`).
-- ---------------------------------------------------------------------
-- v40: `name` holds the field's designation (`AF <field_class>-<sector>-<nn>`,
-- `asteroidFieldData.asteroid_field_designation`). Asteroid fields can sit
-- inside a nebula (inside_nebula_id, v39) but never contain anything.
CREATE TABLE IF NOT EXISTS asteroid_fields (
    id                    BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    -- v18: the nearest already-generated sector to center_x/y/z_pc below --
    -- see nebulae's identical "v18" column comment above.
    sector_id             BIGINT UNSIGNED,
    name                  VARCHAR(255) NOT NULL,
    -- v38: letter (composition and density) plus size digit, e.g. 'C3'
    -- (program_constants.ASTEROID_FIELD_COMPOSITIONS).
    field_class           VARCHAR(4) NOT NULL,
    composition_family    VARCHAR(16) NOT NULL,
    density               VARCHAR(16) NOT NULL CHECK (density IN ('dense', 'sparse', 'typical')),
    radius_ly             DOUBLE NOT NULL,
    composition_summary   TEXT NOT NULL,
    galactic_orbital_period_gy           DOUBLE NOT NULL,
    galactic_orbital_phase_deg           DOUBLE NOT NULL,
    galactic_min_update_interval_years   DOUBLE NOT NULL,
    epoch_unix                  DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due             DOUBLE,

    -- v18: this field's own galaxy-frame center -- see schema.sql's "v18"
    -- header note. NULL together: never placed in the galaxy.
    center_x_pc         DOUBLE,
    center_y_pc         DOUBLE,
    center_z_pc         DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    velocity_x_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_y_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    velocity_z_kms              DOUBLE NULL,  -- v102 (DB.55): galactic-axes velocity, km/s
    mass_kg                     DOUBLE NULL,  -- v102 (DB.55): point mass, a stated zero where the mass is not known
    formed_at                   DOUBLE NULL,  -- v102 (DB.55)
    origin_event_id             BIGINT UNSIGNED NULL,  -- v102 (DB.55)
    KEY idx_asteroid_fields_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    galactic_radius_pc  DOUBLE,
    -- v41: the sector octant the center sits in (star_systems.quadrant's labels).
    quadrant            VARCHAR(4) CHECK (quadrant IN ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')),

    -- v27: row timestamps -- see the header comment's "v27" note.
    -- v39: the innermost nebula or supernova remnant this sits inside
    -- (at most one set; _db.refresh_containment).
    inside_nebula_id    BIGINT UNSIGNED,
    inside_remnant_id   BIGINT UNSIGNED,

    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    chondrite_group     VARCHAR(4) NULL,  -- v190 (GEN.475): the nearest chondrite group of the field's mean composition, `none` far from every group, NULL before v190
    -- v58: this object's unique ID, see the header comment's "v58" note. NULL for a row saved before v58.
    uid                   BINARY(12),

    -- Named explicitly -- see `nebulae`'s identical "v18" CHECK comment above.
    UNIQUE KEY uq_asteroid_fields_uid (uid),
    KEY idx_asteroid_fields_name (name),  -- v178 (PERF.205)
    KEY idx_asteroid_fields_placed_name (((center_x_pc IS NOT NULL)), name),   -- v186 (PERF.211): the Phenomena table's sorts read in index order
    KEY idx_asteroid_fields_descriptor_name (density, name),
    KEY idx_asteroid_fields_radius_name (radius_ly, name),
    CONSTRAINT chk_asteroid_fields_placement CHECK (
        (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
        (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
        (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
    ),

    CONSTRAINT fk_asteroid_fields_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE SET NULL,
    KEY idx_asteroid_fields_sector_id (sector_id),
    KEY idx_asteroid_fields_next_update_due (next_update_due),
    KEY idx_asteroid_fields_galactic_radius_pc (galactic_radius_pc),
    -- v26: see black_holes' identical "v26" index comment above.
    KEY idx_asteroid_fields_center (center_x_pc, center_y_pc, center_z_pc),
    -- v27: see the header comment's "v27" note.
    KEY idx_asteroid_fields_modified_at (modified_at),
    -- v39: containment -- see the header comment's "v39" note.
    CONSTRAINT fk_asteroid_fields_inside_nebula
        FOREIGN KEY (inside_nebula_id) REFERENCES nebulae(id) ON DELETE SET NULL,
    CONSTRAINT fk_asteroid_fields_inside_remnant
        FOREIGN KEY (inside_remnant_id) REFERENCES supernova_remnants(id) ON DELETE SET NULL
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
-- sector_name_registry / system_name_registry -- name-uniqueness
-- bookkeeping (v24, see the header comment's "v24" note and
-- `planetgen/names/uniqueness.py`'s own module docstring for the
-- sector > system decoration hierarchy; v34 dropped the planet/moon
-- registry, since those names now derive from the system's). One row per
-- distinct base name (the name with every decoration this project could
-- have added stripped back off -- `nameUniqueness.strip_decoration`)
-- that has ever collided at least once; a base name that's only ever
-- been used by a single row anywhere has no row here at all.
--
-- `first_*_id` always points at the row that *first* used this base
-- name -- the one that gets renamed as collisions happen (bare -> Alpha,
-- later Alpha -> Alpha ... I). It never changes to a different row once
-- set. `ON DELETE CASCADE` means a deleted "first"
-- row's bookkeeping simply disappears -- correct, since the base name is
-- genuinely free again once nothing live still uses any decorated form
-- of it descended from that row.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS sector_name_registry (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    base_name         VARCHAR(255) NOT NULL,
    occurrence_count  INT NOT NULL,
    first_sector_id   BIGINT UNSIGNED NOT NULL,

    UNIQUE (base_name),
    CONSTRAINT fk_sector_name_registry_first_sector
        FOREIGN KEY (first_sector_id) REFERENCES sectors(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS system_name_registry (
    id                     BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    base_name              VARCHAR(255) NOT NULL,
    occurrence_count       INT NOT NULL,
    -- The first holder of the name: a star system, or (v40) a uniquely
    -- named phenomenon -- `first_object_table`/`first_object_id` (e.g.
    -- 'nebulae', 12), with first_star_system_id NULL. A phenomenon holder
    -- has no foreign key; deleting it leaves the count one high, which
    -- only means the next holder gets a decoration it could have skipped.
    first_star_system_id   BIGINT UNSIGNED,
    first_object_table     VARCHAR(32),
    first_object_id        BIGINT UNSIGNED,

    -- Index into names.DIMINUTIVE_PREFIXES already applied to this base
    -- name's system side against a colliding sector -- NULL until the
    -- first such cross-level hit (nameUniqueness.resolve_diminutive).
    diminutive_index       INT,

    UNIQUE (base_name),
    CONSTRAINT fk_system_name_registry_first_star_system
        FOREIGN KEY (first_star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- nearest_systems (v41): the 3 star systems nearest each placed star
-- system and phenomenon, searched across sector boundaries in the galaxy
-- frame (`_db.refresh_nearest_systems`). `object_table`/`object_id` name
-- the object (`'star_systems'` or a phenomenon table); `star_system_id`
-- repeats `object_id` for a system so deleting it takes its rows along,
-- and `sector_id` is the object's sector, so deleting a sector does too.
-- Only systems within `_db.NEAREST_SYSTEMS_SEARCH_PC` count, so an
-- isolated object can have fewer than 3 rows. Filled when a sector is
-- generated (a new sector also updates its neighbors' lists) and by the
-- correlative update (`planetgen.cli.orbits`), which moves everything.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS nearest_systems (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id            BIGINT UNSIGNED NOT NULL,
    object_table         VARCHAR(32) NOT NULL,
    object_id            BIGINT UNSIGNED NOT NULL,
    star_system_id       BIGINT UNSIGNED,
    neighbor_rank        TINYINT NOT NULL,
    neighbor_system_id   BIGINT UNSIGNED NOT NULL,
    distance_pc          DOUBLE NOT NULL,

    UNIQUE KEY uq_nearest_systems_object_rank (object_table, object_id, neighbor_rank),
    KEY idx_nearest_systems_neighbor (neighbor_system_id),
    CONSTRAINT fk_nearest_systems_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    CONSTRAINT fk_nearest_systems_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_nearest_systems_neighbor
        FOREIGN KEY (neighbor_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- bright_stars (v43): every star at least
-- `galaxy_shape.bright_star_min_luminosity_sol` bright, generated and
-- placed galaxy-wide by `planetgen plan`'s scatter step before any
-- sector is filled. Each row is a finished star (no planets or companion
-- yet, and no name: it's named when its system is made) at a galaxy-frame
-- position (milliparsecs) inside its sector's cell. Filling that sector
-- builds a full system around it and sets `star_system_id`; the row is
-- kept as the Galaxy Map's index of bright stars. A plan re-run or reset
-- empties the table. See
-- /mnt/project-files/galaxy-studies/bright-star-preplacement-plan.md.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS bright_stars (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    shell_index           INT NOT NULL,
    band_index            INT NOT NULL,
    wedge_index           INT NOT NULL,
    position_x_mpc       DOUBLE NOT NULL,
    position_y_mpc       DOUBLE NOT NULL,
    position_z_mpc       DOUBLE NOT NULL,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(position_x_mpc / 1000.0) < 262144 AND ABS(position_y_mpc / 1000.0) < 262144 AND ABS(position_z_mpc / 1000.0) < 262144, ((((((((((((LEAST(CAST(FLOOR(((position_x_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((position_x_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((position_x_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((position_y_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((position_y_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((position_y_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((position_z_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((position_z_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((position_z_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    KEY idx_bright_stars_cell (cell_key, position_x_mpc, position_y_mpc, position_z_mpc),
    population           VARCHAR(24) NOT NULL,
    star_type            VARCHAR(64) NOT NULL,
    yerkes_class         VARCHAR(16) NOT NULL,
    mass_kg              DOUBLE NOT NULL,
    radius_km            DOUBLE NOT NULL,
    temperature_k        DOUBLE NOT NULL,
    luminosity_w         DOUBLE NOT NULL,
    age_gy               DOUBLE NOT NULL,
    lifespan_gy          DOUBLE,
    initial_mass_sol     DOUBLE NOT NULL,
    phase_end_age_gy     DOUBLE,
    seed                 BIGINT UNSIGNED NOT NULL,
    star_system_id       BIGINT UNSIGNED,
    created_at           TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    off_plane            TINYINT AS (ABS(position_z_mpc) >= 250000) VIRTUAL,
    velocity_x_kms       DOUBLE,  -- v93 (GEN.223): the star's point in space, see the header comment
    velocity_y_kms       DOUBLE,
    velocity_z_kms       DOUBLE,
    epoch_unix           DOUBLE,
    next_update_due      DOUBLE,
    member_count         SMALLINT UNSIGNED,  -- v93: companions, see the header comment
    system_mass_solar    DOUBLE,
    metallicity_fe_h     DOUBLE,  -- v95 (GEN.158): [Fe/H] in dex, NULL for a cluster star until GEN.159 and for older rows
    cluster_id           BIGINT UNSIGNED,  -- v127 (DB.61): the cluster this star belongs to, NULL for a field star
    rotation_period_hours DOUBLE,  -- v185 (GEN.846): the spin period in hours, NULL on a row made before v185
    nuclear_arrival_age_gy DOUBLE,  -- v197 (GEN.1174): the age MIST's tracks reach the main sequence, see the header comment

    KEY idx_bright_stars_address (shell_index, band_index, wedge_index),
    KEY idx_bright_stars_luminosity (luminosity_w),
    KEY idx_bright_stars_next_update_due (next_update_due),
    KEY idx_bright_stars_population (population, off_plane, luminosity_w),
    CONSTRAINT chk_bright_stars_population CHECK (population IN ('young', 'intermediate', 'old', 'bulge', 'halo', 'cluster', 'satellite', 'nuclear_cluster', 'nuclear_disc', 'satellite_young', 'satellite_intermediate', 'open_cluster_00', 'open_cluster_01', 'open_cluster_02', 'open_cluster_03', 'open_cluster_04', 'open_cluster_05', 'open_cluster_06', 'open_cluster_07', 'open_cluster_08', 'open_cluster_09', 'open_cluster_10', 'open_cluster_11')),
    KEY idx_bright_stars_cluster (cluster_id),
    CONSTRAINT fk_bright_stars_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE SET NULL,
    CONSTRAINT fk_bright_stars_cluster
        FOREIGN KEY (cluster_id) REFERENCES galaxy_components(id) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- ---------------------------------------------------------------------
-- sector_stats (v53, GEN.44 and PERF.11): one row per sector address
-- (shell, band, wedge) that a backfill reached or a fill generated --
-- unfilled sectors have no `sectors` row, so this is keyed by address.
--
-- `bright_level_sol` is how deep the sector's bright stars go (GEN.44):
-- -1 untouched (it follows `galaxy_shape.bright_star_min_luminosity_sol`,
-- and holding stars with no level there means a run failed part way, so
-- they are wiped and drawn again), a positive L_sun for the dimmest a
-- backfill drew every star down to, 0 once the sector is generated. It is
-- written in the same transaction as the stars it describes. A fill keeps
-- the level it found in `level_before_fill_sol`, which deleting the
-- sector puts back. A new scatter sets every positive level back to -1.
--
-- `relative_density` and `expected_systems` are the galaxy model's
-- density at the sector's center and the systems expected there
-- (PERF.11); `actual_systems`, `actual_stars` and the mean temperature
-- (K) and luminosity (L_sun) of its stars are what a fill gave it, NULL
-- before. A fill also keeps the facts the Galaxy Map colors from (v55,
-- DB.14): `mean_age_gy`, the mean age of its stars, and
-- `total_luminosity_sol`, their summed luminosity (both NULL for a
-- sector with no stars). The map works the color out from these; none
-- is baked into the row.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS sector_stats (
    shell_index           INT NOT NULL,
    band_index            INT NOT NULL,
    wedge_index           INT NOT NULL,
    bright_level_sol       DOUBLE NOT NULL DEFAULT -1,
    level_before_fill_sol  DOUBLE,
    -- v79 (GEN.187): see header comment.
    bright_mass_sol        DOUBLE,
    relative_density       DOUBLE,
    expected_systems       DOUBLE,
    actual_systems         INT,
    actual_stars           INT,
    mean_temperature_k     DOUBLE,
    mean_luminosity_sol    DOUBLE,
    mean_age_gy            DOUBLE,
    total_luminosity_sol   DOUBLE,
    habitable_worlds       INT,   -- v139 (MAP.143): planets and moons in the sector's systems that are Blue or Green overall (PHI-6)
    total_mass_kg          DOUBLE,   -- v168 (DB.108): the mass of the sector's stars, planets, moons, belts, comets and phenomena (not nebulae), NULL when a body has no stored mass
    filled_at              TIMESTAMP(3) NULL,
    updated_at             TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,

    PRIMARY KEY (shell_index, band_index, wedge_index),
    KEY idx_sector_stats_level (bright_level_sol)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- ---------------------------------------------------------------------
-- sector_paths and sector_path_knots (v62, GEN.123): the path a body with
-- no closed orbit takes through its sector (`physics/sector_path.py`), as
-- cubic Hermite spline knots. One `sector_paths` row per body
-- (`object_table`/`object_id`: a star system, rogue planet or interstellar
-- comet), in the sector that holds it, with `exited` saying whether the
-- body leaves the sector (rather than the path being cut off) and
-- `duration_years` how long the crossing takes; its `sector_path_knots`
-- hold the spline: `t_years` after the body's current position and
-- velocity, position in galactic parsecs, velocity in km/s. There is no
-- foreign key to the body, so a deleted body's path lives until its
-- sector's paths are next rewritten (`db.sector_paths.compute_sector_paths`).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS sector_paths (
    id              BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    sector_id       BIGINT UNSIGNED NOT NULL,
    object_table    VARCHAR(24) NOT NULL CHECK (object_table IN ('star_systems', 'rogue_planets', 'interstellar_comets')),
    object_id       BIGINT UNSIGNED NOT NULL,
    exited          TINYINT(1) NOT NULL,
    duration_years  DOUBLE NOT NULL,
    computed_at     TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    UNIQUE (object_table, object_id),
    KEY idx_sector_paths_sector (sector_id),
    CONSTRAINT fk_sector_paths_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS sector_path_knots (
    id        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    path_id   BIGINT UNSIGNED NOT NULL,
    position  INT NOT NULL,
    t_years   DOUBLE NOT NULL,
    x_pc      DOUBLE NOT NULL,
    y_pc      DOUBLE NOT NULL,
    z_pc      DOUBLE NOT NULL,
    vx_kms    DOUBLE NOT NULL,
    vy_kms    DOUBLE NOT NULL,
    vz_kms    DOUBLE NOT NULL,

    UNIQUE (path_id, position),
    CONSTRAINT fk_sector_path_knots_path
        FOREIGN KEY (path_id) REFERENCES sector_paths(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- ---------------------------------------------------------------------
-- generation_runs (v52, DB.6): one row per `planetgen` run that changes
-- the galaxy, written when it starts and finished when it ends: its
-- subcommand and command line (one `generation_run_arguments` row per
-- argument, without the --mysql-* and --debug options), the run's own 128-bit seed (what a random start or a
-- one-off system draws from) and the galaxy seed it ran against, the
-- code's version key and versions, and the outcome (`ok`, `failed`,
-- `interrupted`; NULL while it runs or when it died without saying).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS generation_runs (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    command              VARCHAR(32) NOT NULL,
    run_seed             BINARY(16),
    galaxy_seed          BINARY(16),
    version_key          CHAR(22) NOT NULL,
    planetgen_version    VARCHAR(32) NOT NULL,
    python_version       VARCHAR(16) NOT NULL,
    platform             VARCHAR(64) NOT NULL,
    started_at           TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    finished_at          TIMESTAMP(3) NULL,
    outcome              VARCHAR(16),
    generator_epoch      SMALLINT UNSIGNED NULL,  -- v131 (OPS.58): the epoch the run's code generated under, NULL before v131

    KEY idx_generation_runs_started_at (started_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v63 (DB.13): a run's command line, one row per argument in order (it was
-- one JSON text column).
CREATE TABLE IF NOT EXISTS generation_run_arguments (
    run_id       BIGINT UNSIGNED NOT NULL,
    position     INT NOT NULL,
    value        VARCHAR(1024) NOT NULL,
    PRIMARY KEY (run_id, position),
    CONSTRAINT fk_generation_run_arguments_run FOREIGN KEY (run_id) REFERENCES generation_runs(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- v164 (GEN.764): the generator epoch each kind of scatter was drawn under. A scatter is galaxy-wide, and its rows carry no epoch of their own, so
-- the draw that finishes a kind writes its row and a clear deletes it; the sectors (v121), the plan (v131) and the species (v164) carry theirs.
CREATE TABLE IF NOT EXISTS scatter_epochs (
    kind             VARCHAR(16) NOT NULL PRIMARY KEY,
    generator_epoch  SMALLINT UNSIGNED NULL,
    run_id           BIGINT UNSIGNED NULL,
    recorded_at      TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),

    CONSTRAINT chk_scatter_epochs_kind CHECK (kind IN ('mass', 'luminosity', 'band', 'phenomena', 'nebula-field')),
    CONSTRAINT fk_scatter_epochs_run FOREIGN KEY (run_id) REFERENCES generation_runs(id) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- ---------------------------------------------------------------------
-- phenomenon_scatter (v64, GEN.100): every black hole, neutron star,
-- planetary nebula and supernova remnant of the galaxy, its hypervelocity
-- stars and its nucleus, placed galaxy-wide by `planetgen plan` just
-- before the bright-star scatter and before any sector is filled. A row is
-- a point (milliparsecs, galaxy frame) in its sector's cell and the 63-bit
-- seed the object is built from; filling the sector builds the object at
-- that point and sets `built_at`, so what was placed stays where it was
-- put. `subtype` is `supermassive` for the nucleus black hole and NULL for
-- the rest (the object draws its own class from `seed`); the velocity
-- columns hold a hypervelocity star's galactic velocity (km/s) and are
-- NULL otherwise. A plan re-run or reset empties the table.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS phenomenon_scatter (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    shell_index           INT NOT NULL,
    band_index            INT NOT NULL,
    wedge_index           INT NOT NULL,
    kind                 VARCHAR(24) NOT NULL,
    subtype              VARCHAR(16),
    position_x_mpc       DOUBLE NOT NULL,
    position_y_mpc       DOUBLE NOT NULL,
    position_z_mpc       DOUBLE NOT NULL,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(position_x_mpc / 1000.0) < 262144 AND ABS(position_y_mpc / 1000.0) < 262144 AND ABS(position_z_mpc / 1000.0) < 262144, ((((((((((((LEAST(CAST(FLOOR(((position_x_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((position_x_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((position_x_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((position_y_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((position_y_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((position_y_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((position_z_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((position_z_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((position_z_mpc / 1000.0) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    KEY idx_phenomenon_scatter_cell (cell_key, position_x_mpc, position_y_mpc, position_z_mpc),
    velocity_x_kms       DOUBLE,
    velocity_y_kms       DOUBLE,
    velocity_z_kms       DOUBLE,
    seed                 BIGINT UNSIGNED NOT NULL,
    built_at             TIMESTAMP NULL,
    epoch_unix           DOUBLE,  -- v69 (GEN.137): when the position holds, see the header comment
    mass_solar           DOUBLE,  -- v80 (MAP.165): see header comment
    map_luminosity_sol   FLOAT,   -- v84 (DB.23): see header comment
    next_update_due      DOUBLE,  -- v93 (GEN.223): see header comment
    member_count         SMALLINT UNSIGNED,  -- v93: companions, see the header comment
    system_mass_solar    DOUBLE,

    KEY idx_phenomenon_scatter_address (shell_index, band_index, wedge_index),
    KEY idx_phenomenon_scatter_next_update_due (next_update_due),
    KEY idx_phenomenon_scatter_class (kind, subtype, mass_solar),  -- v82 (PERF.68): see header comment
    CONSTRAINT chk_phenomenon_scatter_kind CHECK (kind IN (
        'black-hole', 'neutron-star', 'planetary-nebula', 'supernova-remnant', 'hypervelocity-star', 'quasar'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Companions of a scatter point (v93, docs/design/systems/multi-star-systems.md section 5): the point is the system's
-- barycentre, `member_count` and `system_mass_solar` say it has members, and these two tables hold them, one column
-- per value. `source` 1 is `bright_stars`, 2 `phenomenon_scatter`; no foreign key (the parent is polymorphic), so
-- the build deletes a built row's members and nodes in its own transaction.
CREATE TABLE IF NOT EXISTS scatter_members (
    source               TINYINT NOT NULL,
    source_id            BIGINT UNSIGNED NOT NULL,
    member_no            SMALLINT NOT NULL,
    parent_node_no       SMALLINT NOT NULL,
    slot                 SMALLINT NOT NULL,
    kind                 VARCHAR(20) NOT NULL,
    state                VARCHAR(16),
    flags                SMALLINT UNSIGNED,
    offset_x_au          DOUBLE,
    offset_y_au          DOUBLE,
    offset_z_au          DOUBLE,
    velocity_x_kms       DOUBLE,
    velocity_y_kms       DOUBLE,
    velocity_z_kms       DOUBLE,
    decay_age_yr         DOUBLE,
    mass_solar           DOUBLE NOT NULL,
    initial_mass_solar   DOUBLE NOT NULL,
    seed                 BIGINT UNSIGNED NOT NULL,
    epoch_unix                  DOUBLE NULL,  -- v102 (DB.55): when the position and velocity hold
    radius_au                   DOUBLE NULL,  -- v184 (GEN.952): the extent of the member in AU, NULL on a row made before v184
    radius_max_au               DOUBLE NULL,  -- v189 (GEN.1009): the largest radius the member has had, in AU, NULL on a row made before v187
    mass_max_solar              DOUBLE NULL,  -- v189 (GEN.1009): the largest mass the member has had, in solar masses, NULL on a row made before v187
    PRIMARY KEY (source, source_id, member_no),
    CONSTRAINT chk_scatter_members_source CHECK (source IN (1, 2)),
    CONSTRAINT chk_scatter_members_kind CHECK (kind IN (
        'main-sequence', 'giant', 'white-dwarf', 'neutron-star', 'black-hole', 'brown-dwarf',
        'pre-main-sequence', 'protostar', 'stripped-helium', 'thorne-zytkow', 'wolf-rayet'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS scatter_nodes (
    source               TINYINT NOT NULL,
    source_id            BIGINT UNSIGNED NOT NULL,
    node_no              SMALLINT NOT NULL,
    parent_node_no       SMALLINT,
    slot                 SMALLINT,
    mass_solar           DOUBLE NOT NULL,
    separation_km        DOUBLE NOT NULL,
    eccentricity         DOUBLE NOT NULL,
    inclination_deg      DOUBLE NOT NULL,
    ascending_node_deg   DOUBLE NOT NULL,
    arg_periapsis_deg    DOUBLE NOT NULL,
    mean_anomaly_deg     DOUBLE NOT NULL,
    epoch_unix           DOUBLE NOT NULL,
    period_years         DOUBLE NOT NULL,
    kind                 VARCHAR(16) NOT NULL,
    interaction          VARCHAR(20),
    flags                TINYINT UNSIGNED NOT NULL DEFAULT 0,
    PRIMARY KEY (source, source_id, node_no),
    CONSTRAINT chk_scatter_nodes_source CHECK (source IN (1, 2)),
    CONSTRAINT chk_scatter_nodes_kind CHECK (kind IN ('close', 'wide', 'compact', 'group')),
    CONSTRAINT chk_scatter_nodes_interaction CHECK (interaction IS NULL OR interaction IN (
        'detached', 'semi-detached', 'contact', 'accreting', 'post-common-envelope', 'merger-product'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- ---------------------------------------------------------------------
-- orbit_nodes (v104, DB.53): the tree of orbits of a built multi-star system, one row per node, with the same names and
-- units as `scatter_nodes` so a fill copies column for column. A node is one orbit (kind `close`, `wide`, `compact` or
-- `group`) of the masses below it; `parent_node_no` is the node it orbits in (NULL for the root), `slot` its place among
-- its parent's children. A star names its node in `stars.parent_node_id` (the `node_no` here), and a circumbinary planet
-- in `planets.host_node_id`. No foreign key joins them: a node is deleted with its system. `disk_*` is a disk of gas and
-- dust about the node's centre (no bodies in it).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS orbit_nodes (
    star_system_id       BIGINT UNSIGNED NOT NULL,
    node_no              SMALLINT NOT NULL,
    parent_node_no       SMALLINT,
    slot                 SMALLINT,
    mass_solar           DOUBLE NOT NULL,
    separation_km        DOUBLE NOT NULL,
    eccentricity         DOUBLE NOT NULL,
    inclination_deg      DOUBLE NOT NULL,
    ascending_node_deg   DOUBLE NOT NULL,
    arg_periapsis_deg    DOUBLE NOT NULL,
    mean_anomaly_deg     DOUBLE NOT NULL,
    epoch_unix           DOUBLE NOT NULL,
    period_years         DOUBLE NOT NULL,
    kind                 VARCHAR(16) NOT NULL,
    interaction          VARCHAR(20),
    flags                TINYINT UNSIGNED NOT NULL DEFAULT 0,
    disk_kind            VARCHAR(16),
    disk_mass_solar      DOUBLE,
    disk_inner_au        DOUBLE,
    disk_outer_au        DOUBLE,
    PRIMARY KEY (star_system_id, node_no),
    CONSTRAINT fk_orbit_nodes_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT chk_orbit_nodes_kind CHECK (kind IN ('close', 'wide', 'compact', 'group')),
    CONSTRAINT chk_orbit_nodes_interaction CHECK (interaction IS NULL OR interaction IN (
        'detached', 'semi-detached', 'contact', 'accreting', 'post-common-envelope', 'merger-product'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- ---------------------------------------------------------------------
-- dark_members (v194, GEN.1013): the members of a built system's tree that are not stars of it, one row each: a brown dwarf, a
-- neutron star, a black hole or a stripped helium star. `parent_node_no` and `slot` name the node it hangs in and its place among
-- that node's children (the numbering of `orbit_nodes`); the node masses above it count its `mass_solar`. It has no light and no
-- planets. A row goes with its system.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS dark_members (
    star_system_id       BIGINT UNSIGNED NOT NULL,
    member_no            SMALLINT NOT NULL,
    parent_node_no       SMALLINT NOT NULL,
    slot                 SMALLINT NOT NULL,
    kind                 VARCHAR(16) NOT NULL,
    mass_solar           DOUBLE NOT NULL,
    radius_au            DOUBLE NULL,
    PRIMARY KEY (star_system_id, member_no),
    CONSTRAINT fk_dark_members_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT chk_dark_members_kind CHECK (kind IN ('brown-dwarf', 'neutron-star', 'black-hole', 'stripped-helium'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- ---------------------------------------------------------------------
-- sector_system_counts (v82, PERF.74; v181, PERF.208): each sector's system count, system density and Quadrant as of
-- the last refresh (`query.refresh_sector_system_counts`), one row for each sector (a sector made since has none yet),
-- so the Sectors list reads these instead of counting every system per request and its sorts by systems, density and
-- position page by key on the three indexes.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS sector_system_counts (
    sector_id      BIGINT UNSIGNED NOT NULL PRIMARY KEY,
    system_count   INT UNSIGNED NOT NULL,
    system_density DOUBLE NULL,         -- v181: systems per cubic milliparsec of the sector, NULL for a sector with no galaxy position
    quadrant       VARCHAR(8) NULL,     -- v181: 'I' to 'IV' as `spherical.sector_quadrant` says, NULL for a sector with no galaxy position
    KEY idx_sector_system_counts_count (system_count, sector_id),
    KEY idx_sector_system_counts_density (system_density, sector_id),
    KEY idx_sector_system_counts_quadrant (quadrant, sector_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- system_interests (v165, MAP.305, MAP.303): one row for each reason a star system is worth a mark on the Galaxy Map
-- (`planetgen.db.interests.REASONS`: a planet or moon rated Yellow or better overall on the PHI-6 score, a scored barren one, a gas giant, a
-- belt, and more as they are added), with how many bodies give it (`n`), the largest score among them (`best_score`, NULL for a reason
-- with no score; the `civilization` reason's is a technology level, v174) and where the system is, so the map's marker layer finds the systems in a box (`query.galaxy_system_interests_in_box`)
-- without reading every planet. Derived data: `interests.record` rebuilds a system's rows when its sector is filled and whenever its
-- scores are worked out again; the position is the system's `galaxy_*_pc`, which never changes. It replaces `habitable_systems` (v148).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS system_interests (
    star_system_id    BIGINT UNSIGNED NOT NULL,
    reason            VARCHAR(24) NOT NULL,
    n                 SMALLINT UNSIGNED NOT NULL,
    best_score        DOUBLE NULL,
    galaxy_x_pc       DOUBLE NOT NULL,
    galaxy_y_pc       DOUBLE NOT NULL,
    galaxy_z_pc       DOUBLE NOT NULL,
    PRIMARY KEY (star_system_id, reason),
    KEY idx_system_interests_position (galaxy_x_pc, galaxy_y_pc, galaxy_z_pc),
    CONSTRAINT fk_system_interests_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- nebula_field (v83, DB.24): the galaxy's field clouds of classes M and N
-- (giant molecular clouds and dark clouds; P and Q stay implicit), one row
-- each, made by `planetgen plan` from the galaxy seed. (`cell_i`, `cell_j`,
-- `cell_k`) is the 50 pc field cell and `n` the cloud's place in that cell's
-- draw (`nebula_field.cell_draw`), so a row is found again from the seed
-- alone. `uid` and `is_stored` stay empty until a sector's save stores the cloud
-- as a `nebulae` row, which gives it its object ID (working it out for every
-- cloud costs a nebula shape each); `idx_nebula_field_uid` finds it by that ID.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS nebula_field (
    cell_i        INT NOT NULL,
    cell_j        INT NOT NULL,
    cell_k        INT NOT NULL,
    n             INT UNSIGNED NOT NULL,
    nebula_class  CHAR(1) NOT NULL,
    radius_ly     DOUBLE NOT NULL,
    center_x_pc   DOUBLE NOT NULL,
    center_y_pc   DOUBLE NOT NULL,
    center_z_pc   DOUBLE NOT NULL,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    KEY idx_nebula_field_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    uid           BINARY(12),
    is_stored     TINYINT(1) NOT NULL DEFAULT 0,
    PRIMARY KEY (cell_i, cell_j, cell_k, n),
    KEY idx_nebula_field_uid (uid)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- nebula_block_stats (v116, MAP.172): what the nebulae add up to in one block of the Galaxy Map's block
-- ladder, `(m, shell_index, band_index, wedge_index)` as `galaxy_block_stats` names a block, for the five
-- levels above a sector. A nebula counts in the block that holds its centre (`nebula_count`, `volume_pc3`,
-- the volume-weighted `centroid_*_pc`, `bound_radius_pc` around it, and the volume of each family in
-- `volume_<family>_pc3`; the class with the most volume is `dominant_class`); the part of its sphere that
-- lies in a block is what `cover_volume_pc3` adds up, so a block a cloud only reaches has a row with a count
-- of 0, and `cover` = 1 - exp(-cover_volume_pc3 / the block's volume), the share of the block inside a
-- nebula. The nebulae are the rows of `nebula_field` and the stored `nebulae` that are not field clouds
-- (matched by `uid`). No foreign keys: a block is not a row of its own.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS nebula_block_stats (
    m                         SMALLINT UNSIGNED NOT NULL,
    shell_index               INT NOT NULL,
    band_index                INT NOT NULL,
    wedge_index               INT NOT NULL,
    nebula_count              INT UNSIGNED NOT NULL,
    volume_pc3                DOUBLE NOT NULL,
    centroid_x_pc             DOUBLE NULL,
    centroid_y_pc             DOUBLE NULL,
    centroid_z_pc             DOUBLE NULL,
    bound_radius_pc           DOUBLE NULL,
    volume_diffuse_pc3        DOUBLE NOT NULL,
    volume_emission_pc3       DOUBLE NOT NULL,
    volume_reflection_pc3     DOUBLE NOT NULL,
    volume_planetary_pc3      DOUBLE NOT NULL,
    volume_dark_pc3           DOUBLE NOT NULL,
    volume_remnant_pc3        DOUBLE NOT NULL,
    dominant_class            CHAR(1) NULL,
    cover_volume_pc3          DOUBLE NOT NULL,
    cover                     DOUBLE NOT NULL,
    PRIMARY KEY (m, shell_index, band_index, wedge_index)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- nebula_block_classes (v116, MAP.172): the volume of each nebula class in a block of `nebula_block_stats`
-- (nebulae by the block of their centre), so a stored nebula can be added to a block and the dominant class
-- found again without reading the nebulae.
CREATE TABLE IF NOT EXISTS nebula_block_classes (
    m                         SMALLINT UNSIGNED NOT NULL,
    shell_index               INT NOT NULL,
    band_index                INT NOT NULL,
    wedge_index               INT NOT NULL,
    nebula_class              CHAR(1) NOT NULL,
    volume_pc3                DOUBLE NOT NULL,
    PRIMARY KEY (m, shell_index, band_index, wedge_index, nebula_class)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- nebula_block_state (v116, MAP.172): one row, `nebulae_through_id`, the highest `nebulae.id` the two tables above have
-- folded in (`db/nebula_blocks.catch_up` adds the nebulae stored after it). Absent until the tables are first built.
CREATE TABLE IF NOT EXISTS nebula_block_state (
    id                        TINYINT UNSIGNED NOT NULL,
    nebulae_through_id        BIGINT UNSIGNED NOT NULL,
    PRIMARY KEY (id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- orbital_events (v118, DB.56): one row per orbital event. `kind` is `collision` (an outcome of the table in
-- docs/design/orbits/collisions-and-mergers.md section 5.1), `merger`, `split`, `ejection`, `capture`, `unbind` or `warning`;
-- `outcome_code` is the Record column of that table (`merger`, `hit_and_run`, `erosion`, `disruption`, `absorbed`,
-- `tidal_disruption`, `giant_merger`, `brown_dwarf_formed`, `star_created`, `stellar_merger`) and `regime` the impact
-- regime. `event_time_unix` is when it happened (t_hit), `detected_time_unix` when the update found it, `seed_string` the
-- string its draws were seeded with, so the outcome can be replayed. Units: km, km/s, kg, J/kg, seconds; a column an event
-- kind has no number for is NULL. The numbers (separations, b and b_crit, v_inf, v_esc, v_imp, Q_R, Q*_RD, x, the
-- interacting-mass fraction, the pair's catastrophic velocity, the window and whether it was capped, whether a focusing
-- radius was used), the conservation check (`mass_error` and `momentum_error` relative, against `conservation_tolerance`,
-- the impact kinetic energy and the share assumed deposited as heat), the products (`lost_mass_kg`, `merged_mass_kg`,
-- `new_radius_km`, `spin_breakup_ratio`, the fusion lines used and `borderline` within 3 percent of one) and the context
-- (`containing_uid` the nebula or remnant the event happened inside, the bodies inside the influence radius, the sectors
-- marked dirty, suppressed duplicates, a pair already overlapping at the start of the window, `numeric_warning`, and
-- `flags`, a comma-separated list of the flags an outcome carries). `run_id` is the `generation_runs` row that made it
-- (no foreign key: a run row may be purged and the event stays).
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS orbital_events (
    id                          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    run_id                      BIGINT UNSIGNED NULL,
    kind                        VARCHAR(16) NOT NULL,
    outcome_code                VARCHAR(24) NULL,
    regime                      VARCHAR(24) NULL,
    event_time_unix             DOUBLE NOT NULL,
    detected_time_unix          DOUBLE NOT NULL,
    seed_string                 VARCHAR(255) NOT NULL,
    separation_start_km         DOUBLE NULL,
    separation_impact_km        DOUBLE NULL,
    impact_parameter_km         DOUBLE NULL,
    critical_impact_parameter_km DOUBLE NULL,
    v_inf_kms                   DOUBLE NULL,
    v_esc_kms                   DOUBLE NULL,
    v_imp_kms                   DOUBLE NULL,
    q_r_j_per_kg                DOUBLE NULL,
    q_star_rd_j_per_kg          DOUBLE NULL,
    erosion_x                   DOUBLE NULL,
    interacting_mass_fraction   DOUBLE NULL,
    catastrophic_velocity_kms   DOUBLE NULL,
    window_s                    DOUBLE NULL,
    window_capped               TINYINT NOT NULL DEFAULT 0,
    focusing_radius_used        TINYINT NOT NULL DEFAULT 0,
    mass_error                  DOUBLE NULL,
    momentum_error              DOUBLE NULL,
    conservation_tolerance      DOUBLE NULL,
    impact_energy_j             DOUBLE NULL,
    heat_share                  DOUBLE NULL,
    lost_mass_kg                DOUBLE NULL,
    merged_mass_kg              DOUBLE NULL,
    new_radius_km               DOUBLE NULL,
    spin_breakup_ratio          DOUBLE NULL,
    deuterium_limit_mj          DOUBLE NULL,
    hydrogen_limit_sol          DOUBLE NULL,
    borderline                  TINYINT NOT NULL DEFAULT 0,
    containing_uid              BINARY(12) NULL,
    influence_body_count        INT UNSIGNED NULL,
    sectors_dirtied             INT UNSIGNED NULL,
    suppressed_duplicates       INT UNSIGNED NOT NULL DEFAULT 0,
    already_overlapping         TINYINT NOT NULL DEFAULT 0,
    numeric_warning             TINYINT NOT NULL DEFAULT 0,
    flags                       VARCHAR(255) NOT NULL DEFAULT '',
    recorded_at                 TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    KEY idx_orbital_events_time (event_time_unix, id),
    KEY idx_orbital_events_kind (kind, event_time_unix, id),
    KEY idx_orbital_events_run (run_id, id),
    CONSTRAINT chk_orbital_events_kind CHECK (kind IN
        ('collision', 'merger', 'split', 'ejection', 'capture', 'unbind', 'warning')),
    CONSTRAINT chk_orbital_events_window_capped CHECK (window_capped IN (0, 1)),
    CONSTRAINT chk_orbital_events_focusing CHECK (focusing_radius_used IN (0, 1)),
    CONSTRAINT chk_orbital_events_borderline CHECK (borderline IN (0, 1)),
    CONSTRAINT chk_orbital_events_overlapping CHECK (already_overlapping IN (0, 1)),
    CONSTRAINT chk_orbital_events_numeric_warning CHECK (numeric_warning IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- orbital_event_bodies (v118, DB.56): the bodies of an event, one row each. `role` is `participant` (a body that went into it)
-- or `product` (one that came out: a survivor, a field, a merged giant, a star), `ordinal` the order within the role, `uid` the
-- 96-bit object ID (no foreign key to an object table: a destroyed body's row is retired and its ID never reused), `fate` what
-- happened to it (`destroyed`, `survived`, `absorbed` or `created`), and the rest the body's state at the event: name, class,
-- mass, radius, composition, its sector and galactic position (pc) and velocity (km/s), and whether an admin had edited it.
CREATE TABLE IF NOT EXISTS orbital_event_bodies (
    event_id                    BIGINT UNSIGNED NOT NULL,
    role                        VARCHAR(12) NOT NULL,
    ordinal                     SMALLINT UNSIGNED NOT NULL,
    uid                         BINARY(12) NOT NULL,
    object_kind                 VARCHAR(24) NOT NULL,
    fate                        VARCHAR(12) NULL,
    name                        VARCHAR(255) NULL,
    object_class                VARCHAR(64) NULL,
    mass_kg                     DOUBLE NULL,
    radius_km                   DOUBLE NULL,
    composition                 VARCHAR(64) NULL,
    sector_id                   BIGINT UNSIGNED NULL,
    x_pc                        DOUBLE NULL,
    y_pc                        DOUBLE NULL,
    z_pc                        DOUBLE NULL,
    vx_kms                      DOUBLE NULL,
    vy_kms                      DOUBLE NULL,
    vz_kms                      DOUBLE NULL,
    admin_edited                TINYINT NOT NULL DEFAULT 0,
    note                        VARCHAR(255) NULL,
    PRIMARY KEY (event_id, role, ordinal),
    KEY idx_orbital_event_bodies_uid (uid, event_id),
    CONSTRAINT fk_orbital_event_bodies_event FOREIGN KEY (event_id) REFERENCES orbital_events(id) ON DELETE CASCADE,
    CONSTRAINT chk_orbital_event_bodies_role CHECK (role IN ('participant', 'product')),
    CONSTRAINT chk_orbital_event_bodies_fate CHECK (fate IS NULL OR fate IN ('destroyed', 'survived', 'absorbed', 'created')),
    CONSTRAINT chk_orbital_event_bodies_edited CHECK (admin_edited IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- object_timing_stats (v85, PERF.92): one row per object kind a run built in a
-- sector: `object_count` objects, `total_us` microseconds of their own time (the
-- time spent building objects inside them taken out), the fastest and slowest,
-- and `histogram`, 24 comma-separated counts in doubling-microsecond buckets
-- (`object_timing.HISTOGRAM_BUCKETS`). A run's figures for a kind are the sum of
-- its sector rows. `sector_id` is NULL for work not tied to one sector.
-- object_timing_samples: the objects kept with their full timing, `sampled` = 1
-- for the repeatable 1-in-1,000 sample (a hash of the unit, kind and `ordinal`,
-- the object's number within its kind in that sector) and 0 for the extra rows a
-- benchmark run asks for.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS object_timing_stats (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    run_id        BIGINT UNSIGNED NOT NULL,
    sector_id     BIGINT UNSIGNED NULL,
    kind          VARCHAR(24) NOT NULL,
    object_count  INT UNSIGNED NOT NULL,
    total_us      BIGINT UNSIGNED NOT NULL,
    min_us        BIGINT UNSIGNED NOT NULL,
    max_us        BIGINT UNSIGNED NOT NULL,
    histogram     VARCHAR(255) NOT NULL,
    KEY idx_object_timing_stats_run (run_id, kind),
    KEY idx_object_timing_stats_sector (sector_id),
    CONSTRAINT fk_object_timing_stats_run FOREIGN KEY (run_id) REFERENCES generation_runs(id) ON DELETE CASCADE,
    CONSTRAINT fk_object_timing_stats_sector FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS object_timing_samples (
    id         BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    run_id     BIGINT UNSIGNED NOT NULL,
    sector_id  BIGINT UNSIGNED NULL,
    kind       VARCHAR(24) NOT NULL,
    ordinal    INT UNSIGNED NOT NULL,
    own_us     BIGINT UNSIGNED NOT NULL,
    total_us   BIGINT UNSIGNED NOT NULL,
    sampled    TINYINT(1) NOT NULL,
    KEY idx_object_timing_samples_run (run_id, kind),
    KEY idx_object_timing_samples_sector (sector_id),
    CONSTRAINT fk_object_timing_samples_run FOREIGN KEY (run_id) REFERENCES generation_runs(id) ON DELETE CASCADE,
    CONSTRAINT fk_object_timing_samples_sector FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- phenomenon_scatter_classes (v77): how many scattered phenomena of each
-- class were placed and how many a sector's save has built since, so the
-- Phenomena table counts and pages the scatter without scanning its
-- hundred million rows. One row per (`kind`, `subtype`); `subtype` is ''
-- where the class has none. `placed` is written when the scatter ends
-- (`store.record_phenomenon_scatter_classes`) and `built` goes up as
-- `store.mark_phenomena_built` stamps rows; a plan re-run or reset empties
-- the table with the scatter.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS phenomenon_scatter_classes (
    kind     VARCHAR(24) NOT NULL,
    subtype  VARCHAR(16) NOT NULL DEFAULT '',
    placed   BIGINT UNSIGNED NOT NULL DEFAULT 0,
    built    BIGINT UNSIGNED NOT NULL DEFAULT 0,
    PRIMARY KEY (kind, subtype)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- ---------------------------------------------------------------------
-- id_blocks (v45, PERF.13): the next unreserved id for each table in
-- `_db.ID_BLOCK_TABLES`. `_db._allocate_id` moves `next_id` up by a
-- block at a time (never below the table's own MAX(id) + 1) and the
-- writer then names those ids in its INSERTs. Reserved on a separate
-- autocommitted connection, so writers don't queue behind each other.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS id_blocks (
    table_name  VARCHAR(64) NOT NULL PRIMARY KEY,
    next_id     BIGINT UNSIGNED NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- id_counters (v78, DB.20): the next number of each run-time counter of
-- object IDs. `kind` 'sector': `scope` is a birth sector's address (7
-- bytes) and the number the next run-time serial there. `kind` 'system':
-- `scope` is a system's ID (12 bytes) and the number the next body number
-- in it. `kind` 'sdel' (PERF.38), empty `scope`: the deletion epoch, one more
-- with every sector deleted, so the Galaxy Map's stamp tells a deletion
-- without counting the placed sectors. A row only moves up, on its own autocommitted connection
-- (`store._reserve_counter`), and survives the deletion of the object or
-- sector, so a number is never given twice. Kept by `planetgen reset`
-- like `id_blocks`.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS id_counters (
    kind        VARCHAR(8) NOT NULL,
    scope       VARBINARY(16) NOT NULL,
    next_value  BIGINT UNSIGNED NOT NULL,
    PRIMARY KEY (kind, scope)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- Population and politics (v44, POP.1 to POP.4): see
-- docs/design/habitability/population-and-politics.md. Filled by `planetgen
-- population` (planetgen/population/model.py) from what is already
-- stored; nothing in system generation writes these.
--
-- species: one dominant species per life world (a planet whose
-- evolutionary milestone is multicellularity or a technological
-- civilization). civilization_age_years and era are NULL without a
-- civilization; spacefaring = era Interstellar or later.
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS species (
    id                       BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    name                     VARCHAR(255) NOT NULL,
    homeworld_planet_id      BIGINT UNSIGNED NOT NULL,
    star_system_id           BIGINT UNSIGNED NOT NULL,
    life_class               VARCHAR(2),
    life_stage               VARCHAR(32) NOT NULL,
    build                    VARCHAR(16) NOT NULL,
    climate                  VARCHAR(16) NOT NULL,
    size                     VARCHAR(16) NOT NULL,
    civilization_age_years   DOUBLE,
    era                      VARCHAR(16),
    spacefaring              TINYINT(1) NOT NULL DEFAULT 0 CHECK (spacefaring IN (0, 1)),
    -- v88 (POP.9): the six tech indexes, the level and where it came from, see header comment.
    energy_index             TINYINT UNSIGNED CHECK (energy_index BETWEEN 0 AND 7),
    materials_index          TINYINT UNSIGNED CHECK (materials_index BETWEEN 0 AND 7),
    information_index        TINYINT UNSIGNED CHECK (information_index BETWEEN 0 AND 7),
    medical_index            TINYINT UNSIGNED CHECK (medical_index BETWEEN 0 AND 7),
    propulsion_index         TINYINT UNSIGNED CHECK (propulsion_index BETWEEN 0 AND 7),
    defense_index            TINYINT UNSIGNED CHECK (defense_index BETWEEN 0 AND 7),
    tech_level               DECIMAL(3,2),
    tech_source              VARCHAR(16) CHECK (tech_source IN ('generated', 'admin')),
    generator_epoch          SMALLINT UNSIGNED NULL,  -- v164 (GEN.764): the epoch the population pass made this species under
    created_at               TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    UNIQUE KEY uq_species_name (name),
    UNIQUE KEY uq_species_homeworld (homeworld_planet_id),
    KEY idx_species_system (star_system_id),
    KEY idx_species_spacefaring_name (spacefaring, name),     -- v177 (PERF.205): the Species table's sorts read in index order
    KEY idx_species_tech_level_name (tech_level, name),       -- (the tech level search by range uses the leading column)
    KEY idx_species_era_name (era, name),
    CONSTRAINT chk_species_life_stage CHECK (life_stage IN ('multicellularity', 'technological_civilization')),
    CONSTRAINT fk_species_planet
        FOREIGN KEY (homeworld_planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    CONSTRAINT fk_species_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- polities: one government per spacefaring species, its capital the
-- homeworld's system. reach_ly is how far its territory extends.
CREATE TABLE IF NOT EXISTS polities (
    id                  BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    name                VARCHAR(255) NOT NULL,
    species_id          BIGINT UNSIGNED NOT NULL,
    capital_system_id   BIGINT UNSIGNED NOT NULL,
    government          VARCHAR(32) NOT NULL,
    color               CHAR(7) NOT NULL,
    reach_ly            DOUBLE NOT NULL,
    created_at          TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    UNIQUE KEY uq_polities_name (name),
    UNIQUE KEY uq_polities_species (species_id),
    KEY idx_polities_capital (capital_system_id),
    CONSTRAINT fk_polities_species
        FOREIGN KEY (species_id) REFERENCES species(id) ON DELETE CASCADE,
    CONSTRAINT fk_polities_capital
        FOREIGN KEY (capital_system_id) REFERENCES star_systems(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- facilities (v42): starbases, colonies and outposts. Exactly one host,
-- named by `host_type`: a star (`star_system_id`, plus `star_id` unless
-- it's a close pair's shared orbit), a planet, a moon, an asteroid belt,
-- an asteroid field, or open space (`sector_id` plus a galaxy-frame
-- center, like a phenomenon). `star_system_id` is set for every host in a
-- system so a system's facilities list in one lookup. `placement` is
-- terrestrial, orbital, asteroid or standalone, and
-- `program_constants.FACILITY_RULES` says which kinds may go where
-- (`facilities.check_facility`). An orbital facility's circular orbit
-- (`orbit_distance_km`, `orbit_period_years`, `orbital_speed_kms`,
-- `orbit_phase_deg`) comes from its host's mass the way a moon's does.
-- Every host foreign key cascades, so a deleted host takes its
-- facilities along. "Exactly one host" is checked by `_db.add_facility`,
-- not a CHECK: MySQL refuses a CHECK on a column with a cascading
-- foreign key.
-- ---------------------------------------------------------------------
-- facility_types (v89, DB.28): the admin-defined refinements of a facility kind, see header comment.
CREATE TABLE IF NOT EXISTS facility_types (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    name                 VARCHAR(255) NOT NULL,
    kind                 VARCHAR(16) NOT NULL,
    icon                 VARCHAR(64),
    description          TEXT,
    sort_order           INT NOT NULL DEFAULT 0,
    active               TINYINT(1) NOT NULL DEFAULT 1,

    UNIQUE KEY uq_facility_types_name (name),
    KEY idx_facility_types_kind (kind),
    CONSTRAINT chk_facility_types_kind CHECK (kind IN ('colony', 'outpost', 'mining-colony', 'station', 'starbase')),
    CONSTRAINT chk_facility_types_active CHECK (active IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- facility_type_fields (v89, DB.28): the typed fields a facility type carries.
CREATE TABLE IF NOT EXISTS facility_type_fields (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    type_id              BIGINT UNSIGNED NOT NULL,
    `key`                VARCHAR(32) NOT NULL,
    label                VARCHAR(255) NOT NULL,
    datatype             VARCHAR(8) NOT NULL,
    unit                 VARCHAR(32),
    required             TINYINT(1) NOT NULL DEFAULT 0,
    min_value            DOUBLE,
    max_value            DOUBLE,
    choices              TEXT,
    default_text         VARCHAR(255),
    sort_order           INT NOT NULL DEFAULT 0,

    UNIQUE KEY uq_facility_type_fields_key (type_id, `key`),
    CONSTRAINT chk_facility_type_fields_datatype CHECK (datatype IN ('int', 'decimal', 'text', 'bool', 'choice')),
    CONSTRAINT chk_facility_type_fields_required CHECK (required IN (0, 1)),
    CONSTRAINT fk_facility_type_fields_type
        FOREIGN KEY (type_id) REFERENCES facility_types(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS facilities (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    name                 VARCHAR(255) NOT NULL,
    kind                 VARCHAR(16) NOT NULL,
    placement            VARCHAR(16) NOT NULL,
    host_type            VARCHAR(16) NOT NULL,
    star_system_id       BIGINT UNSIGNED,
    star_id              BIGINT UNSIGNED,
    planet_id            BIGINT UNSIGNED,
    moon_id              BIGINT UNSIGNED,
    asteroid_belt_id     BIGINT UNSIGNED,
    asteroid_field_id    BIGINT UNSIGNED,
    sector_id            BIGINT UNSIGNED,
    center_x_pc          DOUBLE,
    center_y_pc          DOUBLE,
    center_z_pc          DOUBLE,
    -- v93 (region queries): where the object is now as one integer, the Morton key of its 0.25 pc cube; MySQL computes it,
    -- NULL outside the box of plus or minus 262,144 pc. `planetgen/galaxy/cellkey.py` has the same function.
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    mass_kg                     DOUBLE NULL,  -- v102 (DB.55): point mass, a stated zero where the mass is not known
    KEY idx_facilities_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    galactic_radius_pc   DOUBLE,
    velocity_x_kms       DOUBLE NOT NULL DEFAULT 0,
    velocity_y_kms       DOUBLE NOT NULL DEFAULT 0,
    velocity_z_kms       DOUBLE NOT NULL DEFAULT 0,
    orbit_distance_km    DOUBLE,
    orbit_period_years   DOUBLE,
    orbital_speed_kms    DOUBLE,
    orbit_phase_deg      DOUBLE,
    epoch_unix           DOUBLE,  -- v66 (GEN.106), see header comment
    next_update_due      DOUBLE,
    description          TEXT,
    uid                  BINARY(12),  -- v78 (DB.20), see header comment
    -- v89 (DB.28, POP.10): type, affiliation and ratings, see header comment.
    type_id              BIGINT UNSIGNED,
    affiliation_polity_id BIGINT UNSIGNED,
    affiliation_name     VARCHAR(255),
    rating_crime         TINYINT UNSIGNED,
    rating_housing       TINYINT UNSIGNED,
    rating_resources     TINYINT UNSIGNED,
    rating_maintenance   TINYINT UNSIGNED,
    rating_health        TINYINT UNSIGNED,
    ratings_source       VARCHAR(16),
    ratings_note         TEXT,
    ratings_set_at       TIMESTAMP NULL,

    created_at           TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    modified_at          TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),

    KEY idx_facilities_type (type_id),
    KEY idx_facilities_affiliation (affiliation_polity_id),
    KEY idx_facilities_rating_crime (rating_crime),
    KEY idx_facilities_rating_housing (rating_housing),
    KEY idx_facilities_rating_resources (rating_resources),
    KEY idx_facilities_rating_maintenance (rating_maintenance),
    KEY idx_facilities_rating_health (rating_health),
    KEY idx_facilities_star_system (star_system_id),
    KEY idx_facilities_next_update_due (next_update_due),
    KEY idx_facilities_sector (sector_id),
    KEY idx_facilities_name (name),
    UNIQUE KEY uq_facilities_uid (uid),
    CONSTRAINT chk_facilities_kind CHECK (kind IN ('colony', 'outpost', 'mining-colony', 'station', 'starbase')),
    CONSTRAINT chk_facilities_placement CHECK (placement IN ('terrestrial', 'orbital', 'asteroid', 'standalone')),
    CONSTRAINT chk_facilities_rating_crime CHECK (rating_crime BETWEEN 1 AND 3),
    CONSTRAINT chk_facilities_rating_housing CHECK (rating_housing BETWEEN 1 AND 3),
    CONSTRAINT chk_facilities_rating_resources CHECK (rating_resources BETWEEN 1 AND 3),
    CONSTRAINT chk_facilities_rating_maintenance CHECK (rating_maintenance BETWEEN 1 AND 3),
    CONSTRAINT chk_facilities_rating_health CHECK (rating_health BETWEEN 1 AND 3),
    CONSTRAINT fk_facilities_star_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_facilities_star
        FOREIGN KEY (star_id) REFERENCES stars(id) ON DELETE CASCADE,
    CONSTRAINT fk_facilities_planet
        FOREIGN KEY (planet_id) REFERENCES planets(id) ON DELETE CASCADE,
    CONSTRAINT fk_facilities_moon
        FOREIGN KEY (moon_id) REFERENCES moons(id) ON DELETE CASCADE,
    CONSTRAINT fk_facilities_asteroid_belt
        FOREIGN KEY (asteroid_belt_id) REFERENCES asteroid_belts(id) ON DELETE CASCADE,
    CONSTRAINT fk_facilities_asteroid_field
        FOREIGN KEY (asteroid_field_id) REFERENCES asteroid_fields(id) ON DELETE CASCADE,
    CONSTRAINT fk_facilities_sector
        FOREIGN KEY (sector_id) REFERENCES sectors(id) ON DELETE CASCADE,
    CONSTRAINT fk_facilities_type
        FOREIGN KEY (type_id) REFERENCES facility_types(id) ON DELETE RESTRICT,
    CONSTRAINT fk_facilities_affiliation_polity
        FOREIGN KEY (affiliation_polity_id) REFERENCES polities(id) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- facility_field_values (v89, DB.28): one typed value per facility and field.
CREATE TABLE IF NOT EXISTS facility_field_values (
    facility_id          BIGINT UNSIGNED NOT NULL,
    field_id             BIGINT UNSIGNED NOT NULL,
    value_int            BIGINT,
    value_dec            DOUBLE,
    value_text           VARCHAR(255),
    value_bool           TINYINT(1),

    PRIMARY KEY (facility_id, field_id),
    KEY idx_facility_field_values_field (field_id),
    CONSTRAINT fk_facility_field_values_facility
        FOREIGN KEY (facility_id) REFERENCES facilities(id) ON DELETE CASCADE,
    CONSTRAINT fk_facility_field_values_field
        FOREIGN KEY (field_id) REFERENCES facility_type_fields(id) ON DELETE RESTRICT
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- system_owners: the polity with the strongest claim (reach / distance)
-- on each system inside some polity's reach. Recomputed from scratch by
-- every population pass.
CREATE TABLE IF NOT EXISTS system_owners (
    star_system_id   BIGINT UNSIGNED PRIMARY KEY,
    polity_id        BIGINT UNSIGNED NOT NULL,
    distance_ly      DOUBLE NOT NULL,

    KEY idx_system_owners_polity (polity_id),
    CONSTRAINT fk_system_owners_system
        FOREIGN KEY (star_system_id) REFERENCES star_systems(id) ON DELETE CASCADE,
    CONSTRAINT fk_system_owners_polity
        FOREIGN KEY (polity_id) REFERENCES polities(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- population_state: one row (id = 1), the highest planets.id the
-- population pass has scanned.
CREATE TABLE IF NOT EXISTS population_state (
    id                     TINYINT UNSIGNED PRIMARY KEY,
    scanned_planet_id      BIGINT UNSIGNED NOT NULL DEFAULT 0,
    CHECK (id = 1)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- galaxy_components (v97, GEN.159): the discrete stellar components of a galaxy beyond the disk, one row each, of every
-- kind in one table: the globular clusters, and later the satellites and streams (GEN.237). A component is a density
-- (a profile of its mass about a centre), not a list of stars; its stars come from the profile when a sector is
-- generated. `profile` 'king' carries the King (1966) concentration `king_w0` (3 to 12) and the `tidal_radius_pc` (`Component` checks that, not the table);
-- `planetgen/physics/king.py` turns a row into a density. `source` is 'catalogue' for a real object and 'generated' for
-- one drawn from the galaxy seed (`planetgen/generation/clusters.py`); `planetgen plan` draws the clusters (GEN.160).
CREATE TABLE IF NOT EXISTS galaxy_components (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    kind                 VARCHAR(12) NOT NULL,
    name                 VARCHAR(255) NULL,  -- v127 (DB.61): NULL for a cluster nobody has named, the unique key is over the named rows
    -- The catalogue's own name for the object ("NGC 104"), so a later import of a real catalogue can match its rows.
    catalogue_id         VARCHAR(40),
    source               VARCHAR(12) NOT NULL,
    -- A cluster inside a satellite, or a satellite of a larger galaxy; NULL for a component of the galaxy itself.
    parent_component_id  BIGINT UNSIGNED,
    -- The component's centre in the galaxy frame (parsecs) and its velocity there, at `epoch_unix`.
    center_x_pc          DOUBLE NOT NULL,
    center_y_pc          DOUBLE NOT NULL,
    center_z_pc          DOUBLE NOT NULL,
    cell_key             BIGINT UNSIGNED GENERATED ALWAYS AS (IF(ABS(center_x_pc) < 262144 AND ABS(center_y_pc) < 262144 AND ABS(center_z_pc) < 262144, ((((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_x_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 2) | ((((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_y_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)) << 1) | (((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 0) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 0) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 7) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 21) | ((((((((((LEAST(CAST(FLOOR(((center_z_pc) + 262144) / 0.25) AS UNSIGNED), 2097151)) >> 14) & 127) * 257) & 28687) * 17) & 274627) * 5) & 299593) << 42)), NULL)) STORED,
    KEY idx_galaxy_components_cell (cell_key, center_x_pc, center_y_pc, center_z_pc),
    velocity_x_kms       DOUBLE NOT NULL,
    velocity_y_kms       DOUBLE NOT NULL,
    velocity_z_kms       DOUBLE NOT NULL,
    epoch_unix           DOUBLE NOT NULL,
    mass_solar           DOUBLE NOT NULL,
    half_mass_radius_pc  DOUBLE NOT NULL,
    profile              VARCHAR(12) NOT NULL,
    king_w0              DOUBLE,
    tidal_radius_pc      DOUBLE,
    age_gy               DOUBLE NOT NULL,
    metallicity_fe_h     DOUBLE NOT NULL,
    core_collapsed       TINYINT(1) NOT NULL DEFAULT 0,
    seed                 BIGINT UNSIGNED NOT NULL,
    created_at           TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    UNIQUE KEY uq_galaxy_components_name (kind, name),
    KEY idx_galaxy_components_parent (parent_component_id),
    CONSTRAINT fk_galaxy_components_parent
        FOREIGN KEY (parent_component_id) REFERENCES galaxy_components(id) ON DELETE SET NULL,
    CONSTRAINT chk_galaxy_components_kind CHECK (kind IN ('cluster', 'satellite', 'stream')),
    CONSTRAINT chk_galaxy_components_source CHECK (source IN ('catalogue', 'generated')),
    CONSTRAINT chk_galaxy_components_profile CHECK (profile IN ('king', 'exponential', 'plummer')),
    CONSTRAINT chk_galaxy_components_core_collapsed CHECK (core_collapsed IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- v113 (DB.79): the object index. One row per object of every kind the site lists; see the header note.
CREATE TABLE IF NOT EXISTS object_index (
    kind            VARCHAR(24) NOT NULL,
    uid             BINARY(12) NOT NULL,
    source_id       BIGINT UNSIGNED NOT NULL,
    name            VARCHAR(255) NOT NULL,
    sector_id       BIGINT UNSIGNED NULL,
    system_id       BIGINT UNSIGNED NULL,
    parent_uid      BINARY(12) NULL,
    object_class    VARCHAR(64) NULL,
    mass_kg         DOUBLE NULL,
    radius_km       DOUBLE NULL,
    placed          TINYINT(1) NOT NULL DEFAULT 1,
    indexed_at      TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    PRIMARY KEY (kind, uid),
    KEY idx_object_index_source (kind, source_id),
    KEY idx_object_index_name (kind, name, uid),
    KEY idx_object_index_class (kind, object_class, name, uid),
    KEY idx_object_index_mass (kind, mass_kg, uid),
    KEY idx_object_index_radius (kind, radius_km, uid),
    KEY idx_object_index_sector (sector_id, kind, name, uid),
    KEY idx_object_index_system (system_id, kind, uid),
    FULLTEXT KEY ft_object_index_name (name),
    CONSTRAINT chk_object_index_placed CHECK (placed IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v130 (DB.85): the stored globular clusters and superbubbles; see the header note.
CREATE TABLE IF NOT EXISTS galaxy_containers (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    uid           BINARY(12) NOT NULL,
    kind          VARCHAR(16) NOT NULL,
    name          VARCHAR(255) NULL,
    component_id  BIGINT UNSIGNED NULL,
    bubble_col      INT NULL,
    bubble_row      INT NULL,
    bubble_n      SMALLINT UNSIGNED NULL,
    center_x_pc   DOUBLE NOT NULL,
    center_y_pc   DOUBLE NOT NULL,
    center_z_pc   DOUBLE NOT NULL,
    extent_pc     DOUBLE NOT NULL,
    shell_index   INT NOT NULL,
    band_index    INT NOT NULL,
    wedge_index   INT NOT NULL,
    created_at    TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    UNIQUE KEY uq_galaxy_containers_uid (uid),
    UNIQUE KEY uq_galaxy_containers_bubble (bubble_col, bubble_row, bubble_n),
    KEY idx_galaxy_containers_kind_x (kind, center_x_pc),
    KEY idx_galaxy_containers_birth (shell_index, band_index, wedge_index),
    KEY idx_galaxy_containers_component (component_id),
    CONSTRAINT fk_galaxy_containers_component FOREIGN KEY (component_id) REFERENCES galaxy_components(id) ON DELETE CASCADE,
    CONSTRAINT chk_galaxy_containers_kind CHECK (kind IN ('globular_cluster', 'superbubble')),
    CONSTRAINT chk_galaxy_containers_extent CHECK (extent_pc >= 1)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v129 (DB.63): the ships and stations; see the header note.
CREATE TABLE IF NOT EXISTS vessels (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    uid                  BINARY(12) NOT NULL,
    kind                 VARCHAR(8) NOT NULL,
    name                 VARCHAR(255) NOT NULL,
    origin               VARCHAR(12) NOT NULL,
    owner_id             BIGINT UNSIGNED NULL,
    owner_name           VARCHAR(255) NULL,
    center_x_pc          DOUBLE NULL,
    center_y_pc          DOUBLE NULL,
    center_z_pc          DOUBLE NULL,
    mass_kg              DOUBLE NULL,
    velocity_x_kms       DOUBLE NOT NULL DEFAULT 0,
    velocity_y_kms       DOUBLE NOT NULL DEFAULT 0,
    velocity_z_kms       DOUBLE NOT NULL DEFAULT 0,
    epoch_unix           DOUBLE NULL,
    course               JSON NULL,
    seed                 BIGINT UNSIGNED NULL,
    version_key          CHAR(22) NULL,
    created_at           TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    UNIQUE KEY uq_vessels_uid (uid),
    KEY idx_vessels_owner (owner_id, id),
    CONSTRAINT chk_vessels_kind CHECK (kind IN ('ship', 'station')),
    CONSTRAINT chk_vessels_origin CHECK (origin IN ('generated', 'user'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- v128 (DB.64): the stored constellations; see the header note.
CREATE TABLE IF NOT EXISTS constellations (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    system_uid           BINARY(12) NOT NULL,
    constellation_index  SMALLINT UNSIGNED NOT NULL,
    name_id              BINARY(12) NULL,
    members              JSON NOT NULL,
    edges                JSON NOT NULL,
    seed                 BIGINT UNSIGNED NOT NULL,
    version_key          CHAR(22) NOT NULL,
    created_at           TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    UNIQUE KEY uq_constellations_system_index (system_uid, constellation_index),
    KEY idx_constellations_name (name_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- v126 (DB.39): the edit request comments and notices; see the header note.
CREATE TABLE IF NOT EXISTS request_comments (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    plan_id       BIGINT UNSIGNED NOT NULL,
    author_id     BIGINT UNSIGNED NULL,
    author_name   VARCHAR(255) NOT NULL,
    author_role   VARCHAR(16) NOT NULL,
    body          TEXT NOT NULL,
    created_at    TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    KEY idx_request_comments_plan (plan_id, id),
    CONSTRAINT fk_request_comments_plan FOREIGN KEY (plan_id) REFERENCES change_plans(id) ON DELETE CASCADE,
    CONSTRAINT chk_request_comments_author_role CHECK (author_role IN ('user', 'moderator', 'admin'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS request_notifications (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    plan_id       BIGINT UNSIGNED NOT NULL,
    recipient_id  BIGINT UNSIGNED NOT NULL,
    kind          VARCHAR(32) NOT NULL,
    created_at    TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    sent_at       TIMESTAMP(3) NULL,
    KEY idx_request_notifications_recipient (recipient_id, sent_at, id),
    KEY idx_request_notifications_plan (plan_id, id),
    CONSTRAINT fk_request_notifications_plan FOREIGN KEY (plan_id) REFERENCES change_plans(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- v124 (DB.38): the region density overrides; see the header note.
CREATE TABLE IF NOT EXISTS region_density_overrides (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    shape         VARCHAR(16) NOT NULL,
    params        JSON NOT NULL,
    multiplier    DECIMAL(10,4) NOT NULL,
    set_by_id     BIGINT UNSIGNED NULL,
    set_by_name   VARCHAR(255) NOT NULL,
    set_at        TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    reason        VARCHAR(500) NULL,
    undone_at     TIMESTAMP(3) NULL,
    KEY idx_region_density_overrides_active (undone_at, id),
    CONSTRAINT chk_region_density_overrides_shape CHECK (shape IN ('sphere', 'box', 'shell', 'band', 'wedge', 'sectors')),
    CONSTRAINT chk_region_density_overrides_multiplier CHECK (multiplier >= 0.0001)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- v123 (DB.37): the lock marks; see the header note.
CREATE TABLE IF NOT EXISTS object_locks (
    object_uid       BINARY(12) NOT NULL,
    aspect           VARCHAR(64) NOT NULL DEFAULT '*',
    locked_by_id     BIGINT UNSIGNED NULL,
    locked_by_name   VARCHAR(255) NOT NULL,
    locked_at        TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    reason           VARCHAR(500) NULL,
    PRIMARY KEY (object_uid, aspect),
    KEY idx_object_locks_locked_at (locked_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- v122 (DB.36): change plans, their operations and the edit journal; see the header note.
CREATE TABLE IF NOT EXISTS change_plans (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    requester_id     BIGINT UNSIGNED NULL,
    requester_name   VARCHAR(255) NOT NULL,
    requester_role   VARCHAR(16) NOT NULL,
    request_text     TEXT NULL,
    state            VARCHAR(16) NOT NULL DEFAULT 'draft',
    approver_role    VARCHAR(16) NULL,
    decided_by_id    BIGINT UNSIGNED NULL,
    decided_by_name  VARCHAR(255) NULL,
    decided_at       TIMESTAMP(3) NULL,
    decision_note    TEXT NULL,
    apply_job_id     BIGINT UNSIGNED NULL,
    created_at       TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    updated_at       TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    KEY idx_change_plans_state (state, id),
    KEY idx_change_plans_requester (requester_id, id),
    CONSTRAINT chk_change_plans_state CHECK (state IN
        ('draft', 'requested', 'proposed', 'approved', 'applied', 'undone', 'stale', 'rejected', 'withdrawn')),
    CONSTRAINT chk_change_plans_requester_role CHECK (requester_role IN ('user', 'moderator', 'admin')),
    CONSTRAINT chk_change_plans_approver_role CHECK (approver_role IS NULL OR approver_role IN ('moderator', 'admin'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS change_plan_ops (
    id          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    plan_id     BIGINT UNSIGNED NOT NULL,
    op_order    INT NOT NULL,
    kind        VARCHAR(48) NOT NULL,
    target_uid  BINARY(12) NULL,
    before_json JSON NULL,
    after_json  JSON NULL,
    status      VARCHAR(16) NOT NULL DEFAULT 'pending',
    UNIQUE KEY uq_change_plan_ops_order (plan_id, op_order),
    KEY idx_change_plan_ops_target (target_uid),
    CONSTRAINT fk_change_plan_ops_plan FOREIGN KEY (plan_id) REFERENCES change_plans(id) ON DELETE CASCADE,
    CONSTRAINT chk_change_plan_ops_status CHECK (status IN ('pending', 'applied', 'undone', 'skipped', 'failed'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS edit_journal (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    plan_id          BIGINT UNSIGNED NOT NULL,
    op_id            BIGINT UNSIGNED NULL,
    object_uid       BINARY(12) NOT NULL,
    sector_id        BIGINT UNSIGNED NULL,
    aspect           VARCHAR(64) NOT NULL,
    before_json      JSON NULL,
    after_json       JSON NULL,
    generator_epoch  SMALLINT UNSIGNED NOT NULL,
    applied_at       TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    KEY idx_edit_journal_object (object_uid, id),
    KEY idx_edit_journal_sector (sector_id, id),
    KEY idx_edit_journal_plan (plan_id, id),
    CONSTRAINT fk_edit_journal_plan FOREIGN KEY (plan_id) REFERENCES change_plans(id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- v120 (DB.49): the digests above the sector leaves; see the header note.
CREATE TABLE IF NOT EXISTS sector_tree_nodes (
    level         TINYINT UNSIGNED NOT NULL,
    shell_index   INT NOT NULL,
    band_index    INT NOT NULL DEFAULT 0,
    digest        BINARY(32) NULL,
    child_count   INT UNSIGNED NOT NULL DEFAULT 0,
    dirty         TINYINT(1) NOT NULL DEFAULT 1,
    updated_at    TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),
    PRIMARY KEY (level, shell_index, band_index),
    KEY idx_sector_tree_nodes_dirty (dirty, level),
    CONSTRAINT chk_sector_tree_nodes_level CHECK (level IN (1, 2)),
    CONSTRAINT chk_sector_tree_nodes_shell_band CHECK (level = 1 OR band_index = 0),
    CONSTRAINT chk_sector_tree_nodes_dirty CHECK (dirty IN (0, 1))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;


-- v119 (DB.30): the per-species environment limits; see the header note.
CREATE TABLE IF NOT EXISTS species_limit_values (
    species_id   BIGINT UNSIGNED NOT NULL,
    limit_name   VARCHAR(48) NOT NULL,
    value        DOUBLE NOT NULL,
    PRIMARY KEY (species_id, limit_name),
    CONSTRAINT fk_species_limit_values_species
        FOREIGN KEY (species_id) REFERENCES species(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ---------------------------------------------------------------------
-- sector_objects -- unified "every stellar object in a sector" search.
--
-- Standard relational choice here: a VIEW that UNIONs the five typed
-- tables above, not a sixth physical table duplicating their rows. A
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
    JOIN star_systems ss ON ss.id = ab.star_system_id

    UNION ALL

    SELECT
        'comet'                         AS object_type,
        c.id                            AS object_id,
        c.star_system_id                AS star_system_id,
        ss.sector_id                    AS sector_id,
        c.name                          AS name,
        c.orbit_type                    AS summary,
        NULL                            AS orbital_index
    FROM comets c
    JOIN star_systems ss ON ss.id = c.star_system_id;

SET FOREIGN_KEY_CHECKS = 1;
