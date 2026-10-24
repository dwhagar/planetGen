# Numeric precision audit at galactic scale (GEN.372)

Boss (6:54 PM and 6:55 PM EDT, 2026-10-16): are 64-bit numbers enough for positions at galactic scale, or do we need 128 or 256 bits? This page is the one place the answer, the numbers behind it and every precision concern are collected. It audits; it builds nothing. Owner: Core 1 (point columns, IDs). The camera and render rows were read by Core 2. A lane that meets a precision bug or concern sends it to the TODO thread, which attaches it to GEN.372.

## Short answer

**64-bit floats are enough, provided no code takes the difference of two large absolute coordinates where the small difference matters.** A 64-bit float has a relative precision of 2.2e-16. Out to the 262 kpc reach that is 28 km at the Sun's distance and 898 km at the edge, against the smallest move the galaxy stores (0.01 milliparsec, 3.1e8 km) a factor of 3.4e5 (at the edge) to 1.1e7 (at the Sun) larger. The design that holds this together already exists: a coarse integer address (the sector, 96-bit object IDs, the Morton cell key) plus a local offset, with systems and moons in their own frames. 128-bit or 256-bit floats are not needed for storage. Where a wider sum is justified is narrow (section 4).

The audit found nine places where the code or its documents do not match that design (section 2) and eight cases that need work (section 5, one of them, the documents, fixed with this page), and no case where a float wider than 64 bits is the right fix.

## 1. The numbers

The spacing of float64 is `2^(e-52)` for a value in `[2^e, 2^(e+1))`. A value of 262,144 pc (the reach, exactly 2^18) is the first of its band, so the spacing just below it is 2.9e-11 pc and just above it 5.8e-11 pc (1,796 km). The same distance in other units differs by up to a factor of two, which is why older figures disagree.

| Distance from origin | float64 in parsecs | in km | float64 in metres | float32 in parsecs |
|---|---|---|---|---|
| 1 pc | 2.2e-16 | 0.007 | 4 m | 1.2e-7 pc |
| 100 pc | 1.4e-14 | 0.44 | 512 m | 7.6e-6 pc |
| 1 kpc | 1.1e-13 | 3.5 | 4.1 km | 6.1e-5 pc |
| 8.128 kpc (the Sun) | 9.1e-13 | 28.1 | 32.8 km | 4.9e-4 pc (about 100 AU) |
| 15 kpc | 1.8e-12 | 56.1 | 65.5 km | 9.8e-4 pc |
| 30 kpc | 3.6e-12 | 112 | 131 km | 2.0e-3 pc |
| 65.5 kpc | 1.5e-11 | 449 | 262 km | 7.8e-3 pc |
| 131 kpc to just under 262,144 pc | 2.9e-11 | 898 | 1,049 km at the edge | 3.1e-2 pc (about 6,400 AU) |

Stored coordinates are parsecs or milliparsecs, not metres, so the figure for the Sun's position is 28.1 km (parsecs) or 28.7 km (milliparsecs); the 32.8 km quoted in `orbital-updates.md` is the spacing of a coordinate in metres. The 14 km quoted for GEN.142 is half that spacing (the largest error of a rounded value), and `collisions-and-mergers.md` quotes 3.3e4 m for metres and 28 km for parsec doubles; both are right for their unit.

What has to survive:

| Smallest change that must survive | Size | Spacings at 8.128 kpc | Spacings at 262 kpc |
|---|---|---|---|
| The galactic move threshold of a star (`THRESHOLDS_M["galactic"]`, 0.01 mpc) | 3.1e8 km | 1.1e7 | 3.4e5 |
| A planet's orbit stored in the system frame | 1.5e8 km, against 0.002 km spacing at 1.5e13 km | 7.5e10 | 7.5e10 |
| A planet added to a star's galactic position | 1.5e8 km | 5.3e6 | 1.7e5 |
| A moon's move threshold (`THRESHOLDS_M["planetary"]`, 100,000 km) added to a galactic position | 1e5 km | 3.6e3 | 111 |
| A scatter member's offset of 1 AU added to a galactic position | 1.5e8 km | 5.3e6 | 1.7e5 |

A moon's threshold and a planet's offset are why systems and moons keep their own frames, and why the points below matter.

## 2. Where the code and the documents disagree

| # | Statement | What is true | Evidence |
|---|---|---|---|
| 1 | "Positions are stored in integer milliparsecs relative to the sector" (`orbital-updates.md` sections 7 and 9) | Star positions in a sector are DOUBLE milliparsecs, not integers. Scatter points (`bright_stars`, `phenomenon_scatter`) are DOUBLE milliparsecs in the **galaxy** frame; they are whole milliparsecs only at birth and arbitrary doubles after the first orbit update. | `db/schema.sql` `star_systems.position_x_mpc`, `bright_stars.position_x_mpc`; `generation/bright_stars.py` `round(value * 1000)` |
| 2 | "32.8 km at 8.128 kpc" | That is metres. Stored parsec coordinates are 28.1 km and milliparsec coordinates 28.7 km (section 1). | numpy spacing |
| 3 | The sector grid is an integer cube grid | It is spherical shells (4 pc deep), bands and wedges: three integers (`galaxy/spherical.py`), with a local offset along the sector's own axes (`galaxy/frame.py`). | `docs/design/galaxy/spherical-shell-sectors.md` |
| 4 | The map tile root covers the galaxy | The tile root is plus or minus 32,768 pc (65,536 pc edge); `tile_keys_containing` returns nothing outside it, so the 262 kpc reach is not tiled. | `galaxy/viewport.py` 205-209, `galaxymap3d.js` `TILE_ROOT` |
| 5 | `point_columns.FRAMES`: a system-frame position is "from the object's star" | A moon is measured from its parent planet. Belt velocity is in galactic axes while planet, moon and comet velocity is relative to the anchor. | `db/point_columns.py`, `db/schema.sql` |
| 6 | Fingerprints quantise to 9 significant digits (`generation-determinism.md`) | The code does not quantise: it hashes the shortest round-trip text and folds -0.0 to 0.0. | `db/fingerprint.py` `canonical_value` |
| 7 | One phase formula | The SQL phase updates and `_galactic_turn` reduce first (`MOD(elapsed, period) / period * 360`); `physics/body_positions.py`, `system_tree.node_relative_au` and the browser's `orbitpositions.js` add `360 * years / period` to the phase and reduce after. For a stored phase and an elapsed time of many periods they differ in the last digits that the reduce-first form keeps. | `db/store.py`, `physics/body_positions.py` |
| 8 | "A printed object ID is 22" | A 96-bit ID prints as 24 hex digits, 26 characters with hyphens. The JavaScript comment `shell << 28 | (band+2048) << 16` is stale: the code packs 19, 18 and 19 bits with a bias of 131072. | `web/lib/tilepack.py` 30; `galaxyprisms.js` 475 |
| 9 | Object-ID resolution of 0.0008 pc at 500 pc and 0.022 pc at 15 kpc (`object-ids.md`) | That describes the rejected positional ID. The 96-bit ID is a birth identity (sector address, serial, body), with no positional resolution. | `docs/design/architecture/object-id-options.md` |

## 3. Every stored or computed float quantity

Verdict key: **OK** float64 is enough as the code stands. **WATCH** enough, but one rule must hold (stated). **WORK** needs a change (item in section 5). **PIN** a test fixes the number.

### 3.1 Positions

| Quantity | Stored as | Origin | Largest distance | float64 spacing there | Must survive | First place error shows | Verdict | Test |
|---|---|---|---|---|---|---|---|---|
| Sector centre `sectors.center_*_pc` | DOUBLE pc | galaxy | 262,144 pc | 898 km | derived from the integer address, never stored as truth | none | OK | `test_sector_centre_round_trips_the_address` at the reach |
| Sector address `shell, band, wedge`; sector `uid` (56 bits) | integers | none | shell 65,535 at 262 kpc (limit 83,037) | exact | all | none | OK | `check_reach` boundary test |
| Star system in its sector `position_*_mpc` | DOUBLE mpc | sector centre, local axes | about 3.7 pc half-diagonal | 1e-13 mpc | whole stored move of 0.01 mpc | none | OK | covered by existing sector tests |
| Star system in the galaxy `galaxy_*_pc` | DOUBLE pc | galactic centre | 262,144 pc | 898 km | 0.01 mpc move (3.1e8 km) | the nightly update rotates this absolute point and recomputes the local offset as `absolute - centre` (both about 8 kpc, difference about a parsec), so the local offset carries the absolute coordinate's error (28 km at the Sun), and the two columns can disagree by an ulp of the absolute coordinate | WATCH: the move is 1e7 spacings, so the loss is invisible, but the local offset must not be treated as exact to its own precision | `test_local_offset_after_a_rotation_is_within_one_ulp_of_the_absolute_point` at 8 kpc and at 262 kpc |
| Free objects (black holes, neutron stars, nebulae and remnants, quasar, rogue planets, interstellar comets, asteroid fields, facilities) `center_*_pc` | DOUBLE pc | galactic centre | 262,144 pc | 898 km | the 0.01 mpc move | none | OK | `test_free_object_move_survives_at_the_reach` |
| `bright_stars`, `phenomenon_scatter` `position_*_mpc` | DOUBLE mpc | galactic centre | 262,144,000 mpc | 3e-8 mpc (920 km) | the 0.01 mpc move; conversion `mpc / 1000` before the key | none for the move; see the cell key (3.6) | OK | `test_scatter_point_survives_a_small_move_at_the_reach` |
| Scatter members `offset_*_au` | DOUBLE AU | the scatter centre | the system extent | far below 1 AU | an offset of AU size added to the centre in pc | `member_positions_pc` adds `offset / AU_PER_PARSEC` to the absolute centre: the offset is rounded to the centre's spacing (28 km at the Sun) before it is used as a galactic position. A member 0.1 AU (1.5e7 km) from its partner is 5e5 spacings, so fine; a contact binary at 0.01 AU is 5e4 spacings, fine for display, not for dynamics | WATCH: members are only ever placed relative to each other (`member_offsets_au`); never difference two absolute member positions | `test_member_offsets_difference_in_the_local_frame` |
| Stars inside a system `stars.offset_*_au` (DB.53) | DOUBLE AU | system centre | system extent (about 1e5 AU for a Sun at 8 kpc, computed) | 1.5e-11 AU at 1e5 AU | an orbit of 0.01 AU | none | OK | `test_star_offset_survives_at_the_system_edge` |
| Planet, moon, comet, belt `position_*_km` | DOUBLE km | parent star or planet | long-period comet aphelion about 7.5e12 km (computed from `tuning.COMET_PERIOD_CLASSES`) | 0.002 km at 1.5e13 km | a moon's 100,000 km move, a planet's orbit | none | OK | `test_comet_aphelion_position_resolves_a_kilometre` |
| Position object `SpatialPosition3D` | metres (or ly/AU by unit) | any frame; the truth frame is kept as set | 8.09e21 m at 262 kpc | 1,049 km | frame conversions | converting a moon's metres to galactic metres loses its digits below 1,049 km at the reach; the class keeps the truth frame so that is not done to stored data | OK | `test_position_keeps_the_moon_in_its_own_frame` |
| Sector paths `sector_path_knots` | DOUBLE pc, km/s | galaxy | 262,144 pc | 898 km | neighbour pulls | built as `(centre_pc + offset / 1000) * PARSEC_M` in metres; the kernel then differences absolute metres (3.4) | WATCH | see 3.4 |
| Hypervelocity straight-line motion `p + v * elapsed / KM_PER_PC` | Python float | galaxy | 262,144 pc | 898 km | elapsed times far below a Gyr | none | OK | `test_straight_line_motion_over_a_gigayear` |

### 3.2 Velocities, epochs, phases and elements

| Quantity | Stored as | Spacing | Must survive | First place error shows | Verdict | Test |
|---|---|---|---|---|---|---|
| Velocities `velocity_*_kms` | DOUBLE | 1e-15 km/s relative | 1 km/s over a Gyr is about 1,000 pc, so a relative 1e-16 gives 1e-13 pc, far below the position spacing | none | OK | covered by the motion tests |
| `epoch_unix` (about 1.8e9 s) | DOUBLE | 2.4e-7 s | a nightly step | none | OK | none |
| `next_update_due` near the cap 3.2e16 s | DOUBLE | 4 s | scheduling only | none | OK | none |
| The orbit clock `orbit_simulation_state.last_updated_at` | TIMESTAMP, whole seconds | 1 s | `finish_orbit_update` floors it so the next run never reads a clock that went backwards | none | OK | `test_orbit_clock_never_runs_backwards` (TEST.129 area) |
| Orbital phase `mean_anomaly_deg`, `orbit_phase_deg` | DOUBLE degrees | 5.7e-14 deg at 360 | a planet's daily step of about 1 deg | accumulated in place: `MOD(phase + turn, 360)` per step; 1e6 steps cost about 3.5e-11 rad (`orbital-solvers-and-integrators.md` 504) | OK for years; WATCH for runs of 1e9 steps | `test_phase_accumulation_over_a_million_steps` |
| Phase formula (item 7 above) | Python, SQL, JavaScript | | | the add-then-reduce form in Python and JavaScript | WORK | `test_phase_forms_agree_for_a_long_elapsed_time` |
| Orbital elements, modified equinoctial | float | round trip 1e-16 | any conic | none | OK | existing kepler tests |
| Kepler solver tolerance | Python `KEPLER_XTOL` 1e-14, `UNIVERSAL_TOLERANCE` 1e-13; browser Newton 1e-10 | | an orbit position to a kilometre in a 40 AU orbit needs 1e-10 relative | browser tolerance is looser than the server's by 1e4 | WATCH | `test_js_kepler_matches_python_to_a_kilometre` |
| Galactic rotation per nightly step | `_rotate_about_axis`, `_galactic_turn` (fmod first) | | a step that moves 0.01 mpc | when the angle is under 1.05e-8 rad, `cos` rounds to exactly 1, so each step scales the radius by 1 + angle squared over two (7.8e-19 at 8 kpc for the 0.01 mpc step, about 2e-11 per Gyr, always the same sign) | OK, negligible, pinned | `test_rotation_does_not_drift_the_radius` over 1e6 steps |

### 3.3 Fields in the database

| Quantity | Type | Verdict | Note |
|---|---|---|---|
| `phenomenon_scatter.map_luminosity_sol` | FLOAT (7 digits), the only FLOAT column | WATCH | a display luminosity; never used as an input |
| `species.tech_level` | DECIMAL(3,2) | OK | exact by design; the column is on `species`, not `polities` (Decided, Boss, contradiction sweep card, 7:30 AM EDT 2026-10-21: no schema change) |
| Everything else numeric | DOUBLE or integer | OK | |

### 3.4 Sums and the field kernel

| Quantity | How computed | Spacing / error | Verdict | Test |
|---|---|---|---|---|
| Kernel accelerations and potentials (`physics/field_kernel.py`, numba, no fastmath) | `dx = positions[j] - px` in absolute galactic metres; plain sequential `+=` over at most 32 relevant masses | each position has 1,049 km spacing in metres at the reach (32.8 km at the Sun), so `dx` between two stars a parsec apart (3.1e16 m) has a relative error of 1e-12 at the Sun and 3e-11 at the reach; between two bodies of one system (1e11 m) it would be 3e-7 at the Sun and 1e-5 at the reach | WATCH: correct for the inter-star pulls it is used for (neighbour reach 2.5 sector edges); **must not be used for bodies of one system without subtracting a local origin first** | `test_kernel_differenced_in_a_local_frame_matches_absolute_within_the_spacing` at 8 kpc and 262 kpc |
| Kernel integration `t = t + h` | accumulated | a step of 1e9 s over 2,000 steps: relative 1e-16 per add | OK | none |
| `integrated_masses` | Python `sum` over 216,000 terms (bulge 40 cubed, halo steps 60 cubed) | left-to-right, relative error up to about 2e-11 | OK for normalising a density (a ratio); WATCH if the sum is ever compared for equality | none |
| Cumulative sums in the scatter draw (`np.cumsum`) | numpy | pairwise-ish, relative 1e-15 | OK | none |
| `block_stats`, `galaxy_block_stats` sums (`age_sum`, `luminosity_sol`) | SQL aggregate | order not fixed by the query | WATCH: totals can differ in the last digits between runs of the same data | `test_block_totals_do_not_depend_on_row_order` |
| Component density (`component_field.py`) | numpy pairwise sum over 6 cubed or 12 cubed samples | 1e-16 | OK | none |
| Neighbour search (`galaxy/sector.py` float64 `einsum`, `_CLEAR_MARGIN` 1e-9, exact `math.dist` fallback; `nav_graph.py` KD-tree on absolute pc) | float64 | an absolute pc difference of two stars 1 pc apart at the reach has 2.9e-11 relative to 1 pc spacing (898 km out of 3e13 km) | OK | none |
| Region test (`galaxy/region.py`, law of cosines `cos_delta`, `_EPS` 1e-9) | float64 | `acos` near 1 loses half its digits: an angle under 4e-5 rad is lost to the 1e-9 epsilon (about 10 pc at 262 kpc) | WATCH: a region smaller than about 10 pc at the reach should use the haversine form | `test_region_radius_of_one_parsec_at_the_reach` |
| Sector containment (`SectorCell.contains`, tolerance 1e-12 times the outer radius, 1e-12 rad) | float64 | 2.6e-7 pc at the reach | OK | existing sector tests |
| Mass and number sums over a sector in `sectors` rows | integer or DOUBLE | | OK | |

### 3.5 Camera, render path and wire formats

| Item | Representation | Spacing | First place error shows | Verdict | Test |
|---|---|---|---|---|---|
| Browser stars (`galaxymap3d.js` MAP.102) | `star - Math.round(target)` into a Float32Array; the points object sits at that frame | float32 on the offset only | none | OK | existing |
| Browser territory dots, superbubble lines, "made" ring points, block prisms (`galaxyprisms.js`, `galaxyblocks.js`, `makeBlockMesh`), `galaxystageview.js`, system view dots and belt particles | Float32Array of **absolute** coordinates | 4.9e-4 pc (100 AU) at 8 kpc, 3.1e-2 pc (6,400 AU) at 262 kpc | the vertices jitter when the camera is near them; whether a mesh translation cancels the absolute vertices was not verified for the prisms | WORK: item 5.3 | `test_browser_float32_attributes_hold_camera_relative_values` (a browser test that reads the attribute range) |
| Depth buffer | `logarithmicDepthBuffer: true`; far plane about 6 times the view radius | | z-fighting at the far plane | OK | none |
| GLSL precision | no `precision` or `highp` declarations in `nearfield.js` | float32 on most GPUs, mediump on some phones | star shader at near range | WATCH | none (not testable without a device) |
| Tile wire format (`web/lib/tilepack.py`) | per axis a float64 min and step plus 16-bit codes: step = (max - min) / 65,535; luminosity 16 bits over log10 -8 to 8; temperature and radius 8 bits | at the finest 16 pc tile 2.4e-4 pc (computed) | star position shown 0.0002 pc off | OK for display | `test_tile_quantisation_error_below_half_a_step` |
| Tile root | plus or minus 262,144 pc (524,288 pc edge, MAP.239) | | every address in the reach is in a tile | OK |  `test_tiles_cover_the_reach` |
| Tile cache `_trim_star` | x, y, z to 3 decimals (1 mpc), luminosity 4 digits, radius 3 digits, temperature 10 K | 1 mpc is 3.1e10 km, larger than any threshold | stars closer than 1 mpc merge in the cache | OK (display) | none |
| API rounding (`api/routes.py`, `db/query.py`, `web/maps`) | 1 to 5 decimals, `%.3g`, `%.4g` | | a displayed number, never stored | OK | none |
| JSON numbers | Flask `jsonify` writes the shortest round-trip text; a JavaScript Number is the same double | lossless | none | OK | none |
| `repr + "e0"` bulk literals (`db/store.py` `_row_literal`) | shortest round-trip | lossless | none | OK | |
| MySQL DOUBLE read-back text | not verified (17 digits or shortest) | | | NOT VERIFIED | `test_double_round_trips_through_the_driver` |
| Browser twins of Python | `galaxyprisms.js`, `galaxystages.js`, `galaxytilecut.js`, `orbitpositions.js`, `orbitclock.js`, `starmagnitude.js`, `starlight.js`, `nearfield.js`, `tilepack.js`, `numberformat.js`, `unitladder.js` | libm trig may differ from Python's in the last 1 to 11 ulp | an angle off by 1e-15 rad is 1e-8 pc at the reach; not visible | OK | `tests/test_js_unit.py` already compares some; add the phase and Kepler pair (item 5.4) |
| Browser clock `Date.now()` in milliseconds | | 1 ms | a phase moves by 1e-9 deg | OK | none |

### 3.6 Identifiers and keys

| Item | Representation | Verdict | Note |
|---|---|---|---|
| Object ID | 96-bit integer (shell 19, biased band 18, wedge 19, serial 28 with the top 2 bits the kind, body 12), stored BINARY(12) big-endian, Python int | OK | exact; the shell field holds 524,287, the reach needs 65,535 |
| Sector uid | 56-bit BIGINT UNSIGNED; unplaced is `1 << 56` plus the row id | OK | |
| Cell key (Morton, 0.25 pc cube, 21 bits per axis, 63 bits) | `floor((c + 262144) / 0.25)`; SQL generated column and `galaxy/cellkey.py` | WORK at one edge: the double just under 262,144 (`2^18 - 2^-35`) passes `ABS(x) < 262144`, but `x + 262144` rounds (ties to even) to 524,288.0, so the index is `2^21`, which overflows the 21-bit field and wraps to the opposite edge's cell. Python `cell_index` has the same arithmetic. No object is generated that close to the boundary today (the wrap needs a coordinate within 2.9e-11 pc of the reach), so the cost is zero; the fix is a one-line clamp | `test_cell_key_at_the_reach_edge` (item 5.7, done in DB.65: the index is clamped at `2^21 - 1`, migration 0106) |
| Scatter cell key | computed from `position_mpc / 1000.0` | WATCH: an inexact divide can move a point within one ulp of a cell boundary into the neighbouring cell compared with a key computed from parsecs directly; the effect is a point listed in the next 0.25 pc cell, which the region query already widens by a margin (4 mpc) | `test_scatter_key_agrees_with_the_parsec_key_off_boundary` |
| Cell versus spacing | 0.25 pc cell against 2.9e-11 pc spacing | OK | 8.6e9 spacings per cell; the division by 0.25 is exact |

## 4. Where something wider than float64 would be used

Nothing is stored wider than 64 bits and nothing needs to be. The cases a wider sum would help are all computation, not storage:

1. **Long integrations from the galactic centre** (**decided: wait for GEN.109**; the accumulator is built with its first user) (an orbit integrated over many Gyr in one go, GEN.109 and the galactic-potential work): `galactic-potential.md` already compares float64 roundoff over 1e6 daily steps (9e-12 pc) with 80-bit arithmetic. The recommended tool is a compensated (Kahan or `math.fsum`) accumulation of the position increment, with the position kept as a coarse part (integer parsecs) plus a fractional part, so the accumulated sum keeps a 28 km spacing from being the limit. `mpmath` is a diagnostic check (run once offline on a few orbits), not a production path.
2. **Ill-conditioned encounter solves** (a swept-sphere discriminant `B^2 - 4AC` with terms 1e42 and 1e23): `collisions-and-mergers.md` already specifies exact arithmetic; keep it and do not move the rest.
3. **Reproducibility across machines** (`generation-determinism.md`): a `detmath` module with `fsum` is proposed there and does not exist; it belongs to the determinism work, not this audit.

None of these changes what is stored. There is no case for 128-bit or 256-bit columns.

## 5. Cases that need work (proposed items; TODO files them)

Boss decided four of them on cards at 7:36 PM to 7:38 PM EDT on 2026-10-16: 5.1 document and test only, 5.6 raise the tile root to the reach, 5.3 fix the browser buffers now, and section 4's compensated sums wait for GEN.109. The rest are defaults.

| # | Case | What to change | Test |
|---|---|---|---|
| 5.1 | The nightly rotation recomputes a local offset as `absolute - centre` (section 3.1) | **Decided: document and test.** The offset is good only to the absolute spacing (28 km at the Sun, 898 km at the reach), 1e7 and 3e5 times below the smallest stored move; revisit if generating anywhere needs finer offsets | `test_local_offset_after_a_rotation_is_within_one_ulp_of_the_absolute_point` |
| 5.2 | The kernel differences absolute metres | Subtract a local origin (the sector centre or the heaviest mass) before the loop, in `physics/sector_path.py`, before any code uses the kernel for bodies of one system | `test_kernel_differenced_in_a_local_frame_matches_absolute_within_the_spacing` |
| 5.3 | Browser float32 attributes hold absolute coordinates for territory dots, superbubbles, the "made" rings, block prisms, stage views and the system view | **Decided: fix now (map lane).** Render them camera-relative like the stars (subtract a rounded camera target in double precision when the buffer is built, move the object by that target) | browser test reading the attribute range |
| 5.4 | Phase and Kepler differ between Python, SQL and JavaScript (section 2, item 7; tolerances 1e-10 against 1e-14) | One reduce-first form everywhere, and a JavaScript tolerance of 1e-12 | `test_phase_forms_agree_for_a_long_elapsed_time`, `test_js_kepler_matches_python_to_a_kilometre` |
| 5.5 | `orbital-updates.md` and `orbital-solvers-and-integrators.md` state 32.8 km at 8.128 kpc and integer milliparsecs | Done in this change: both now state the figures of section 1 and point here | |
| 5.6 | The map tile root covers only 32.8 kpc | **Done (MAP.239): the root is the reach, `TILE_ROOT_EDGE_PC` = 2 × `cellkey.REACH_PC` in `viewport.py`; the browser takes it from the page data** (map lane) | `test_tiles_cover_the_reach` |
| 5.7 | The cell key wraps just under the reach | **Done (DB.65, migration 0106):** the index is clamped to `2^21 - 1` in `cell_index`, `cell_keys` and the SQL generated expression | `test_cell_key_at_the_reach_edge` |
| 5.8 | The region test uses the law of cosines | Haversine form for angles under 1e-3 rad | `test_region_radius_of_one_parsec_at_the_reach` |

Each item is small and independent. None blocks GEN.294: its point columns are DOUBLE in the frames the audit confirms (section 3.1), and GEN.336 and GEN.337 add no column of another type.

## 6. Stale statements fixed or to fix

The two documents in 5.5 were corrected with this page. The stale comments in `web/lib/tilepack.py` line 30 and `galaxyprisms.js` line 475, the `point_columns.FRAMES` wording and the `generation-determinism.md` quantisation paragraph are listed as part of items 5.4 and 5.6 for whoever touches those files; the audit changes no code.

## 7. What was not verified

The text form MySQL returns for a DOUBLE and whether it round-trips exactly (3.5); the order of SQL aggregates (3.4); whether a mesh translation cancels the absolute prism vertices (3.5); GLSL precision on a real device; `physics/mathcheck.py` tolerances; the figure of 1.8e-12 pc in `galactic-potential.md` lines 373-378 (the spacing at 8.128 kpc is 9.1e-13 pc, so that figure is an error bound, not a spacing); the largest planet distance from its star. Spacings in section 1, the comet aphelion, the system extent, the tile step, the radial growth of the rotation and the cell-key wrap were worked out with numpy for this audit and not read from repository output.
