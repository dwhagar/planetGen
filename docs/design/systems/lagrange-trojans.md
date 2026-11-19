# Lagrange trojans (GEN.1149 module; GEN.1150, GEN.1151 and MAP.354 wire it)

Boss, 10-30: trojan swarms in the scatter, whole bodies at Lagrange points where the mass rule allows, and an icon at the centre of each point with objects, labelled `[secondary]-L#` (for example `Jupiter-L4`). The research and sources are in `reports/lagrange-trojans-2026-10-30.md`; this note is what the code does.

## What GEN.1149 builds (lock-free, no generator change)

`physics/trojans.py` and `physics/data/trojans.json`, stand-in rows 6.121 TJ1 to TJ5:

- `stable`, `max_trojan_mass`, `allows`: L4 and L5 hold a trojan when the secondary plus the trojan stay under 0.0385 of the three masses (Gascheau-Routh; [R3984]). L1, L2, L3 hold nothing.
- `libration_period`, `motion`: the slow libration of a tadpole, and tadpole against horseshoe by amplitude.
- `swarm_split`, `leading_share`: the scatter's `trojans` share into L4 and L5 regions.
- `point_label`, `pair_text`: `Jupiter-L4` and `L4 of Sun and Jupiter`, for the map.
- `longitude_offset_deg`: +60 for L4, -60 for L5, so a stored trojan has the secondary's `a`, `e` and `i` near it and `mean_anomaly_deg` = the secondary's plus the offset.

## Wiring (needs an epoch; not in this item)

- **GEN.1150** (train 8, after GEN.1006): the generator turns the scatter `trojans` share into a leading and a trailing region at the L4/L5 of a giant or brown dwarf; new columns (host id, point, libration amplitude and phase; `db/column_units.py`); the orbit updater adds the libration as an overlay; GEN.1006 treats a trojan as a member, not an overlap.
- **GEN.1151** (train 8 or 9): whole planets and moons (Class M and Class J allowed by the mass rule) at L4/L5, only from the disc formation model; no rate is set (the observed searches only bound it); W7 folds in.
- **MAP.354**: the System Map icon, one per point that has objects, at its centre, labelled with `point_label`.

## GEN.1150 slice 1: the swarm as two regions of belt rows

- **Mass.** `swept_hand_on.trojan_part`: the bound share of the mass an outer disk lost (what its outermost giant cleared and scattered, less what the giant threw out by its Safronov number), times the swept split's `trojans` row (0.002), times what the neighbours leave (`trojans.surviving_share`, GEN.1455). It is part of `other`, not added to it. `belt_birth.hand_on` puts it on the disk as `swarm_kg`.
- **Regions.** `generation/trojan_region.build` runs once the system is final (after the spacing pass, the irregular moons and the comets, `StarSystem._add_trojan_regions`), so the pass never sees a body sharing a giant's orbit. `trojan_swarm.draw_regions` splits the mass between L4 and L5 and draws each region's libration, mean eccentricity and mean inclination; the draw and the makeup are on the swarm's own streams (keyed by the giant), so no other draw moves. A region is a `TrojanRegion` (an `AsteroidBelt`): mass after the star's age of collisions, largest body, composition drawn at the giant's distance, orbit radii `a (1 -/+ e)`. Under `MIN_PRESENT_MASS_EARTH` a region is not kept. No members yet.
- **Storage.** A row of `asteroid_belts` whose parent is the giant (`parent_type` `planet`), with `lagrange_point`, `libration_amplitude_deg`, `libration_phase_deg`, `mean_eccentricity`, `mean_inclination_deg` (migration 0269). It is saved after its giant in the same list and read back beside it. Its name is `<giant> L4 trojan region`, worked out, never stored, and it takes no number among the star's belts.
- **How often.** No rate is set. A swarm needs the outermost giant to have swept an outer disk, and a giant of a high Safronov number throws nearly all of it out. In 160 G2V systems 19 held a giant, 2 had an outer disk the giant swept, and the swarms of those two were 1.1e22 kg (0.0019 Earth mass) and 3e7 kg (nothing): about 1 percent of G2V systems hold a swarm worth keeping before the neighbours' share. Jupiter's own Trojans (1e-5 Earth mass) are far below what that giant's sweep would leave, so the swarm law is not tuned to them.
- **Slice 2 (libration overlay and MAP.354).** `physics/trojan_overlay.py`: a region is stored with no position; its centre is the giant's longitude plus 60 (L4) or minus 60 (L5) degrees plus the libration `A sin(2 pi t / P_lib + phase)`, `t` the time since Unix time zero (so the stored phase is the phase at Unix time zero and the orbit updater writes nothing for a region). The system map (`web/maps/systemmap.py`) draws a region as an arc on the giant's orbit, as wide as the libration amplitude, and a diamond at its centre labelled `<giant>-L4` (`trojans.point_label`), on the whole-system scene and on a member or wide secondary's scene; it takes the giant's current direction, so it follows the giant as the updater moves it. The 3D scene still draws a region as a belt.
- **Not yet:** members of a region, regions of the giants of a close pair's members and of a circumbinary node, whole planets and moons at L4/L5 of a moon (GEN.1151 slice B, design only). A system's belt count includes the regions (they are belt rows).

## Whole planets at L4 and L5 of a giant (GEN.1151 slice A)

`generation/trojan_planet.py`, run by `StarSystem._add_trojan_planets` after the belt members are promoted and before the irregular moons, makes at most one planet a giant, from the giant's own stream:

- **Zone solids** (`trojan_bodies.zone_solids_earth`): the disk's solid surface density where the giant formed (`planet_mass.solids_surface_density_gcm2` with the star's formation mass, snow line and metals and the architecture's disk factor, the disk the planets and belts were made from) times the area of the two tadpole zones: `ZONE_HALF_WIDTH_HILL` (2) Hill radii of the giant either side of its orbit, `ZONE_AZIMUTH_SHARE` (one third of the circle). Rows TB5, TB6.
- **Ceiling** (`trojan_bodies.ceiling_earth`): `FORMATION_EFFICIENCY` (TB3, kept at 0.1) of the zone solids times the share of the growth time the gas disk allowed (`disk_life_factor`, TB4: nothing under 0.3 Myr, the whole at 3 Myr, a ramp between; it replaces the flat 1 Myr minimum that sat on the median 0.9 Myr dispersal of a G dwarf's disk), held to the room the mass rule leaves.
- **Mass**: none when the ceiling is under the lightest planet (`PLANET_MIN_MASS_EARTH`), else log-uniform from it to the ceiling. No rate is set.
- **Stay**: the giant's neighbours and companion stars leave a share of the zone stable (`trojans.surviving_share`, the one the swarm uses, GEN.1455), taken as the chance the body stays. A giant that moved more than a quarter of where it formed (`MAX_MIGRATION_SHARE`, TB7) holds none.
- **Point**: L4 with the leading share (`trojans.leading_share`, TJ4), else L5.
- **Orbit**: the planet generator is given `trojan_of=(giant, point)` (reused from the belt promotion: a stream of its own, `Planet(mass=...)`): the giant's plane, node, periapsis and eccentricity, and the mean anomaly of the giant plus or minus 60 degrees (`trojan_bodies.mean_anomaly_rad`). The state vector is the body, so the orbit updater keeps the pair 60 degrees apart. A planet has no libration (the swarm regions do).
- **Stored**: `planets.lagrange_point` and `planets.lagrange_host_id` (migration 0270 on this branch, 0272 at the train). The giant's row is saved first, the planet names it, and the loader links them.

The measure that set TB3, TB4 and the gate is in `reports/gen1151/` (G, K and M stars: the zone holds 1 to 3 Earth masses of solids, so a whole trojan is rocky or icy and under about 7 Earth masses; with a flat 1 Myr minimum and no ramp over a third of the giants beyond 0.5 AU could hold a planet-mass trojan, because the median G dwarf's disk sat on the minimum and the giant's core had taken solids of the zone that the measure counted as free; the ramp cuts the median ceiling to a fifth).

Slice B (a moon at L4 or L5 of a moon of a giant) needs a gate for the solids of a co-orbital zone in a circumplanetary disk, which nothing holds yet; it stays design only.

## Limits

Eccentricity and inclination limits (TJ5) and the leading share (TJ4) are stand-ins at L4 with no source read; an N-body run replaces them. The finite-trojan rule (TJ2) is the 4 percent figure TROY III uses; its own text is not read.

## Pure parts built ahead of the wiring (GEN.1151, MAP.354)

- `physics/trojan_bodies.py` (rows 6.125 TB3 to TB7): `case`, `room_earth`, `fits`, and since slice A the zone, the ceiling and the mass draw. The mass is tested alone, no class letter is read (class is decided last). A Sun and a Jupiter admit a giant-mass and a rocky-mass trojan; a 0.1 solar-mass dwarf admits a 35 Earth-mass trojan only if its giant is under about 4 Jupiter masses; two stars admit nothing; a Jupiter-mass planet holds about 11.7 Earth masses of trojan beside a 1 Earth-mass moon. The formation gate is `zone_solids_earth`, `disk_life_factor`, `ceiling_earth` and `draw_mass_earth` (the section below).
- `physics/lagrange_geometry.py` (no stand-ins): `points` gives L1 to L5 for a pair in the system frame (L1 to L3 by bisection of the force balance, L4 and L5 equilateral, L4 ahead along the motion), `occupied` keeps those with objects and labels them `Jupiter-L4`. Shared with the gravity layer MAP.169.

## Neighbours take the swarm (GEN.1455; system-interactions.md rule 5)

A giant's swarm of the swept disk used to be a flat 0.2 percent of the scatter whatever else was in the system. It now keeps only the part that every other giant and every companion star leaves stable (`trojans.surviving_share`, rows 6.121 TJ6 to TJ8):

- **Another giant:** the separation of the two over `HILL_SEPARATION_K` (3.46, the Hill stability of a pair) mutual Hill radii is `x`. The share is 0 at `x <= 1` and `1 - exp(-STABLE_RAMP (x - 1))` past it. The swarm keeps the product of the share against every other giant, so a second giant at twice the distance cuts a swarm that a lone giant keeps whole.
- **A companion star:** `Star.a_crit_au`, the stability limit the clearance ledger already holds, over `STAR_LIMIT_MARGIN` times the giant's distance is `x`, and the same ramp applies. A star with no companion has no limit and keeps the swarm whole.
- **What is lost** is not made again: `scatter_shares(trojan_survival=...)` adds the lost part to the ejected share, so the shares still add to one. A survival of 1.0 gives the shares of before to the last digit.

The census (`tests/audit/trojan_census.py`) generates each seed with and without the rule and counts the swarms that change; a system with one giant and no companion never does.

In 960 seeds the generator drew 22 outer-disk swarms; neighbours cut 17 of them (none of the 5 of a lone giant at a lone star), but all 22 had a swept share of 0.0 already, since the belt lay beyond the giant's zone. The swarm law is therefore ready for the swarm draw of GEN.1150 rather than moving stored numbers today.
