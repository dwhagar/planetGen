# Lagrange trojans (GEN.1149 module; GEN.1150, GEN.1151 and MAP.354 wire it)

Boss, 10-30: trojan swarms in the scatter, whole bodies at Lagrange points where the mass rule allows, and an icon at the centre of each point with objects, labelled `[secondary]-L#` (for example `Jupiter-L4`). The research and sources are in `reports/lagrange-trojans-2026-10-30.md`; this note is what the code does.

## What GEN.1149 builds (lock-free, no generator change)

`physics/trojans.py` and `physics/data/trojans.json`, stand-in rows 6.121 TJ1 to TJ5:

- `stable`, `max_trojan_mass`, `allows`: L4 and L5 hold a trojan when the secondary plus the trojan stay under 0.0385 of the three masses (Gascheau-Routh; Schwarz et al. 2016). L1, L2, L3 hold nothing.
- `libration_period`, `motion`: the slow libration of a tadpole, and tadpole against horseshoe by amplitude.
- `swarm_split`, `leading_share`: the scatter's `trojans` share into L4 and L5 regions.
- `point_label`, `pair_text`: `Jupiter-L4` and `L4 of Sun and Jupiter`, for the map.
- `longitude_offset_deg`: +60 for L4, -60 for L5, so a stored trojan has the secondary's `a`, `e` and `i` near it and `mean_anomaly_deg` = the secondary's plus the offset.

## Wiring (needs an epoch; not in this item)

- **GEN.1150** (train 8, after GEN.1006): the generator turns the scatter `trojans` share into a leading and a trailing region at the L4/L5 of a giant or brown dwarf; new columns (host id, point, libration amplitude and phase; `db/column_units.py`); the orbit updater adds the libration as an overlay; GEN.1006 treats a trojan as a member, not an overlap.
- **GEN.1151** (train 8 or 9): whole planets and moons (Class M and Class J allowed by the mass rule) at L4/L5, only from the disc formation model; no rate is set (the observed searches only bound it); W7 folds in.
- **MAP.354**: the System Map icon, one per point that has objects, at its centre, labelled with `point_label`.

## Limits

Eccentricity and inclination limits (TJ5) and the leading share (TJ4) are stand-ins at L4 with no source read; an N-body run replaces them. The finite-trojan rule (TJ2) is the 4 percent figure TROY III uses; its own text is not read.

## Pure parts built ahead of the wiring (GEN.1151, MAP.354)

- `physics/trojan_bodies.py` (rows 6.125 TB1 to TB4): `case`, `room_earth`, `classes_allowed`, `fits`, `can_form`. A Sun and a Jupiter admit Class J and Class M; a 0.1 solar-mass dwarf admits a Class J trojan only if its giant is under about 4 Jupiter masses; two stars admit nothing; a Jupiter-mass planet holds about 11.7 Earth masses of trojan beside a 1 Earth-mass moon. The formation gate names its inputs for the oligarch chain of `planet_architecture.py`: `solids_in_zone_earth` (the solids it leaves in the secondary's co-orbital zone) and `disc_lifetime_myr`.
- `physics/lagrange_geometry.py` (no stand-ins): `points` gives L1 to L5 for a pair in the system frame (L1 to L3 by bisection of the force balance, L4 and L5 equilateral, L4 ahead along the motion), `occupied` keeps those with objects and labels them `Jupiter-L4`. Shared with the gravity layer MAP.169.
