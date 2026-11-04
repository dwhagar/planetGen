# planetgen/db/orbital_events.py

"""
The orbital event writer (GEN.197; docs/design/orbits/collisions-and-mergers.md sections 5 and 7, tables of DB.56)
=================================================================================================================

Three small parts that the orbital update calls in this order when a body's update crosses another body's stored path:

1. `find_contact` asks `physics.closest_approach` the one question (do the two paths come within `R1 + R2` in the pair's window?) and gives
   back a `Contact`: the time of the hit, the miss distance (the impact parameter), the separation at the start of the window, the window used
   and whether the curvature cap cut it.
2. `resolve_contact` hands the two `Party`s and the contact to `physics.collision_outcome.resolve` and gives back its `Outcome`.
3. `record_event` writes the outcome as one `orbital_events` row and its `orbital_event_bodies` rows (the participants, then the products), works out
   the conservation check from the outcome itself, and writes one `GEN orbital_event` line to the activity log. `event_row` is the same without the
   database, for a caller that wants the row only.

**The ID rules** (`product_uids`, Boss 2026-10-09: the 80-bit birth-location part of the 96-bit object ID stays with the body it was born with):

- a body that survives keeps its ID, and a merger keeps the heavier body's (the outcome names the survivor by the larger parent's key), so the lighter
  one is retired;
- an ejected planet keeps its ID (`ejected_uid` is the identity, kept as a name so that a caller and a test say what they mean) and becomes a rogue
  planet row, which is `store.insert_rogue_planet`'s to write when the unbind of GEN.291 calls this module;
- a split, and the field a collision makes, give every new body a new run-time serial of the sector the event happened in
  (`store._runtime_uids`), and the counters only move up, so the ID of a destroyed body is never handed out again.

Nothing here draws a random number: the same bodies and window always give the same outcome and the same row. The orbital update that calls
`find_contact` for every pair of a body's influence set (the broad phase of section 5.4) is not here.
"""

import dataclasses
import math

from planetgen.admin import activity_log
from planetgen.galaxy import object_uid
from planetgen.physics import closest_approach, collision_outcome, constants

M_PER_KM = 1000.0
M_PER_PC = constants.PARSEC_M
"""float: Metres in a parsec; the bodies of an event are stored in pc and km/s, the outcome works in metres and m/s."""

EVENT_KINDS = ("collision", "merger", "split", "ejection", "capture", "unbind", "warning")
"""tuple: The `kind` values `orbital_events` accepts (its check constraint)."""

MERGER_CODES = ("merger", "graze_and_merge", "absorbed", "giant_merger", "stellar_merger")
"""tuple: The outcome codes in which the lighter body becomes part of the heavier: the event is a `merger` and the lighter body's fate is `absorbed`."""

CONSERVATION_TOLERANCE = 1e-9
"""float: The relative mass and momentum error `record_event` stores as acceptable. The outcome conserves both exactly bar rounding, so a larger error is a bug
and the row carries it (and the `numeric_warning` flag) for the admin report to show."""

HEAT_SHARE = 1.0
"""float: The share of the impact kinetic energy the report assumes is deposited as heat (section 7); a stand-in until the heating of an impact is modelled,
so the report states that it is an assumption, not a result."""

NUMERIC_RESOLUTION_M = 28.0e3
"""float: The coarsest position resolution of a stored point in a sector (14 to 28 km, section 7); a pair whose `R1 + R2` is under 100 times this carries the
`numeric_warning` flag."""


@dataclasses.dataclass(frozen=True)
class Party:
    """A body that takes part in an event: its `body` for the outcome, its stored 96-bit `uid` (bytes) and the words and place the report shows."""

    body: collision_outcome.Body
    uid: bytes
    name: str = ""
    object_class: str = ""
    sector_id: int | None = None
    admin_edited: bool = False


@dataclasses.dataclass(frozen=True)
class Contact:
    """A hit found in a pair's window: when (`t_hit_s`, the caller's clock), how near the paths pass (`impact_parameter_m`, `0 <= b <= R1 + R2`) and the window."""

    t_hit_s: float
    impact_parameter_m: float
    separation_start_m: float
    window_start_s: float
    window_end_s: float
    window_capped: bool
    already_overlapping: bool


def key_of(uid):
    """The orderable name of a body for `collision_outcome` (the printed object ID of its 12 bytes)."""
    return object_uid.format_id(object_uid.from_bytes(bytes(uid)))


def _norm(vector):
    return math.sqrt(sum(c * c for c in vector))


def _miss_distance(dp, dv):
    """How near the straight relative path `dp + dv t` passes, `t >= 0`: the distance at the closest approach, or `|dp|` when the pair is moving apart."""
    speed2 = sum(c * c for c in dv)
    if speed2 == 0.0:
        return _norm(dp)
    t = -sum(p * v for p, v in zip(dp, dv)) / speed2
    if t <= 0.0:
        return _norm(dp)
    return _norm(tuple(p + v * t for p, v in zip(dp, dv)))


def find_contact(x, y, radius_sum_m, previous_epoch_x_s, now_x_s, omega_s=0.0, edge_m=math.inf):
    """
    Whether body X's update from `previous_epoch_x_s` to `now_x_s` carries it within `radius_sum_m` of body Y's stored path, and the `Contact` if so (`None` if not,
    or if Y's segment starts after the update's end). `x` and `y` are `closest_approach.Segment`s. The window is `[max(e_X_prev, e_Y), e_X_now]`; where the galaxy's
    curvature (`omega_s`, the angular speed in 1/s of the galactic orbit, and `edge_m` the sector edge) caps the straight window, the window is cut and the
    contact's `window_capped` says so. A pair that already overlaps at the window's start is a hit at that time with `already_overlapping` set: the caller
    leaves out pairs made in the same run.
    """
    window = closest_approach.pair_window(y, previous_epoch_x_s, now_x_s)
    if window is None:
        return None
    start, end = window
    state_at = closest_approach.linear_state_at(x, y)
    dp0, dv0 = state_at(start)
    v_rel = _norm(dv0)
    cap = closest_approach.window_cap(radius_sum_m, v_rel, omega_s, edge_m) if math.isfinite(edge_m) else math.inf
    approach = closest_approach.first_hit(state_at, start, end, radius_sum_m, cap_s=cap)
    if not approach.hit:
        return None
    # The miss distance of the straight path: from the relative state at the window's start, not at the hit (the hit is on the way in).
    b = min(_miss_distance(dp0, dv0), radius_sum_m)
    return Contact(approach.t_hit_s, b, _norm(dp0), start, end, bool(end - start > cap), _norm(dp0) <= radius_sum_m)


def resolve_contact(a, b, contact, always_destroy_terrestrial=None):
    """The `collision_outcome.Outcome` of parties `a` and `b` for a `contact` (the seed of its draws starts at the window's start)."""
    for party in (a, b):
        if party.body.key != key_of(party.uid):
            raise ValueError(f"the body's key {party.body.key!r} is not its object ID {key_of(party.uid)!r}")
    return collision_outcome.resolve(a.body, b.body, contact.impact_parameter_m, contact.window_start_s, always_destroy_terrestrial)


def ejected_uid(uid):
    """An ejected planet keeps its ID (it becomes a rogue planet row; the ID is the same 96 bits)."""
    return bytes(uid)


def product_uids(outcome, parties, new_uids=()):
    """
    The ID of every product of `outcome`: `{key: uid}` for the survivors (the larger parent's ID in a merger, each parent's own where both survive) and, as the
    last entry, `"field"` for the debris field if there is one. `new_uids` is a list of the new run-time IDs the caller reserved
    (`store._runtime_uids`), one for the field; it is read only when the outcome has a field.
    """
    by_key = {key_of(party.uid): bytes(party.uid) for party in parties}
    uids = {}
    for survivor in outcome.survivors:
        if survivor.key not in by_key:
            raise ValueError(f"the outcome names a survivor {survivor.key} that is not one of the parties")
        uids[survivor.key] = by_key[survivor.key]
    if outcome.debris is not None:
        if not new_uids:
            raise ValueError("the outcome has a field, which needs a new object ID")
        uids["field"] = bytes(new_uids[0])
    return uids


def _participant_fate(outcome, key):
    if key in {survivor.key for survivor in outcome.survivors}:
        return "survived"
    return "absorbed" if outcome.code in MERGER_CODES else "destroyed"


def _body_row(role, ordinal, uid, kind, fate, name, object_class, mass_kg, radius_m, composition, sector_id, position_m, velocity_ms, admin_edited, note=None):
    return {
        "role": role, "ordinal": ordinal, "uid": bytes(uid), "object_kind": kind, "fate": fate, "name": name or None,
        "object_class": object_class or None, "mass_kg": mass_kg, "radius_km": None if radius_m is None else radius_m / M_PER_KM,
        "composition": composition, "sector_id": sector_id,
        "x_pc": position_m[0] / M_PER_PC, "y_pc": position_m[1] / M_PER_PC, "z_pc": position_m[2] / M_PER_PC,
        "vx_kms": velocity_ms[0] / M_PER_KM, "vy_kms": velocity_ms[1] / M_PER_KM, "vz_kms": velocity_ms[2] / M_PER_KM,
        "admin_edited": 1 if admin_edited else 0, "note": note,
    }


def _conservation(outcome, parties):
    """`(mass_error, momentum_error)`, relative: the products and what left with no field against the two bodies, mass and momentum."""
    mass_in = sum(party.body.mass_kg for party in parties)
    momentum_in = [sum(party.body.mass_kg * party.body.velocity_ms[i] for party in parties) for i in range(3)]
    mass_out = sum(s.mass_kg for s in outcome.survivors) + outcome.lost_mass_kg
    momentum_out = [sum(s.mass_kg * s.velocity_ms[i] for s in outcome.survivors) for i in range(3)]
    if outcome.debris is not None:
        mass_out += outcome.debris.mass_kg
        momentum_out = [momentum_out[i] + outcome.debris.mass_kg * outcome.debris.velocity_ms[i] for i in range(3)]
    scale = max(_norm(momentum_in), mass_in * max(_norm(party.body.velocity_ms) for party in parties), 1e-300)
    return abs(mass_out - mass_in) / mass_in, _norm(tuple(o - i for o, i in zip(momentum_out, momentum_in))) / scale


def event_row(outcome, parties, contact, detected_time_unix, kind=None, run_id=None, containing_uid=None, influence_body_count=None,
              sectors_dirtied=None, suppressed_duplicates=0, focusing_radius_used=False, product_ids=None, window_unix_offset_s=0.0):
    """
    `(event, bodies)` for `outcome`: the `orbital_events` columns by name and the list of `orbital_event_bodies` rows by name (less `event_id`).
    `contact` times are on the caller's clock; `window_unix_offset_s` is added to turn them into unix seconds (0 when the caller's clock is unix time).
    `product_ids` is `product_uids`' answer (needed only for the bodies' rows; the event row does not need it).
    """
    if kind is None:
        kind = "merger" if outcome.code in MERGER_CODES else "collision"
    if kind not in EVENT_KINDS:
        raise ValueError(f"unknown event kind {kind!r}")
    numbers = outcome.numbers
    mass_error, momentum_error = _conservation(outcome, parties)
    r_sum = sum(party.body.radius_m for party in parties)
    heavier = max(parties, key=lambda party: (party.body.mass_kg, key_of(party.uid)))
    impact_energy = 0.5 * (parties[0].body.mass_kg * parties[1].body.mass_kg / (parties[0].body.mass_kg + parties[1].body.mass_kg)) * numbers["v_imp"] ** 2
    merged = outcome.survivors[0] if outcome.code in MERGER_CODES and outcome.survivors else None
    resolution_bad = NUMERIC_RESOLUTION_M > 0.01 * r_sum
    event = {
        "run_id": run_id, "kind": kind, "outcome_code": outcome.code, "regime": regime_of(numbers),
        "event_time_unix": contact.t_hit_s + window_unix_offset_s, "detected_time_unix": detected_time_unix, "seed_string": outcome.seed,
        "separation_start_km": contact.separation_start_m / M_PER_KM, "separation_impact_km": r_sum / M_PER_KM,
        "impact_parameter_km": contact.impact_parameter_m / M_PER_KM,
        "critical_impact_parameter_km": numbers["b_crit"] * r_sum / M_PER_KM,
        "v_inf_kms": numbers["v_inf"] / M_PER_KM, "v_esc_kms": numbers["v_esc"] / M_PER_KM, "v_imp_kms": numbers["v_imp"] / M_PER_KM,
        "q_r_j_per_kg": numbers["q_r"], "q_star_rd_j_per_kg": numbers["q_star"], "erosion_x": numbers["x"],
        "interacting_mass_fraction": numbers["lens_fraction"],
        "catastrophic_velocity_kms": None,
        "window_s": contact.window_end_s - contact.window_start_s, "window_capped": 1 if contact.window_capped else 0,
        "focusing_radius_used": 1 if focusing_radius_used else 0,
        "mass_error": mass_error, "momentum_error": momentum_error, "conservation_tolerance": CONSERVATION_TOLERANCE,
        "impact_energy_j": impact_energy, "heat_share": HEAT_SHARE, "lost_mass_kg": outcome.lost_mass_kg,
        "merged_mass_kg": None if merged is None else merged.mass_kg,
        "new_radius_km": None if merged is None or merged.radius_m is None else merged.radius_m / M_PER_KM,
        "spin_breakup_ratio": None, "deuterium_limit_mj": None, "hydrogen_limit_sol": None,
        "borderline": 1 if "borderline" in outcome.flags else 0,
        "containing_uid": None if containing_uid is None else bytes(containing_uid),
        "influence_body_count": influence_body_count, "sectors_dirtied": sectors_dirtied, "suppressed_duplicates": int(suppressed_duplicates),
        "already_overlapping": 1 if contact.already_overlapping else 0,
        "numeric_warning": 1 if (resolution_bad or mass_error > CONSERVATION_TOLERANCE or momentum_error > CONSERVATION_TOLERANCE) else 0,
        "flags": ",".join(outcome.flags),
    }
    bodies = []
    for ordinal, party in enumerate(sorted(parties, key=lambda p: key_of(p.uid))):
        body = party.body
        bodies.append(_body_row("participant", ordinal, party.uid, body.kind, _participant_fate(outcome, key_of(party.uid)), party.name,
                                party.object_class, body.mass_kg, body.radius_m, body.composition, party.sector_id, body.position_m,
                                body.velocity_ms, party.admin_edited, "heavier" if party is heavier else None))
    ids = product_ids or {}
    products = [(s.key, s) for s in outcome.survivors]
    ordinal = 0
    for key, survivor in products:
        origin = next(party for party in parties if key_of(party.uid) == key)
        bodies.append(_body_row("product", ordinal, ids.get(key, origin.uid), survivor.kind, "survived", origin.name, origin.object_class, survivor.mass_kg,
                                survivor.radius_m, survivor.composition, origin.sector_id, survivor.position_m, survivor.velocity_ms,
                                origin.admin_edited, ",".join(survivor.flags) or None))
        ordinal += 1
    if outcome.debris is not None:
        field = outcome.debris
        bodies.append(_body_row("product", ordinal, ids["field"], "field", "created", "", "asteroid field", field.mass_kg, field.radius_m, field.composition,
                                heavier.sector_id, field.position_m, field.velocity_ms, False,
                                None if not field.hero_mass_kg else f"hero {field.hero_mass_kg:.6g} kg"))
    return event, bodies


def regime_of(numbers):
    """The impact regime named in the report, from the numbers the outcome was decided from: `merging` under the disruption energy at a head-on pass, `grazing`
    from the critical impact parameter, `catastrophic` from `x = 1`, else `erosive`."""
    x, grazing = numbers["x"], numbers["b"] >= numbers["b_crit"]
    if x >= 1.0:
        return "catastrophic"
    if grazing:
        return "grazing"
    return "erosive" if x >= 0.1 else "merging"


def record_event(conn, outcome, parties, contact, detected_time_unix, new_field_uid=None, **context):
    """
    Writes `outcome` as one `orbital_events` row and its `orbital_event_bodies` rows, one activity-log line `GEN orbital_event`, and returns the event's id. The caller
    reserves the new ID of the field (`store._runtime_uids`) when the outcome has one and passes it as `new_field_uid`. `context` is `event_row`'s keyword arguments
    (`kind`, `run_id`, `containing_uid`, `influence_body_count`, `sectors_dirtied`, `suppressed_duplicates`, `focusing_radius_used`, `window_unix_offset_s`).
    The event and its bodies are written in the caller's transaction; the caller commits.
    """
    ids = product_uids(outcome, parties, [] if new_field_uid is None else [new_field_uid])
    event, bodies = event_row(outcome, parties, contact, detected_time_unix, product_ids=ids, **context)
    columns = list(event)
    cursor = conn.execute(
        f"INSERT INTO orbital_events ({', '.join(columns)}) VALUES ({', '.join('?' for _ in columns)})", tuple(event[c] for c in columns))
    event_id = cursor.lastrowid
    for row in bodies:
        names = ["event_id", *row]
        conn.execute(f"INSERT INTO orbital_event_bodies ({', '.join(names)}) VALUES ({', '.join('?' for _ in names)})", (event_id, *row.values()))
    activity_log.event("GEN", "orbital_event", id=event_id, kind=event["kind"], code=event["outcome_code"], regime=event["regime"],
                       bodies=",".join(key_of(p.uid) for p in parties))
    return event_id
