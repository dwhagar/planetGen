# planetgen/physics/effect_ledger.py

"""
The effect ledger (GEN.1465)
============================

For every body of a system and every effect that acts on it (light, UV, tide, heating, Kozai and secular forcing, precession and the perturbation of the
eccentricity), the size of each source's contribution as a share of the largest source of that effect on that body (the largest has 1.0). The ledger
changes no result: it records what the effect code computed so that "when a tiny effect starts to matter" can be read from the stored galaxy
(docs/design/generation/system-interactions.md, `reports/gravity-audit-2026-11-08.md` gap 2, `reports/gen1465-effect-ledger-design-2026-11-08.md`).

**Every source is computed, whatever its size** (even 1e-12 of the largest). Only the writing is limited by a floor, `tuning.EFFECT_LEDGER_FLOOR`
(1e-3, 0.1 percent, Boss's choice): a source whose share is under it is not a row of its own; its share is added to the body's summary row (`source_type` 0, `source_id` 0), which also
holds how many sources fell under the floor and the largest source's contribution in the effect's own unit. So the sum of the tail is never lost. A floor
under 1e-37 stores every source.

**In memory** a system holds an `EffectLedger` (`system.effect_ledger`): effect code by effect code, each receiving body with its list of
`(source, value)` pairs, bodies and sources being the generator's own objects. **On save** `store.insert_star_system` gives every inserted body's `(table, id)`
(`object_rows.insert_row` notes them) and `EffectLedger.rows(star_system_id, ids)` turns the contributions into the rows of the `effect_ledger` table.

**The five gravity effects are on demand** (`ON_DEMAND_EFFECTS`, Boss 8:37 PM EST 11-08): tide, heating, secular, precession and eccentricity store only the largest
source of each (body, effect) and the summary row (the rest's summed share, how many, the largest's value). Every source's share is recomputed from the stored
masses and orbits when a page or a census needs it (a point mass sum, a few milliseconds a system, `reports/gen1465-ledger-store-vs-compute-2026-11-09.md`);
only light and UV, whose sources are stars across the neighbourhood, keep rows above the floor.

**A lone source is not a row** (light and UV): a body with one star sending it light has no row, since that source is the body's host star and the flux is on
the body; the rows begin where a second source matters (a binary or a multiple). The sources of light and UV are the stars of the body's own system
(`radiation_field`), not the neighbouring stars of the galaxy, so "the largest" is the largest of the system and not of the sky.

**Edits and deletes:** deleting a system deletes its rows (`ON DELETE CASCADE`); an edit of a system (`edits.save_system_edits`) clears the system's rows, because the edited
body's shares are no longer the ones generation noted, and they stay empty until the system's ledger is refilled; a standalone phenomenon deleted alone takes its rows as
receiver and source (`edits.forget_ledger_rows`).

The effects and the unit each is stored in are `EFFECTS`; a lane that builds an effect calls `EffectLedger.add` for every (receiver, source) pair it
computes, with the contribution in that unit (any monotone measure of the effect's size on the receiver works, the share only compares sources of one
effect on one body).
"""

import math

from planetgen import tuning

EFFECTS = (
    ("light", 1, "W_m2", "the flux a source puts on the body"),
    ("uv", 2, "W_m2", "the ultraviolet flux a source puts on the body"),
    ("tide", 3, "m_s2", "the tide-raising acceleration a source gives the body"),
    ("heating", 4, "W", "the tidal heating power the source's forcing puts into the body"),
    ("secular", 5, "1/s", "the rate at which a source forces the body's eccentricity and inclination (Kozai and secular terms)"),
    ("precession", 6, "rad_s", "the rate at which a source turns the body's orbit (apsidal and nodal precession)"),
    ("eccentricity", 7, None, "the eccentricity a source adds to the body's orbit"),
)
"""tuple: `(name, code, unit key of the unit registry or None for a plain number, what it measures)` for each effect; the code is what the table stores."""

ON_DEMAND_EFFECTS = frozenset({"tide", "heating", "secular", "precession", "eccentricity"})
"""frozenset: The effects whose sources are recomputed from stored masses and orbits and not stored: their (body, effect) keeps the largest source and the summary row only."""

EFFECT_CODES = {name: code for name, code, *_rest in EFFECTS}
"""dict: Each effect's name to its code."""

EFFECT_UNITS = {name: unit for name, _, unit, _ in EFFECTS}
"""dict: Each effect's name to its unit key (None for a plain number)."""

BODY_TABLES = ("stars", "planets", "moons", "asteroid_belts", "comets", "black_holes", "neutron_stars")
"""tuple: The kind tables whose rows are receivers or sources of the ledger, in the order of their codes (1 to 7); every one is in `db/kind_tables.py`."""

BODY_CODES = {table: index + 1 for index, table in enumerate(BODY_TABLES)}
"""dict: A body's kind table to its code in the ledger (0 is the summary row's source)."""

ROW_FIELDS = ("star_system_id", "body_type", "body_id", "effect", "source_type", "source_id", "share", "n_below", "largest_value")
"""tuple: The columns of an `effect_ledger` row, in table order."""

SUMMARY_SOURCE = (0, 0)
"""tuple: The `(source_type, source_id)` of a (body, effect) summary row."""


def split_sources(contributions, floor=None):
    """
    The shares of the sources of one effect on one body.

    Args:
        contributions (iterable): `(source, value)` pairs, value the contribution in the effect's unit (finite, zero or more); a source whose value is zero adds nothing.
        floor (float, optional): Share under which a source is not its own row; default `tuning.EFFECT_LEDGER_FLOOR`.

    Returns:
        dict: `n_sources` (those with a value over zero), `largest` (the source with the largest value, the first of equals), `largest_value`,
        `stored` (`[(source, share)]` in the order given, the largest always in it with share 1.0), `n_below` (how many fell under the floor) and
        `rest_share` (the sum of their shares), or `None` when no source has a value over zero.
    """
    floor = tuning.EFFECT_LEDGER_FLOOR if floor is None else floor
    pairs = []
    for source, value in contributions:
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"a contribution must be finite and not negative, got {value!r}")
        if value > 0.0:
            pairs.append((source, value))
    if not pairs:
        return None
    top = 0
    for index, (_, value) in enumerate(pairs):
        if value > pairs[top][1]:
            top = index
    largest, largest_value = pairs[top]
    stored, n_below, rest = [], 0, 0.0
    for index, (source, value) in enumerate(pairs):
        share = value / largest_value
        if index == top:
            stored.append((source, 1.0))
        elif share >= floor:
            stored.append((source, share))
        else:
            n_below += 1
            rest += share
    return {"n_sources": len(pairs), "largest": largest, "largest_value": largest_value, "stored": stored, "n_below": n_below, "rest_share": rest}


class EffectLedger:
    """
    The contributions a system's effect code computed, by effect and receiving body.

    `add(effect, receiver, source, value)` notes one contribution; `rows(star_system_id, ids)` turns them into the rows of the `effect_ledger` table.
    A receiver or source is the generator's own object; `ids` maps `id(object)` to its `(table, row id)`. An object with no stored row (a body the
    save does not file, such as a pseudo-source) is skipped, and `missing` counts how many were.
    """

    def __init__(self):
        self._by_effect = {}
        self.missing = 0

    def add(self, effect, receiver, source, value):
        """Notes that `source` contributes `value` (the effect's unit) to `receiver`; the same pair noted twice adds up."""
        if effect not in EFFECT_CODES:
            raise KeyError(f"unknown effect {effect!r}")
        if receiver is source:
            return
        per_receiver = self._by_effect.setdefault(effect, {})
        sources = per_receiver.setdefault(id(receiver), (receiver, {}))[1]
        entry = sources.setdefault(id(source), [source, 0.0])
        entry[1] += value

    def effects(self):
        """The names of the effects that have any contribution, in code order."""
        return [name for name, *_rest in EFFECTS if name in self._by_effect]

    def contributions(self, effect):
        """`[(receiver, [(source, value)])]` for one effect, in the order receivers were first noted."""
        return [(receiver, [(source, value) for source, value in sources.values()]) for receiver, sources in self._by_effect.get(effect, {}).values()]

    def __bool__(self):
        return bool(self._by_effect)

    def rows(self, star_system_id, ids, floor=None):
        """
        The `effect_ledger` rows of the contributions, as dicts of the table's columns.

        Args:
            star_system_id (int): The system the receivers are in.
            ids (dict): `id(object)` to `(table, row id)` for the system's inserted bodies.
            floor (float, optional): As in `split_sources`, for every effect; left out, `tuning.EFFECT_LEDGER_FLOOR` for light and UV and "largest only" for the
                on-demand effects.

        Returns:
            list: Nothing for a (body, effect) of light or UV that has one source. Otherwise one row per stored source of each (body, effect) and one summary row for it (`source_type` and `source_id` 0, `share` the sum of the shares under the floor,
            `n_below` and `largest_value` set). A (body, effect) whose receiver has no stored row is left out.
        """
        out = []
        for effect, code, *_rest in EFFECTS:
            effect_floor = math.inf if floor is None and effect in ON_DEMAND_EFFECTS else floor
            for receiver, pairs in self.contributions(effect):
                where = ids.get(id(receiver))
                if where is None:
                    self.missing += 1
                    continue
                known = []
                for source, value in pairs:
                    there = ids.get(id(source))
                    if there is None:
                        self.missing += 1
                    else:
                        known.append((there, value))
                split = split_sources(known, effect_floor)
                if split is None:
                    continue
                if split["n_sources"] == 1 and effect not in ON_DEMAND_EFFECTS:
                    continue      # a lone source is implied: it is the body's host star and the flux is on the body, so a lone-star system stores no light or UV row (Lane 13)
                body = {"star_system_id": star_system_id, "body_type": BODY_CODES[where[0]], "body_id": where[1], "effect": code}
                for (table, row_id), share in split["stored"]:
                    out.append({**body, "source_type": BODY_CODES[table], "source_id": row_id, "share": share, "n_below": None, "largest_value": None})
                out.append({**body, "source_type": SUMMARY_SOURCE[0], "source_id": SUMMARY_SOURCE[1], "share": split["rest_share"], "n_below": split["n_below"],
                            "largest_value": split["largest_value"]})
        return out
