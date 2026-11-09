# planetgen/db/galaxy_counts.py

"""
Running galaxy totals and stored facet counts (DB.121, audit C9, F7 and F9).

A page that counts the galaxy used to count it again whenever a sector was added: `COUNT(*)` over millions of systems takes
about 3 s per million, and the stored answer (`countcache`) is made again each time the content stamp moves. The counts are
now kept as they change:

- `sector_facets` holds, for each generated sector, how many of each thing it holds: a `total` row for each of `TOTALS`
  (systems, stars, planets, moons, belts, habitable worlds) and a row for each value of each facet the lists and the search
  page offer (a star's spectral and Yerkes class, a planet's or moon's class, body type and life class, a system's binary
  state and octant). It is made in the sector's own save (`store.record_sector_stats`) from the sector's rows, made again
  when the sector is edited (`store.refresh_edited_sectors`) and dropped when its cell is emptied or deleted
  (`store.forget_sector_fill`), always in the same transaction as the change.
- `galaxy_totals` is one row, the sum of every sector's `total` rows plus the mass and luminosity of `sector_stats`, moved
  by the difference between a sector's old rows and its new ones (`refresh`, `forget`). The row is locked last in a
  transaction, like `galaxy_shape`'s, so two workers saving at once wait only for the commit.

What no sector's rows hold is counted when it is asked (`outside`): the systems filed in no sector or in a sector with no rows, their
bodies, and the moons of free planets. They are few in a built galaxy; a stored count of them would be one more thing to keep true.
A read is `stored + outside`, one row (or a sum over the sectors' facet rows) and one small count, for any size of galaxy.

A galaxy made before this table, or a cell filled by hand, has no rows for a sector: the systems in a sector with no rows are counted
from the tables with the rest (`outside`), so a read is right whatever the rows hold, and `rebuild` (the migration, `planetgen
fill-block-totals`) only makes the reads fast. `check-counts` compares the rows and the totals with a full count. Nothing here changes
a generated value: the counts are read from the rows a sector holds.
"""

TOTALS = ("systems", "stars", "planets", "moons", "belts", "habitable")
"""tuple: The `value`s of a sector's `total` facet; `galaxy_totals` has a column for each (`habitable` is `habitable_worlds`)."""

FACETS = ("binary", "octant", "spectral", "yerkes", "class", "body", "life", "moon_class", "moon_body", "moon_life")
"""tuple: The facets a sector keeps beside `total`. Their names are the search page's facets (`query.SEARCH_TAG_FACETS`) where
there is one."""

_HELD = ("SELECT s.id FROM sectors s JOIN sector_facets f ON f.shell_index = s.shell_index AND f.band_index = s.band_index"
         " AND f.wedge_index = s.wedge_index AND f.facet = 'total' AND f.value = 'systems'")
"""str: The sectors that have rows: what `sector_facets` holds. Whatever lies in any other sector is `outside`, so a sector with no rows
(a cell filled by hand, a galaxy from before the table, a sector in the middle of a rebuild) is counted from the tables and a read is
right without them."""

_OUTSIDE_WHERE = f"sector_id IS NULL OR sector_id NOT IN ({_HELD})"
_OUTSIDE_SYSTEMS = f"SELECT id FROM star_systems WHERE {_OUTSIDE_WHERE}"


def _in_scope(column, kind, sector_id, free=True):
    """
    The condition that `column` (a body's `star_system_id`) lies in the scope, with its parameters. `outside` also takes a row
    with no system (`free`: the moon of a free planet; a column that is never NULL passes `free=False`).
    """
    if kind == "sector":
        return f"{column} IN (SELECT id FROM star_systems WHERE sector_id = ?)", (sector_id,)
    clause = f"{column} IN ({_OUTSIDE_SYSTEMS})"
    return (f"({column} IS NULL OR {clause})" if free else clause), ()


def counts_in(conn, kind, sector_id=None):
    """
    What the rows of one scope hold, as `{(facet, value): n}`: `kind` is `"sector"` (the sector `sector_id`) or `"outside"`
    (`outside`). Counted from the tables; the `habitable` total is left out of `"outside"`.
    """
    from planetgen.db import query

    found = {}

    def add(facet, value, n):
        if value is not None and n:
            found[(facet, str(value))] = found.get((facet, str(value)), 0) + int(n)

    if kind == "sector":
        where, params = "sector_id = ?", (sector_id,)
    else:
        where, params = f"({_OUTSIDE_WHERE})", ()
    systems = 0
    for row in conn.execute(f"SELECT is_binary AS b, quadrant AS q, COUNT(*) AS n FROM star_systems WHERE {where} GROUP BY b, q", params).fetchall():
        systems += row["n"]
        add("binary", "yes" if row["b"] else "no", row["n"])
        add("octant", row["q"], row["n"])
    stars = 0
    condition, params = _in_scope("st.star_system_id", kind, sector_id, free=False)
    for row in conn.execute(f"SELECT SUBSTR(st.star_type, 1, 1) AS v, st.yerkes_class AS y, COUNT(*) AS n FROM stars st WHERE {condition} GROUP BY v, y", params).fetchall():
        stars += row["n"]
        add("spectral", row["v"], row["n"])
        add("yerkes", row["y"], row["n"])
    totals = {"systems": systems, "stars": stars}
    for table, prefix in (("planets", ""), ("moons", "moon_")):
        condition, params = _in_scope("b.star_system_id", kind, sector_id)
        count = 0
        for row in conn.execute(f"SELECT b.planet_class AS c, b.body_type AS t, COUNT(*) AS n FROM {table} b WHERE {condition} GROUP BY c, t", params).fetchall():
            count += row["n"]
            add(f"{prefix}class", row["c"], row["n"])
            add(f"{prefix}body", row["t"], row["n"])
        totals[table] = count
        for row in query.life_class_counts(conn, table, condition, params):
            add(f"{prefix}life", row["v"], row["c"])
    condition, params = _in_scope("star_system_id", kind, sector_id, free=False)
    totals["belts"] = conn.execute(f"SELECT COUNT(*) AS n FROM asteroid_belts WHERE {condition}", params).fetchone()["n"]
    if kind == "sector":
        from planetgen.db import store
        totals["habitable"] = store.sector_habitable_worlds(conn, sector_id)
    for name in TOTALS:
        if kind == "sector" or name != "habitable":
            found[("total", name)] = int(totals.get(name, 0))
    return found


def _stored(conn, address):
    return {(row["facet"], row["value"]): int(row["n"]) for row in conn.execute(
        "SELECT facet, value, n FROM sector_facets WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", tuple(address)).fetchall()}


def _stats_figures(conn, address):
    """`(mass, luminosity)` of the sector's stats row as it is now (0 where there is none or it is NULL)."""
    row = conn.execute("SELECT total_mass_kg, total_luminosity_sol FROM sector_stats WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                       tuple(address)).fetchone()
    return (0.0, 0.0) if row is None else (float(row["total_mass_kg"] or 0.0), float(row["total_luminosity_sol"] or 0.0))


def _apply(conn, delta, mass, luminosity):
    """Moves the totals row by the differences, making it when it is not there. Last in the transaction (see the module)."""
    if not any(delta.values()) and not mass and not luminosity:
        return
    conn.execute(
        "INSERT INTO galaxy_totals (id, systems, stars, planets, moons, belts, habitable_worlds, mass_kg, luminosity_sol)"
        " VALUES (1, ?, ?, ?, ?, ?, ?, ?, ?) ON DUPLICATE KEY UPDATE systems = systems + ?, stars = stars + ?,"
        " planets = planets + ?, moons = moons + ?, belts = belts + ?, habitable_worlds = habitable_worlds + ?,"
        " mass_kg = mass_kg + ?, luminosity_sol = luminosity_sol + ?",
        (*[delta[name] for name in TOTALS], mass, luminosity, *[delta[name] for name in TOTALS], mass, luminosity))


def _replace(conn, address, rows):
    conn.execute("DELETE FROM sector_facets WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", tuple(address))
    if rows:
        conn.executemany("INSERT INTO sector_facets (shell_index, band_index, wedge_index, facet, value, n) VALUES (?, ?, ?, ?, ?, ?)",
                         [(*address, facet, value, n) for (facet, value), n in sorted(rows.items())])


def refresh(conn, sector_id, address, old_figures=None):
    """
    Makes the sector's rows again from what it holds and moves the totals by the difference. `old_figures` is the `(mass,
    luminosity)` of its stats row before the change that brings this call (`stats_figures`; `None` for a change that left them
    as they are); the new ones are read from the stats row as it stands, after the change and whatever the change brought with it.
    """
    old = _stored(conn, address)
    new = counts_in(conn, "sector", sector_id)
    delta = {name: new.get(("total", name), 0) - old.get(("total", name), 0) for name in TOTALS}
    figures = _stats_figures(conn, address)
    old_figures = figures if old_figures is None else old_figures
    _replace(conn, address, new)
    _apply(conn, delta, figures[0] - old_figures[0], figures[1] - old_figures[1])


def forget(conn, address, old_figures=None):
    """The sector's cell was emptied or deleted: takes its rows out of the totals. `old_figures` as in `refresh`, the stats row's
    `(mass, luminosity)` before the change (read now when `None`)."""
    old = _stored(conn, address)
    if old_figures is None:
        old_figures = _stats_figures(conn, address)
    if not old and not any(old_figures):
        return
    _replace(conn, address, {})
    _apply(conn, {name: -old.get(("total", name), 0) for name in TOTALS}, -old_figures[0], -old_figures[1])


def stats_figures(conn, address):
    """`_stats_figures`, for the hooks that read the figures before they write the stats row."""
    return _stats_figures(conn, address)


def totals(conn):
    """
    The galaxy's counts in full, as `{"systems", "stars", "planets", "moons", "belts"}`: the stored row plus `outside`.
    """
    row = conn.execute("SELECT systems, stars, planets, moons, belts FROM galaxy_totals WHERE id = 1").fetchone()
    stored = {name: int(row[name]) for name in TOTALS if name != "habitable"} if row else {name: 0 for name in TOTALS if name != "habitable"}
    rest = counts_in(conn, "outside")
    return {name: stored[name] + rest.get(("total", name), 0) for name in stored}


def facet(conn, name):
    """
    One facet's counts over the galaxy, as `{value: n}`: the sectors' rows summed plus `outside`.
    """
    found = {row["value"]: int(row["n"]) for row in conn.execute(
        "SELECT value, SUM(n) AS n FROM sector_facets WHERE facet = ? GROUP BY value", (name,)).fetchall()}
    for (key, value), n in counts_in(conn, "outside").items():
        if key == name:
            found[value] = found.get(value, 0) + n
    return found


def rebuild(conn, on_progress=None, batch=200):
    """
    Makes every row again from the tables (the migration, `planetgen fill-block-totals`, and the repair `check-counts` points to):
    the facet rows of every sector whose stats are filled, then the totals row as their sum plus the stats rows' mass and
    luminosity. Commits as it goes. Returns the number of sectors done.
    """
    conn.execute("DELETE FROM sector_facets")
    conn.execute("DELETE FROM galaxy_totals")
    sectors = conn.execute(
        "SELECT s.id AS id, st.shell_index AS shell, st.band_index AS band, st.wedge_index AS wedge FROM sector_stats st"
        " LEFT JOIN sectors s ON s.shell_index = st.shell_index AND s.band_index = st.band_index AND s.wedge_index = st.wedge_index"
        " WHERE st.bright_level_sol = 0 ORDER BY s.id").fetchall()
    # A filled sector whose row was lost still has stats and holds nothing: its rows are zeros.
    done = 0
    for row in sectors:
        _replace(conn, (row["shell"], row["band"], row["wedge"]), counts_in(conn, "sector", row["id"]))
        done += 1
        if done % batch == 0:
            conn.commit()
            if on_progress is not None:
                on_progress(done, len(sectors))
    sums = {row["value"]: int(row["n"]) for row in conn.execute(
        "SELECT value, SUM(n) AS n FROM sector_facets WHERE facet = 'total' GROUP BY value").fetchall()}
    stats = conn.execute("SELECT COALESCE(SUM(total_mass_kg), 0) AS m, COALESCE(SUM(total_luminosity_sol), 0) AS l"
                         " FROM sector_stats WHERE bright_level_sol = 0").fetchone()
    conn.execute(
        "INSERT INTO galaxy_totals (id, systems, stars, planets, moons, belts, habitable_worlds, mass_kg, luminosity_sol)"
        " VALUES (1, ?, ?, ?, ?, ?, ?, ?, ?)",
        (*[sums.get(name, 0) for name in TOTALS], float(stats["m"]), float(stats["l"])))
    conn.commit()
    return done


def problems(conn):
    """
    What differs between the stored counts and a full count, as a list of sentences (empty when they agree): the
    totals row against the sum of the sectors' rows, the mass and luminosity against the stats rows, and the whole galaxy's tables
    against the stored counts and `outside`.
    """
    found = []
    row = conn.execute("SELECT * FROM galaxy_totals WHERE id = 1").fetchone()
    sums = {r["value"]: int(r["n"]) for r in conn.execute("SELECT value, SUM(n) AS n FROM sector_facets WHERE facet = 'total' GROUP BY value").fetchall()}
    for name in TOTALS:
        column = "habitable_worlds" if name == "habitable" else name
        have = int(row[column]) if row else 0
        if have != sums.get(name, 0):
            found.append(f"The galaxy's stored {name} total is {have}, the sectors' rows add up to {sums.get(name, 0)}.")
    stats = conn.execute("SELECT COALESCE(SUM(total_mass_kg), 0) AS m, COALESCE(SUM(total_luminosity_sol), 0) AS l FROM sector_stats WHERE bright_level_sol = 0").fetchone()
    for label, column, want in (("mass", "mass_kg", stats["m"]), ("luminosity", "luminosity_sol", stats["l"])):
        have = float(row[column]) if row else 0.0
        if abs(have - float(want)) > 1e-9 * max(abs(float(want)), 1.0):
            found.append(f"The galaxy's stored {label} total is {have}, the sectors' stats add up to {float(want)}.")
    live = counts_in(conn, "outside")
    for table, name in (("star_systems", "systems"), ("stars", "stars"), ("planets", "planets"), ("moons", "moons"), ("asteroid_belts", "belts")):
        counted = conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
        have = (int(row[name]) if row else 0) + live.get(("total", name), 0)
        if int(counted) != have:
            found.append(f"The galaxy holds {counted} {name} but the stored counts say {have}.")
    return found
