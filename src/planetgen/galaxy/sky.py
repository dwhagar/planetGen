# planetgen/galaxy/sky.py

"""
The sky catalogue of a system (VIEW.23, S5 of reports/system-sky-from-the-system-2026-11-04.md, goal #4907): the stars an observer sees
from its place in the galaxy, from the STORED stars only, each with its direction, apparent magnitude, colour and light-time position.
This is the light version of VIEW.12: it reads what is stored (charted stars, `bright_stars`, scattered phenomena, cluster members) and
says which tiers it reached; VIEW.12 extends it with the built tiers in 9.1.

`sky_from_stars(observer, stars, at_unix, limit)` is the whole of it, on arrays; `sky_from(observer, source, at_unix, limit)` asks a
`source` (the sphere query of VIEW.22) for the stars in shells of distance, each shell above the luminosity floor a star must have to reach
the limit there, and hands them to `sky_from_stars`.

Per star, in this order:

1. **Its velocity.** A star with no stored velocity (a scatter star that is not built yet) gets the stand-in `unbuilt_velocity_kms`: the
   circular speed of the galaxy at its place, turned counterclockwise about the galactic pole, plus a velocity dispersion drawn from the star's own
   key (`draw.Stream` seeded by the galaxy key and the star's uid, so the same star always gets the same velocity; GEN.339's
   dispersion should key on the scatter seed so the fill agrees later). Stand-in SKC1.
   A star with no stored epoch (an unbuilt bright star) has its position taken to hold at the sky's time (stand-in SKC2), so it is not moved, only seen back by its light time.
2. **Where it appears.** `light_travel.retarded_offsets` on arrays: the star is moved from its own epoch to the sky's time on a straight line,
   and seen where it was when its light left it. The observer's own epoch is brought to each star's the same way.
3. **Whether it was there.** `light_travel.seen_age`: a star whose light left before it formed is left out (counted in the report). A star whose
   phase ended inside the light time is seen in the phase before it; its stored luminosity and temperature are kept (the stellar track at the
   earlier age is VIEW.12's) and the report counts them (`other_phase`). A star with no stored age is never left out.
4. **How bright.** `sky_photometry.apparent_magnitude` at the retarded distance (the light's travel, `c tau`), through the dust of the disk along
   the line of sight (`galactic_extinction_av`, worked only for the stars that pass without it, since dust only dims). A star fainter than
   `limit` is dropped. Nothing stored for a star's luminosity or temperature leaves it out (counted: `no_photometry`).
5. **What colour.** The blackbody colour of its temperature (`photometry.star_colour`, linear sRGB with the brightest channel 1), from a
   table of 512 temperatures between 1,500 K and 60,000 K, linear between rows (the table is within 0.01 of a direct calculation). Reddening by
   dust and the desaturation of faint stars are VIEW.11's (9.1).
6. **The observer's own system** is left out (every star with the observer's `system_id`, so a companion goes too).

The directions are unit vectors in the galaxy frame (z up); turning them into the system's own frame is the sky route's.
"""

import functools
import math
from typing import NamedTuple

import numpy as np

from planetgen.db import sky_sphere
from planetgen.galaxy import potential as galacticPotential
from planetgen.physics import light_travel, photometry, sky_photometry
from planetgen.physics import unit_registry as units
from planetgen.util import draw

DEFAULT_LIMIT = 7.5
"""float: The limiting magnitude of the drawn stars when the caller gives none (the setting's default, Boss 5:29 PM EST on 2026-11-04)."""

DISPERSION_KMS = (35.0, 25.0, 20.0)
"""tuple: The velocity dispersion `(radial, tangential, vertical)`, km/s, of a star with no stored velocity (stand-in SKC1)."""

COLOUR_TABLE_K = (1500.0, 60000.0)
COLOUR_TABLE_ROWS = 512
"""The temperature range of the colour table (log spaced) and its rows."""

SHELL_RATIO = 2.0
"""float: Each shell `sky_from` asks the source for reaches this many times as far as the one inside."""

FIRST_SHELL_PC = 8.0
"""float: The first shell's outer radius, parsecs."""

MAX_DISTANCE_PC = 40000.0
"""float: The farthest a star is asked for, parsecs."""

DUST_CHUNK = 4096
"""int: Stars whose dust is integrated at once."""

_pc_to_m = units.converter("pc", "m")
_m_to_pc = units.converter("m", "pc")
_km_to_m = units.converter("km", "m")
_yr_to_s = units.converter("yr", "s")


class SkyStars(NamedTuple):
    """
    Stored stars as arrays, what the sphere query of VIEW.22 returns and `sky_from_stars` takes. Every array has one row a star.

    Attributes:
        uid: `(n,)` the star's id (anything that prints the same every time; it seeds the stand-in velocity).
        tier: `(n,)` the layer it came from (`"charted"`, `"bright_star"`, `"phenomenon"`, `"cluster"`, ...).
        system_id: `(n,)` the id of the system it belongs to, `-1` for none (a bright star or a phenomenon has none).
        position_pc: `(n, 3)` galaxy-frame parsecs at `epoch_unix`.
        epoch_unix: `(n,)` seconds, the time the position holds at.
        velocity_kms: `(n, 3)` km/s, `nan` where nothing is stored.
        luminosity_sol: `(n,)` bolometric luminosity, solar (`nan` for none).
        temperature_k: `(n,)` effective temperature (`nan` for none).
        age_years: `(n,)` the star's age at the sky's time (`nan` for none stored).
        phase_end_age_years: `(n,)` the age at which it leaves its present phase (`nan` for none stored).
    """

    uid: np.ndarray
    tier: np.ndarray
    system_id: np.ndarray
    position_pc: np.ndarray
    epoch_unix: np.ndarray
    velocity_kms: np.ndarray
    luminosity_sol: np.ndarray
    temperature_k: np.ndarray
    age_years: np.ndarray
    phase_end_age_years: np.ndarray


class Observer(NamedTuple):
    """
    Where the sky is seen from: `position_pc` (galaxy frame, parsecs) at `epoch_unix`, `velocity_kms` (km/s), and the id of its own system
    (`system_id`, `-1` for none), whose stars are left out.
    """

    position_pc: tuple
    epoch_unix: float
    velocity_kms: tuple = (0.0, 0.0, 0.0)
    system_id: int = -1


class Sky(NamedTuple):
    """
    The sky: one row a visible star. `direction` is a unit vector `(n, 3)` in the galaxy frame, `distance_pc` the length of the apparent offset
    (the light's travel), `light_time_years` the time the light took, `magnitude` the apparent V magnitude, `colour_linear` the linear sRGB `(n, 3)`
    with the brightest channel 1, `uid` and `tier` as in `SkyStars`; `report` is a dict (see `sky_from_stars`).
    """

    uid: np.ndarray
    tier: np.ndarray
    direction: np.ndarray
    distance_pc: np.ndarray
    light_time_years: np.ndarray
    magnitude: np.ndarray
    colour_linear: np.ndarray
    report: dict
    offset_pc: np.ndarray = None
    velocity_kms: np.ndarray = None
    temperature_k: np.ndarray = None
    dimmed_absolute_magnitude: np.ndarray = None


def empty_stars():
    """A `SkyStars` with no rows."""
    return SkyStars(np.empty(0, dtype=object), np.empty(0, dtype=object), np.empty(0, dtype=int), np.empty((0, 3)),
                    np.empty(0), np.empty((0, 3)), np.empty(0), np.empty(0), np.empty(0), np.empty(0))


def concatenate(parts):
    """The rows of several `SkyStars` in one."""
    parts = [p for p in parts if len(p.uid)]
    if not parts:
        return empty_stars()
    return SkyStars(*[np.concatenate([getattr(p, name) for p in parts]) for name in SkyStars._fields])


def luminosity_floor_solar(distance_pc, limit=DEFAULT_LIMIT):
    """
    The bolometric luminosity, in suns, a star at `distance_pc` needs to reach apparent magnitude `limit` at the best: no dust and a bolometric
    correction of 0 (no star is brighter in V for its luminosity). Below it nothing is visible there, so it is the floor a sphere query uses.
    """
    return 10.0 ** (0.4 * (photometry.SUN_ABSOLUTE_BOLOMETRIC_MAG + 5.0 * np.log10(np.asarray(distance_pc, dtype=float) / 10.0) - limit))


def shell_radii_pc(max_distance_pc=MAX_DISTANCE_PC):
    """The outer radii of the shells `sky_from` asks the source for: `FIRST_SHELL_PC`, then each `SHELL_RATIO` times the last, the last at the maximum."""
    radii = [FIRST_SHELL_PC]
    while radii[-1] < max_distance_pc:
        radii.append(min(radii[-1] * SHELL_RATIO, max_distance_pc))
    return radii


def unbuilt_velocity_kms(position_pc, uid, galaxy_key, circular_speed_kms):
    """
    The stand-in velocity of stars with no stored one (SKC1): the galaxy's circular speed at each star's cylindrical radius, along the
    counterclockwise tangent (the sense `galaxy.frame` calls increasing longitude), plus a dispersion `DISPERSION_KMS` in the radial, tangential and
    vertical directions drawn from the star's own seed, so the same star always gets the same velocity.

    Args:
        position_pc: `(n, 3)` galaxy-frame parsecs.
        uid: `(n,)` star ids.
        galaxy_key (str): Names the galaxy (its id, seed and epoch).
        circular_speed_kms: A function of the radius in parsecs (an array) giving the circular speed in km/s.

    Returns:
        numpy.ndarray: `(n, 3)` km/s.
    """
    p = np.asarray(position_pc, dtype=float).reshape(-1, 3)
    radius = np.hypot(p[:, 0], p[:, 1])
    safe = np.where(radius > 0.0, radius, 1.0)
    radial = np.stack([p[:, 0] / safe, p[:, 1] / safe, np.zeros(len(p))], axis=1)
    radial[radius == 0.0] = (1.0, 0.0, 0.0)
    tangent = np.stack([-radial[:, 1], radial[:, 0], np.zeros(len(p))], axis=1)
    speed = np.asarray(circular_speed_kms(radius), dtype=float)
    draws = np.empty((len(p), 3))
    for i, key in enumerate(uid):
        stream = draw.Stream(f"sky-velocity|{galaxy_key}|{key}")
        draws[i] = (stream.gauss(0.0, 1.0), stream.gauss(0.0, 1.0), stream.gauss(0.0, 1.0))
    sr, st, sz = DISPERSION_KMS
    up = np.array([0.0, 0.0, 1.0])
    return (tangent * (speed + st * draws[:, 1])[:, None] + radial * (sr * draws[:, 0])[:, None] + up * (sz * draws[:, 2])[:, None])


@functools.lru_cache(maxsize=1)
def _colour_table():
    temperatures = np.geomspace(COLOUR_TABLE_K[0], COLOUR_TABLE_K[1], COLOUR_TABLE_ROWS)
    rows = np.array([photometry.star_colour(float(t)).linear for t in temperatures])
    return np.log(temperatures), rows


def star_colours(temperature_k):
    """The linear sRGB `(n, 3)` of blackbody stars of `temperature_k` (brightest channel 1), interpolated in log temperature from a table of
    `COLOUR_TABLE_ROWS` rows between `COLOUR_TABLE_K`; a temperature outside takes the end row."""
    log_t, rows = _colour_table()
    x = np.log(np.clip(np.asarray(temperature_k, dtype=float), COLOUR_TABLE_K[0], COLOUR_TABLE_K[1]))
    return np.stack([np.interp(x, log_t, rows[:, k]) for k in range(3)], axis=1)


def extinction_av(observer_pc, direction, distance_pc, dust):
    """The dust `A_V` from `observer_pc` to each of the points `distance_pc` away along `direction` (`(n, 3)` unit vectors),
    `DUST_CHUNK` at a time (`sky_photometry.galactic_extinction_av` with the `dust` overrides)."""
    out = np.empty(len(distance_pc))
    for start in range(0, len(distance_pc), DUST_CHUNK):
        sl = slice(start, start + DUST_CHUNK)
        out[sl] = sky_photometry.galactic_extinction_av(observer_pc, direction[sl], distance_pc[sl], **dust)
    return out


def sky_from_stars(observer, stars, at_unix, limit=DEFAULT_LIMIT, galaxy_key="", circular_speed_kms=None, dust=None):
    """
    The sky of `observer` at the time `at_unix` from the stored `stars` (a `SkyStars`), the stars brighter than `limit`, as the module doc says.

    Args:
        observer (Observer): Where the sky is seen from.
        stars (SkyStars): The stored stars (any number; the caller's sphere query bounds them).
        at_unix (float): The sky's time, seconds.
        limit (float): The faintest star kept (apparent V magnitude).
        galaxy_key (str): Names the galaxy for the stand-in velocities.
        circular_speed_kms: A function of radius in parsecs for the stand-in velocities; defaults to the running galaxy's potential
            (`galaxy.potential.current`).
        dust (dict): Overrides of the dust law (`sky_photometry.dust_av_per_pc` keywords; `profile` is a `physics.dust_model.DustProfile`). `sky_from_database` gives a galaxy the dust of its own gas when this is `None`; the other functions use the fixed Milky Way law.

    Returns:
        Sky: The visible stars, brightest first, with `report`: `{"limit", "stars_in", "tiers": {tier: {"in", "visible", "faintest_magnitude",
            "farthest_pc"}}, "excluded_own_system", "no_photometry", "unknown_epoch", "not_formed_yet", "other_phase", "unbuilt_velocity", "visible"}`.
    """
    dust = dict(dust or {})
    n = len(stars.uid)
    tiers = sorted({str(t) for t in stars.tier})
    report = {"limit": float(limit), "stars_in": n, "tiers": {t: {"in": int(np.sum(stars.tier == t)), "visible": 0, "faintest_magnitude": None,
                                                                   "farthest_pc": None} for t in tiers},
              "excluded_own_system": 0, "no_photometry": 0, "unknown_epoch": 0, "not_formed_yet": 0, "other_phase": 0, "unbuilt_velocity": 0, "visible": 0}
    nothing = Sky(np.empty(0, dtype=object), np.empty(0, dtype=object), np.empty((0, 3)), np.empty(0), np.empty(0), np.empty(0),
                  np.empty((0, 3)), report, np.empty((0, 3)), np.empty((0, 3)), np.empty(0), np.empty(0))
    if n == 0:
        return nothing

    keep = np.ones(n, dtype=bool)
    if observer.system_id >= 0:
        own = np.asarray(stars.system_id) == observer.system_id
        report["excluded_own_system"] = int(own.sum())
        keep &= ~own
    photometry_known = np.isfinite(stars.luminosity_sol) & (stars.luminosity_sol > 0) & np.isfinite(stars.temperature_k) & (stars.temperature_k > 0)
    report["no_photometry"] = int((keep & ~photometry_known).sum())
    keep &= photometry_known
    index = np.nonzero(keep)[0]
    if len(index) == 0:
        return nothing

    def column(name):
        return np.asarray(getattr(stars, name), dtype=float)[index]
    position = column("position_pc")
    epoch = column("epoch_unix").copy()
    no_epoch = ~np.isfinite(epoch)
    report["unknown_epoch"] = int(no_epoch.sum())
    epoch[no_epoch] = at_unix
    velocity = column("velocity_kms").copy()
    unknown = ~np.all(np.isfinite(velocity), axis=1)
    report["unbuilt_velocity"] = int(unknown.sum())
    if unknown.any():
        velocity[unknown] = unbuilt_velocity_kms(position[unknown], np.asarray(stars.uid)[index][unknown], galaxy_key,
                                                 circular_speed_function(circular_speed_kms))
    obs_position = np.asarray(observer.position_pc, dtype=float).reshape(3)
    offsets, light_time_s = apparent_offsets(observer, position, velocity, epoch, at_unix)
    light_time_years = light_time_s / _yr_to_s(1.0)

    age, phase_end = column("age_years"), column("phase_end_age_years")
    seen, formed, other_phase = light_travel.seen_age(age * _yr_to_s(1.0), light_time_s, phase_end * _yr_to_s(1.0))
    formed = formed | ~np.isfinite(age)
    other_phase &= np.isfinite(age)
    report["not_formed_yet"] = int((~formed).sum())
    report["other_phase"] = int((other_phase & formed).sum())

    distance = _m_to_pc(np.sqrt(np.sum(offsets * offsets, axis=1)))
    usable = formed & (distance > 0.0)
    luminosity, temperature = column("luminosity_sol"), column("temperature_k")
    clear = np.full(len(index), np.inf)
    clear[usable] = sky_photometry.apparent_magnitude(luminosity[usable], temperature[usable], distance[usable])
    candidate = np.nonzero(clear <= limit)[0]                      # dust only dims, so a star fainter than the limit without it stays out
    direction = offsets[candidate] / (distance[candidate, None] * _pc_to_m(1.0))
    a_v = extinction_av(obs_position, direction, distance[candidate], dust) if len(candidate) else np.empty(0)
    magnitude = clear[candidate] + a_v
    bright = magnitude <= limit
    chosen = candidate[bright]
    magnitude = magnitude[bright]
    order = np.argsort(magnitude, kind="stable")
    chosen, magnitude = chosen[order], magnitude[order]
    rows = index[chosen]

    result_tier = np.asarray(stars.tier, dtype=object)[rows]
    for tier in tiers:
        mask = result_tier == tier
        if mask.any():
            report["tiers"][tier].update(visible=int(mask.sum()), faintest_magnitude=float(magnitude[mask].max()),
                                         farthest_pc=float(distance[chosen][mask].max()))
    report["visible"] = int(len(chosen))
    advanced = _m_to_pc(offsets[chosen] + _km_to_m(velocity[chosen]) * light_time_s[chosen][:, None])
    return Sky(np.asarray(stars.uid, dtype=object)[rows], result_tier, direction[bright][order], distance[chosen],
               light_time_years[chosen], magnitude, star_colours(temperature[chosen]), report, advanced, velocity[chosen],
               temperature[chosen], magnitude - 5.0 * np.log10(distance[chosen] / 10.0))


def _km_to_pc_per_s(velocity_kms):
    """A velocity in km/s as parsecs a second."""
    return _m_to_pc(_km_to_m(np.asarray(velocity_kms, dtype=float)))
def circular_speed_function(circular_speed_kms=None):
    """`circular_speed_kms` (a function of the radius in parsecs, arrays, giving km/s), or the running galaxy's potential's
    (`galaxy.potential.current`) when `None`: what the stand-in velocities turn with."""
    if circular_speed_kms is not None:
        return circular_speed_kms
    pot = galacticPotential.current()
    return lambda radius: galacticPotential.circular_speed_kms(pot, radius)
def apparent_offsets(observer, position_pc, velocity_kms, epoch_unix, at_unix):
    """`(offsets_m, light_time_s)` of `light_travel.retarded_offsets` for objects at `position_pc` (`(n, 3)`) moving at
    `velocity_kms` (`(n, 3)`) at their own `epoch_unix` (`(n,)`), seen by `observer` (a `sky.Observer`) at `at_unix`: the
    observer is first brought to each object's own epoch, so that object(t) - observer(t) is what the light time is worked
    out from."""
    obs_position = np.asarray(observer.position_pc, dtype=float).reshape(3)
    obs_velocity = np.asarray(observer.velocity_kms, dtype=float).reshape(3)
    obs_at_object = obs_position[None, :] + _km_to_pc_per_s(obs_velocity)[None, :] * (epoch_unix - float(observer.epoch_unix))[:, None]
    return light_travel.retarded_offsets(_pc_to_m(position_pc - obs_at_object), _km_to_m(velocity_kms), _km_to_m(obs_velocity),
                                         at_unix - epoch_unix)


def sky_from(observer, source, at_unix, limit=DEFAULT_LIMIT, galaxy_key="", circular_speed_kms=None, dust=None, max_distance_pc=MAX_DISTANCE_PC):
    """
    The sky of `observer`: asks `source` for the stored stars in shells of distance and works the sky out of them (`sky_from_stars`).

    `source(center_pc, inner_pc, outer_pc, floor_solar)` returns a `SkyStars` of the stored stars between `inner_pc` and `outer_pc` of `center_pc` with a
    luminosity of at least `floor_solar` suns (the sphere query of VIEW.22). The floor of a shell is `luminosity_floor_solar` at its inner edge, the
    smallest a star there may have and still reach `limit`. The report gains `"shells"`: the radii asked for.
    """
    radii = shell_radii_pc(max_distance_pc)
    inner, parts = 0.0, []
    for outer in radii:
        parts.append(source(tuple(observer.position_pc), inner, outer, float(luminosity_floor_solar(max(inner, 0.1), limit))))
        inner = outer
    sky = sky_from_stars(observer, concatenate(parts), at_unix, limit, galaxy_key, circular_speed_kms, dust)
    sky.report["shells"] = radii
    return sky


SOURCE_TIERS = {0: "charted", 1: "bright_star", 2: "phenomenon"}
"""dict: The tier name of each `source` number of `db.sky_sphere.SphereStars` (a star with a `cluster_id` is the tier `"cluster"`)."""

GY_YEARS = 1.0e9
"""float: Years in a gigayear."""


def stars_from_sphere(sphere):
    """
    A `SkyStars` from the stars of a sphere query (`db.sky_sphere.SphereStars`, VIEW.22; read by its attribute names, so this module need not import
    the database layer). The tier comes from the source and the cluster: charted, bright_star, phenomenon, or cluster for a star with a `cluster_id`.
    The uid is `"<source>:<ident>"`. Ages and lifespans become years; a system id is kept for the charted stars only (`-1` for the rest).
    """
    source = np.asarray(sphere.source)
    tier = np.array([("cluster" if c and s != 2 else SOURCE_TIERS[int(s)]) for s, c in zip(source, np.asarray(sphere.cluster_id))], dtype=object)
    uid = np.array([f"{int(s)}:{int(i)}" for s, i in zip(source, np.asarray(sphere.ident))], dtype=object)
    return SkyStars(
        uid, tier, np.where(source == 0, np.asarray(sphere.system_id), -1).astype(int),
        np.column_stack((sphere.x_pc, sphere.y_pc, sphere.z_pc)).astype(float).reshape(-1, 3), np.asarray(sphere.epoch_unix, dtype=float),
        np.column_stack((sphere.vx_kms, sphere.vy_kms, sphere.vz_kms)).astype(float).reshape(-1, 3), np.asarray(sphere.luminosity_sol, dtype=float),
        np.asarray(sphere.temperature_k, dtype=float), np.asarray(sphere.age_gy, dtype=float) * GY_YEARS,
        np.asarray(sphere.lifespan_gy, dtype=float) * GY_YEARS)


def sphere_source(query, levels=None):
    """
    The `source` of `sky_from` for a sphere query `query(center, radius_pc, min_luminosity_sol)` (`db.sky_sphere.stars_in_sphere` with its connection
    bound): a shell is the sphere of its outer radius without the rows nearer than its inner one. When `levels` is a list, each shell appends
    `(inner_pc, outer_pc, floor_solar, complete_above_sol)` to it, `complete_above_sol` being the luminosity above which the cells the sphere reached
    are complete (`CellLevels.complete_above_sol`; `inf` when the query gives none).
    """
    def source(center_pc, inner_pc, outer_pc, floor_solar):
        sphere = query(center_pc, outer_pc, floor_solar)
        if levels is not None:
            reached = getattr(sphere, "levels", None)
            levels.append((float(inner_pc), float(outer_pc), float(floor_solar), float(getattr(reached, "complete_above_sol", math.inf))))
        stars = stars_from_sphere(sphere)
        keep = np.asarray(sphere.distance_pc, dtype=float) >= inner_pc
        return SkyStars(*[getattr(stars, name)[keep] for name in SkyStars._fields])
    return source


def built_to_pc(levels):
    """
    How far the stored stars are complete to the limit, parsecs: the outer radius of the last of the shells in `levels` (from `sphere_source`) that
    are, running outwards without a gap, ones whose cells are drawn down to the shell's luminosity floor (`complete_above_sol` at or under the floor).
    0.0 when even the first shell is not.
    """
    reached = 0.0
    for inner, outer, floor, complete in sorted(levels):
        if complete > floor:
            break
        reached = outer
    return reached


def default_dust(conn):
    """
    The dust a galaxy's sky is seen through when the caller gives none (VIEW.44): the dust of its own gas (`physics.dust_model.dust_for_shape`) for the stored shape;
    the fixed Milky Way law (an empty override) for a database with no shape.
    """
    from planetgen.db import store
    from planetgen.physics import dust_model
    info = store.get_galaxy_shape(conn, light=True)
    return {} if info is None else dust_model.dust_for_shape(info.shape)


def sky_from_database(conn, observer, at_unix, limit=DEFAULT_LIMIT, galaxy_key="", circular_speed_kms=None, dust=None,
                      max_distance_pc=MAX_DISTANCE_PC):
    """
    The sky of `observer` from the stars the database stores (`sky_from` over the sphere query of VIEW.22, `db.sky_sphere.stars_in_sphere`, on `conn`).
    The report gains `"built_to_pc"` (`built_to_pc`: how far the stored stars are complete to `limit`) and `"complete_above_sol"` (the largest level of
    the shells' cells), the "built to magnitude X within Y pc" of the plan.
    """
    if dust is None:
        dust = default_dust(conn)
    levels = []
    sky = sky_from(observer, sphere_source(lambda center, radius, floor: sky_sphere.stars_in_sphere(conn, center, radius, floor), levels), at_unix, limit,
                   galaxy_key, circular_speed_kms, dust, max_distance_pc)
    sky.report["built_to_pc"] = built_to_pc(levels)
    sky.report["complete_above_sol"] = max((level[3] for level in levels), default=math.inf)
    return sky
