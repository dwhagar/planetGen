# planetgen/api/system_sky.py

"""
The sky of a system, built for the route `GET /api/systems/<id>/sky` (VIEW.31, S13 of reports/system-sky-from-the-system-2026-11-04.md)
=========================================================================================================================================

`build_system_sky(conn, system_id, limit)` is the whole answer: the stars the system sees at a limiting magnitude (`galaxy.sky.sky_from`, from the
stored stars only), packed (`web/lib/skypack.py`), with the nebulae and supernova remnants it sees beside them (`galaxy.sky_sprites`), and the report of
how much of the sky is built. The packed list's `meta` carries the sprites and the report, so one file (the cache's file) holds everything the page
needs; `docs/design/interface/system-sky-wire-format.md` is the contract.

The sky is a snapshot at the system's own epoch (`star_systems.epoch_unix`, else the orbit simulation's, else 0): the page's orbit clock moves it on
(`apparent = offset + relative velocity * t - velocity * tau`), so a sky stays right across orbital updates and the cache keys it on the stored data,
not on the clock. The observer's velocity is its system's stored one (the stand-in velocity of a star at its place when it has none).

The `version_key` of an answer is `skypack.version_key(stamp)`: the galaxy's tile stamp (`db.query.galaxy_tile_stamp`: sectors, systems, bright stars
and the generator epoch) with `skypack.SKY_VERSION`, which is raised when the answer changes for the same data.
"""

import hashlib
import math
import time

from planetgen.api import ids
from planetgen.db import sky_sprites, store
from planetgen.galaxy import potential as galacticPotential
from planetgen.galaxy import sky, sky_band, sky_sprites as spriteMaker
from planetgen.generation import sky_tiers
from planetgen.web.lib import skypack

MAX_SPRITES = 600
"""int: The most sprites a sky sends, the widest first (stand-in SKV2: the browser draws a quad each)."""


_bands = sky_band.MemoryBandStore()
"""The glow bands this process has worked out, by cache key: a band is shared by every system of a 50 parsec cell."""

BAND_FORMAT = "rgb16f"
"""str: The glow image's wire form: little-endian half floats, red, green and blue a texel, row by row (`sky_band.to_half_bytes`)."""


class SkyError(ValueError):
    """A system with no place in the galaxy, or a galaxy with no shape."""


def _band_galaxy_key(info, galaxy_key):
    """Names the galaxy for the glow band's cache: the seed and the shape's numbers (not its cluster field)."""
    shape = info.shape._replace(components=None, structure=None)
    return f"{galaxy_key}|{hashlib.sha256(repr(shape).encode('utf-8')).hexdigest()[:16]}"


def band_key(info, galaxy_key, observer, limit):
    """The cache key of the glow image of `observer` (its 50 parsec cell) at `limit`: the band's own key and what turns it into the image; also the name of the file the web layer keeps."""
    inner = sky_band.cache_key(sky_band.cell_of(observer.position_pc), _band_galaxy_key(info, galaxy_key), limit)
    image = f"{inner}|{skypack.SKY_VERSION}|{sky_band.SB_REF}|{sky_band.TINT_TEMPERATURE_K}|{sky_band.EQUIRECT_WIDTH}x{sky_band.EQUIRECT_HEIGHT}|{BAND_FORMAT}"
    return hashlib.sha256(image.encode("utf-8")).hexdigest()[:32]


def _setup(conn, system_id):
    info = store.get_galaxy_shape(conn, light=True)
    if info is None:
        raise SkyError("the galaxy has no shape yet")
    pot = galacticPotential.for_galaxy(info.shape, info.galaxy_seed)
    circular_speed = lambda radius: galacticPotential.circular_speed_kms(pot, radius)
    galaxy_key = "galaxy:" + (info.galaxy_seed.hex() if info.galaxy_seed is not None else "none")
    return info, circular_speed, galaxy_key, _observer(conn, system_id, galaxy_key, circular_speed)


def build_system_band(conn, system_id, limit=None):
    """
    The glow image of the sky of `system_id` at `limit` (`skypack.check_limit`): the band of its 50 parsec cell (`galaxy.sky_band.band_for`, kept in this
    process) as the equirectangular image of the page, in half floats.

    Returns:
        dict: `key` (the band's cache key), `limit`, `cell`, `width`, `height`, `format` (`BAND_FORMAT`), `seconds` and `data` (the bytes).

    Raises:
        SkyError: As `build_system_sky`.
        skypack.SkyLimitError: A limit out of range.
    """
    started = time.perf_counter()
    limit = skypack.check_limit(limit)
    info, _speed, galaxy_key, observer = _setup(conn, system_id)
    key = band_key(info, galaxy_key, observer, limit)
    full = store.get_galaxy_shape(conn)
    band = sky_band.band_for(observer.position_pc, full.shape, _band_galaxy_key(info, galaxy_key), _bands, limit)
    image = sky_band.equirect(band)
    return {"key": key, "limit": limit, "cell": list(sky_band.cell_of(observer.position_pc)), "width": image.shape[1], "height": image.shape[0],
            "format": BAND_FORMAT, "seconds": round(time.perf_counter() - started, 3), "data": sky_band.to_half_bytes(image)}


def _observer(conn, system_id, galaxy_key, circular_speed):
    row = conn.execute(
        "SELECT galaxy_x_pc, galaxy_y_pc, galaxy_z_pc, velocity_x_kms, velocity_y_kms, velocity_z_kms, epoch_unix FROM star_systems WHERE id = ?",
        (system_id,)).fetchone()
    if row is None:
        raise SkyError(f"no star_systems row with id {system_id}")
    if row["galaxy_x_pc"] is None:
        raise SkyError("the system has no place in the galaxy")
    position = (float(row["galaxy_x_pc"]), float(row["galaxy_y_pc"]), float(row["galaxy_z_pc"]))
    epoch = row["epoch_unix"]
    if epoch is None:
        state = conn.execute("SELECT UNIX_TIMESTAMP(last_updated_at) AS unix FROM orbit_simulation_state WHERE id = 1").fetchone()
        epoch = state["unix"] if state is not None and state["unix"] is not None else 0
    velocity = (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])
    if any(v is None for v in velocity):
        velocity = tuple(float(v) for v in sky.unbuilt_velocity_kms([position], [f"system:{system_id}"], galaxy_key, circular_speed)[0])
    return sky.Observer(position, float(epoch), tuple(float(v) for v in velocity), int(system_id))


def _wanted_pc(conn, observer, limit):
    """How far the stored stars must be complete for this observer and limit: the reach of a star of the galaxy scatter's own floor (`sky_tiers.reach_pc`,
    through the dust of the clearest line of sight), past which the pre-placed scatter alone supplies every star that shows. `None` for a galaxy with no
    scatter level. It ignores the budget that drops a band (`sky_tiers.build`), which costs a minute to work out; VIEW.25 sharpens it."""
    scatter = store.bright_star_scatter_settings(conn)
    if scatter is None or not float(scatter[0]) > 0.0:
        return None
    curve = sky_tiers.extinction_curve(observer.position_pc)
    return float(sky_tiers.reach_pc(float(scatter[0]), limit, curve))


def build_system_sky(conn, system_id, limit=None, tile_stamp=""):
    """
    The sky of the system `system_id` at the limiting magnitude `limit`.

    Args:
        conn (store.Connection): An open, read-only connection.
        system_id (int): The `star_systems.id`.
        limit (float, optional): The limiting magnitude (`skypack.check_limit`); the default when `None`.
        tile_stamp (str): The galaxy's tile stamp, for the version key.

    Returns:
        dict: `system` (the printed ID), `limit`, `epoch_unix`, `version_key`, `count` (stars), `sprites` (the number of sprite records), `report`
            (see `docs/design/interface/system-sky-wire-format.md` section 4), `seconds` (what the build took) and `packed` (the bytes of the packed
            list, whose `meta` holds `nebulae` and the report).

    Raises:
        SkyError: As above.
        skypack.SkyLimitError: A limit out of range.
    """
    started = time.perf_counter()
    limit = skypack.check_limit(limit)
    info, circular_speed, galaxy_key, observer = _setup(conn, system_id)
    at_unix = observer.epoch_unix

    stars = sky.sky_from_database(conn, observer, at_unix, limit, galaxy_key, circular_speed)
    built = stars.report["built_to_pc"]
    objects = sky_sprites.objects_in_sphere(conn, observer.position_pc, spriteMaker.SEARCH_RADIUS_PC)
    printed = {}
    for kind in (sky_sprites.NEBULA, sky_sprites.REMNANT):
        rows = objects.ident[objects.kind == kind]
        for row, text in ids.printed_many(conn, kind, rows).items():
            printed[(kind, int(row))] = text
    records, sprite_report = spriteMaker.sprites_from(observer, objects, at_unix, galaxy_key, circular_speed, ids=printed)
    records.sort(key=lambda r: (-r["angularRadiusRad"], r["id"]))
    sprite_report["dropped_for_count"] = max(0, len(records) - MAX_SPRITES)
    records = records[:MAX_SPRITES]
    wanted = _wanted_pc(conn, observer, limit)

    report = {key: stars.report[key] for key in ("limit", "stars_in", "tiers", "excluded_own_system", "no_photometry", "unknown_epoch", "not_formed_yet",
                                                 "other_phase", "unbuilt_velocity", "visible")}
    report["built_to_pc"] = built
    complete = stars.report["complete_above_sol"]
    report["complete_above_sol"] = complete if math.isfinite(complete) else None
    report["wanted_pc"] = wanted
    report["filling"] = bool(wanted is not None and built < wanted)
    report["sprites"] = sprite_report
    key = skypack.version_key(tile_stamp)
    printed_system = ids.printed(conn, "system", system_id) or f"system:{system_id}"
    band = {"key": band_key(info, galaxy_key, observer, limit), "width": sky_band.EQUIRECT_WIDTH, "height": sky_band.EQUIRECT_HEIGHT, "format": BAND_FORMAT,
            "href": f"/system/{printed_system}/sky/band?limit={limit:g}"}
    meta = {"system": printed_system, "version_key": key, "limit": limit, "galaxy": galaxy_key, "nebulae": records, "band": band, "report": report}
    packed = skypack.pack_sky(stars, at_unix, limit, observer.velocity_kms, meta=meta)
    return {"system": printed_system, "limit": limit, "epoch_unix": at_unix, "version_key": key, "count": int(len(stars.uid)), "sprites": len(records),
            "report": report, "seconds": round(time.perf_counter() - started, 3), "packed": packed}
