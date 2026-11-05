# planetgen/galaxy/sky_sprites.py

"""
The nebula and remnant sprites of a system's sky (VIEW.31, the sprite list of the sky route; the physics is VIEW.28's `physics/nebula_sprites.py`)
================================================================================================================================================

`sprites_from(observer, objects, at_unix)` takes the stored nebulae and supernova remnants of a sphere (`db.sky_sprites.SpriteObjects`) and returns the
records the browser draws (`html/static/skysprites.js`, `docs/design/interface/system-sky-wire-format.md` section 3), each:

    {"id", "kind", "family", "direction": [x, y, z], "distanceLy", "angularRadiusRad", "inside", "extinctionMag", "surfaceBrightness", "visible", "seen"}

Per object, as the stars of `galaxy/sky.py` are:

1. **Its velocity.** One with none stored gets the stand-in velocity of a star at its place (`sky.unbuilt_velocity_kms`, the same stand-in SKC1; a
   cloud rides the galaxy's rotation). One with no epoch holds at the sky's time.
2. **Where it appears.** `light_travel.retarded_offsets`: the centre where it was when its light left it. `direction` is the unit vector of that offset
   in the galaxy frame, `distanceLy` its length.
3. **What it looks like.** `nebula_sprites.sprite_data`: the angular radius, whether the observer is inside, the dimming through its middle, the mean surface
   brightness (V mag per square arcsecond, `null` for no light of its own), `visible` (under the naked-eye limit) and `seen` (a remnant's explosion has been
   seen). The surface brightness is dimmed by the dust of the disc between the observer and the object (`sky_photometry.galactic_extinction_av`).
4. **Whether it is worth sending.** A record is kept when it is wider than `MIN_ANGULAR_RADIUS_RAD` and either shines under the limit (`visible`) or dims
   what lies behind it by at least `MIN_EXTINCTION_MAG`, and its explosion has been seen. The counts of what was left out are in the report.
"""

import math

import numpy as np

from planetgen.galaxy import sky
from planetgen.physics import light_travel, nebula_sprites, sky_photometry
from planetgen.physics import unit_registry as units

MIN_ANGULAR_RADIUS_RAD = 0.002
"""float: A sprite narrower than this (radius, radians; 0.115 degrees, a quarter of the Moon's width) is not sent (stand-in SPR1)."""

MIN_EXTINCTION_MAG = 0.02
"""float: A cloud with no light of its own is sent only if it dims what lies behind it by this much (magnitudes; the browser draws none under it, stand-in SPR2)."""

SEARCH_RADIUS_PC = 10000.0
"""float: How far from the observer the sprites are looked for, parsecs (stand-in SPR3: the nebulae and remnants of the stored galaxy are few enough that the
search is bounded by the sphere query, and a sprite narrower than `MIN_ANGULAR_RADIUS_RAD` is left out anyway)."""

_ly_to_pc = units.converter("ly", "pc")
_pc_to_ly = units.converter("pc", "ly")
_pc_to_m = units.converter("pc", "m")
_m_to_pc = units.converter("m", "pc")
_km_to_m = units.converter("km", "m")
_yr_to_s = units.converter("yr", "s")


def _number(value):
    """A float for JSON: a finite number, else `None`."""
    value = float(value)
    return value if math.isfinite(value) else None


def sprites_from(observer, objects, at_unix, galaxy_key="", circular_speed_kms=None, dust=None, ids=None):
    """
    The sprites of the nebulae and remnants `objects` (a `db.sky_sprites.SpriteObjects`) seen from `observer` (a `sky.Observer`) at the time `at_unix`.

    Args:
        observer (sky.Observer): Where the sky is seen from.
        objects (SpriteObjects): The stored objects in the search sphere.
        at_unix (float): The sky's time, seconds.
        galaxy_key (str): Names the galaxy for the stand-in velocities.
        circular_speed_kms: A function of radius in parsecs for them (the running galaxy's by default, as `sky.sky_from_stars`).
        dust (dict): Overrides of the dust law.
        ids (dict, optional): `{(kind, row id): printed ID}`; a record's `id` is the printed ID when there is one, else `"<kind>:<row id>"`.

    Returns:
        tuple: `(records, report)`: the records (a list of dicts, in the order of `objects`) and `{"in", "sent", "too_small", "too_faint", "not_seen", "inside"}`.
    """
    n = len(objects.ident)
    report = {"in": int(n), "sent": 0, "too_small": 0, "too_faint": 0, "not_seen": 0, "inside": 0}
    if n == 0:
        return [], report
    dust = dict(dust or {})
    epoch = np.asarray(objects.epoch_unix, dtype=float).copy()
    epoch[~np.isfinite(epoch)] = at_unix
    position = np.asarray(objects.position_pc, dtype=float)
    velocity = np.asarray(objects.velocity_kms, dtype=float).copy()
    unknown = ~np.all(np.isfinite(velocity), axis=1)
    if unknown.any():
        if circular_speed_kms is None:
            from planetgen.galaxy import potential as galacticPotential
            pot = galacticPotential.current()
            circular_speed_kms = lambda radius: galacticPotential.circular_speed_kms(pot, radius)
        keys = np.array([f"{k}:{i}" for k, i in zip(objects.kind, objects.ident)], dtype=object)
        velocity[unknown] = sky.unbuilt_velocity_kms(position[unknown], keys[unknown], galaxy_key, circular_speed_kms)

    obs_position = np.asarray(observer.position_pc, dtype=float).reshape(3)
    obs_velocity = np.asarray(observer.velocity_kms, dtype=float).reshape(3)
    obs_at_object = obs_position[None, :] + sky._km_to_pc_per_s(obs_velocity)[None, :] * (epoch - float(observer.epoch_unix))[:, None]
    offsets, light_time_s = light_travel.retarded_offsets(
        _pc_to_m(position - obs_at_object), _km_to_m(velocity), _km_to_m(obs_velocity), at_unix - epoch)
    distance_pc = _m_to_pc(np.sqrt(np.sum(offsets * offsets, axis=1)))
    distance_ly = _pc_to_ly(distance_pc)
    direction = offsets / np.where(distance_pc > 0.0, _pc_to_m(distance_pc), 1.0)[:, None]
    foreground = np.zeros(n)
    far = distance_pc > 0.0
    if far.any():
        foreground[far] = np.concatenate([
            sky_photometry.galactic_extinction_av(obs_position, direction[far][start:start + sky.DUST_CHUNK], distance_pc[far][start:start + sky.DUST_CHUNK], **dust)
            for start in range(0, int(far.sum()), sky.DUST_CHUNK)])
    # A remnant's age is at its epoch; the sky's time is later (or earlier), and the explosion is seen when its age at the retarded time is positive.
    age = np.asarray(objects.age_years, dtype=float) + (at_unix - epoch) / _yr_to_s(1.0)
    data = nebula_sprites.sprite_data(objects.family, np.asarray(objects.radius_ly, dtype=float), distance_ly, np.asarray(objects.density_cm3, dtype=float),
                                      np.asarray(objects.extinction_av, dtype=float), foreground, age, light_time_s)
    theta = data["angular_radius_rad"]
    small = theta < MIN_ANGULAR_RADIUS_RAD
    unseen = ~data["seen"]
    shines = data["visible"]
    dims = np.asarray(data["extinction_mag"]) >= MIN_EXTINCTION_MAG
    keep = ~small & ~unseen & (shines | dims)
    report.update(too_small=int((small & ~unseen).sum()), not_seen=int(unseen.sum()), too_faint=int((~small & ~unseen & ~shines & ~dims).sum()),
                  inside=int((keep & data["inside"]).sum()), sent=int(keep.sum()))
    records = []
    for i in np.nonzero(keep)[0]:
        kind, row = str(objects.kind[i]), int(objects.ident[i])
        records.append({
            "id": (ids or {}).get((kind, row)) or f"{kind}:{row}", "kind": kind, "family": str(objects.family[i]),
            "direction": [round(float(v), 7) for v in direction[i]], "distanceLy": round(float(distance_ly[i]), 3),
            "angularRadiusRad": float(theta[i]), "inside": bool(data["inside"][i]), "extinctionMag": round(float(data["extinction_mag"][i]), 4),
            "surfaceBrightness": _number(data["surface_brightness"][i]), "visible": bool(shines[i]), "seen": bool(data["seen"][i])})
    return records, report
