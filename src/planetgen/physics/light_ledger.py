# planetgen/physics/light_ledger.py

"""
The light and ultraviolet rows of the effect ledger (GEN.1465, GEN.1498)
========================================================================

What each star puts on each planet and moon as light and as far ultraviolet, by source: the parts of the field at the distance the body's light was made on
(`radiation_field.field_at`, GEN.1441). `fill` notes them in an `EffectLedger`; it draws nothing and changes no value.

Two callers, one rule. Generation (`StarSystem._fill_light_ledger`) has the distance the body's air was solved on (`_air_light_au`, kept on the body and not stored). The refill job
(`store.refresh_effect_ledger`, GEN.1498) works from a stored system, which keeps only where the orbit settled: `stored_distance_au` is the same circle of mean light worked out from
the stored distance and eccentricity (a moon's from its planet's, GEN.727). The two differ only by the little a later settling moved the orbit, so a refilled row may differ from the
one generation wrote in its last digits, and the sources and their order are the same.
"""

from planetgen.physics import orbit_elements, radiation_field

EFFECTS = (("light", lambda field: field.bolometric_w_m2), ("uv", lambda field: field.far_uv_ratio))
"""tuple: `(effect name, the band of a `radiation_field.Field` it reads)` for each effect this module fills."""

EFFECT_NAMES = tuple(name for name, _band in EFFECTS)
"""tuple: The names of the effects this module fills, which the refill replaces together."""


def stored_distance_au(body, planet=None):
    """
    The distance, AU, of the circle that gives the body the mean light of its stored orbit: the planet's distance and eccentricity, and for a moon its planet's (`planet`).
    `None` for a body with no distance (a belt, which has no air-solved light).
    """
    if getattr(body, "body_type", None) == "a":
        return None
    source = planet if planet is not None else body
    axis = getattr(source, "distance", None)
    if not axis:
        return None
    return orbit_elements.flux_distance_au(float(axis), getattr(source, "eccentricity", None))


def fill(ledger, planets, distance_of):
    """
    Notes in `ledger` the light and ultraviolet each star gives each body of `planets` and their moons.

    Args:
        ledger (EffectLedger): Where the contributions go.
        planets (iterable): The planets and belts of the system (`StarSystem.all_planets`).
        distance_of (callable): `distance_of(body, planet)` gives the distance, AU, the body's light is made on, or a falsy value for a body with none (it is skipped); `planet` is the
            body itself for a planet and its planet for a moon.
    """
    for obj in planets:
        for body in (obj, *getattr(obj, "moons", ())):
            distance = distance_of(body, obj)
            host = getattr(body, "star", None)
            if not distance or host is None:
                continue
            field = radiation_field.field_at(host, distance)
            stars = {companion.member_no: companion.star for companion in radiation_field.companions_of(host)}
            stars[None] = host
            for effect, band_of in EFFECTS:
                for member_no, value in band_of(field).by_source.items():
                    if value > 0.0 and member_no in stars:
                        ledger.add(effect, body, stars[member_no], value)
