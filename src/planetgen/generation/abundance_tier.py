# planetgen/generation/abundance_tier.py

"""
The Class M abundance tier of a galaxy (GEN.1704; Boss, 11:19 PM EST 2026-11-13; docs/design/habitability/class-m-abundance-tiers.md)
========================================================================================================================================

Known science cannot say how common an Earth is: the generator's own laws, with every bug fixed, give under one Class M world in ten
thousand Sun-like systems (reports/class-m-blue-levers-2026-11-14.md). A galaxy therefore carries a tier, one word, that says how far past
the science its makers want its rocky planets in the ecosphere pushed toward Earth's conditions:

| Tier | Target, Class M worlds per Sun-like system |
|---|---|
| `science` (the default) | whatever the laws give: nothing is moved |
| `sparse` | 1 in 500 |
| `scattered` | 1 in 135 |
| `plentiful` | 1 in 37 |
| `abundant` | 1 in 10 |

A tier turns six dials and nothing else (`tuning.ABUNDANCE_TIER_DIALS`, `Dials`):

1. `water_share`: the share of candidates (`candidate`) whose water class is given by the tier and not drawn from the Bern law at the
   distance the planet formed at.
2. `bare_share`: the share of candidates whose starting hydrogen envelope is taken as already lost (a bare core).
3. `request_share`: the share of systems that are given the request for a habitable world (`SystemConfig.HABITABLE_WORLD`), which raises the
   odds and never forces (`util.request_odds`).
4. `dry_part` (GEN.1721; Boss, 8:24 PM EST 2026-11-14: desert, ocean, glaciated and hothouse worlds beside the temperate one, the temperate
   one in front): the part of the planets the water dial reaches that take the dry class (`dry` of `physics.delivered_water`, the water of
   a desert world) and not Earth's (`low`). It splits the water dial's one draw and makes none of its own; the rest, the larger part, take
   Earth's class. A desert world's water is drawn inside `tuning.ABUNDANCE_DESERT_WATER_RANGE`, the part of the dry class that leaves a few small
   seas, and its small sea is an ordinary salt sea, not the concentrated brine of a seep (`desert_sea`; `hydrosphere.salts_concentrated`).
5. `quiet_share` (GEN.1721; Boss, 9:43 PM EST 2026-11-14: "tweak the volcanic activity a bit"): the share of candidates whose interior is quiet, so
   that its volcanic gases take `tuning.ABUNDANCE_QUIET_SINK_FACTOR` of the oxygen a body of the same heat loses (`quiet_interior`;
   `biogenic_oxygen.reductant_flux`). One draw per candidate, made once when the planet is first drawn.
6. `cycle_share` (GEN.1720; Boss, 11:58 PM EST 2026-11-14, "recommendation approved"): the share of candidates whose two carbon-cycle constants, the exponent of the heat flow in the outgassing and the
   factor of the weathering law's thermodynamic limit, are drawn from the part of their published ranges that holds the balance nearest the Earth's (`tuning.ABUNDANCE_CYCLE_OUTGASSING_RANGE`,
   `ABUNDANCE_CYCLE_LAMBDA_RANGE`; `continental_weathering.draw_cycle`). The balance goes as the outgassing to about the fifth power, so over the whole ranges most Earth-shaped planets hold under 2 Pa
   of CO2 or over 70 kPa (reports/gen1720-co2-attribution-2026-11-14.md). One draw per candidate under a tier, none on `science`.

No dial touches a class limit, names a class letter or sets a count: each one moves a draw, and the class of every body is still derived
from its finished conditions. The targets are what the dial values are calibrated to by census; the tier never checks or enforces them.

**Science changes nothing.** On `science` every dial is 0, every function here says so without a draw, and a caller must make no draw of its
own: a seeded system comes out byte for byte as it did before the tier existed.

The tier reaches a system through its config (`SystemConfig.ABUNDANCE_TIER`, a class attribute that is set on a config only for a tier other than science), which a galaxy sector sets from the galaxy's stored choice
(`brightStars.FillContext.apply_population`), as it sets the population and the galaxy's shape.
"""

from collections import namedtuple

from planetgen import tuning

TIERS = ("science", "sparse", "scattered", "plentiful", "abundant")
"""tuple[str]: Every tier, from no adjustment up: the choices of the `class_m` lever (`galaxy.dials`), in its order."""

DEFAULT = TIERS[0]
"""str: The tier of a galaxy that was told none, and of work outside any galaxy."""

TARGET_SYSTEMS_PER_CLASS_M = {"science": None, "sparse": 500, "scattered": 135, "plentiful": 37, "abundant": 10}
"""dict: Tier -> the Sun-like systems per Class M world its dial values are calibrated to (Boss, 9:28 and 9:30 PM EST 2026-11-14, after the laws alone
measured 1 in 2,800: `sparse` 1 in 500, `abundant` 1 in 10 and the two between split evenly on a log scale, each step about 3.7 times; they replace his
first set of 1,000, 100, 30 and 10 of 11:19 PM EST 2026-11-13); `None` for `science`, which has no target."""

Dials = namedtuple("Dials", ["water_share", "bare_share", "request_share", "dry_part", "quiet_share", "cycle_share"])
"""The six dials of a tier, each a share from 0 (off) to 1."""

OFF = Dials(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
"""Dials: Every dial off: `science`."""


def tier_of(subject):
    """
    The tier `subject` is made under: a `SystemConfig`, or anything that carries one as `system_config` (a planet, a moon, a star system).
    `DEFAULT` when it carries none.

    Raises:
        ValueError: The config names a tier that does not exist.
    """
    config = getattr(subject, "system_config", subject)
    tier = getattr(config, "ABUNDANCE_TIER", None) or DEFAULT
    if tier not in TIERS:
        raise ValueError(f"{tier!r} is not one of {', '.join(TIERS)}")
    return tier


def dials_of(subject):
    """The `Dials` of the tier `subject` is made under (`tier_of`); `OFF` on `science`."""
    tier = tier_of(subject)
    return OFF if tier == DEFAULT else Dials(*tuning.ABUNDANCE_TIER_DIALS[tier])


def active(subject):
    """Whether `subject` is made under a tier other than `science`. A caller that draws for a dial checks this (or the dial's own share) first, so `science` draws nothing."""
    return tier_of(subject) != DEFAULT


def candidate(body):
    """
    Whether the dials act on `body`: a rocky planet (`body_type` "t", never a moon) whose orbit lies in its star's ecosphere (`zone` "e").
    The one definition the water and envelope dials share, and the first row of the meter's funnel.
    """
    return getattr(body, "body_type", None) == "t" and not getattr(body, "is_moon", False) and getattr(body, "zone", None) == "e"


def water_share(body):
    """The share of candidates whose water class the tier of `body` gives (Earth's or the dry one, `dry_part`): 0 for a body that is not a `candidate`, and on `science`."""
    return dials_of(body).water_share if candidate(body) else 0.0


def dry_part(body):
    """The part of the planets the water dial reaches that take the dry class under the tier of `body` (GEN.1721): 0 for a body that is not a `candidate`, and on `science`."""
    return dials_of(body).dry_part if candidate(body) else 0.0


def desert_sea(body):
    """Whether `body` is a desert world a tier made (its water class was dialed to the dry one: `body._water_dialed`, set by `bulk_draw` under a tier only), whose small sea is an ordinary salt sea, not a brine seep (GEN.1721; `hydrosphere.salts_concentrated`)."""
    return getattr(body, "_water_dialed", None) == "dry"


def quiet_share(body):
    """The share of candidates whose interior is quiet under the tier of `body` (GEN.1721): 0 for a body that is not a `candidate`, and on `science`."""
    return dials_of(body).quiet_share if candidate(body) else 0.0


def quiet_interior(body):
    """Whether a tier drew `body` with a quiet interior (`body._quiet_interior`, set by `Planet._draw_body` under a tier only; never on `science`)."""
    return bool(getattr(body, "_quiet_interior", False))


def cycle_share(body):
    """The share of candidates whose carbon-cycle constants are drawn from the part of their ranges nearest the Earth's balance under the tier of `body` (GEN.1720): 0 for a body that is not a `candidate`, and on `science`."""
    return dials_of(body).cycle_share if candidate(body) else 0.0


def bare_share(body):
    """The share of candidates whose starting envelope is taken as lost under the tier of `body`: 0 for a body that is not a `candidate`, and on `science`."""
    return dials_of(body).bare_share if candidate(body) else 0.0


def request_share(config):
    """The share of systems given the request for a habitable world under the tier of `config`; 0 on `science`."""
    return dials_of(config).request_share
