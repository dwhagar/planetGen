# planetgen/generation/abundance_tier.py

"""
The Class M abundance tier of a galaxy (Boss, 11:19 PM EST 2026-11-13; docs/design/habitability/class-m-abundance-tiers.md)
==============================================================================================================

Known science cannot say how common an Earth is: the generator's own laws, with every bug fixed, give under one Class M world in ten
thousand Sun-like systems (reports/class-m-blue-levers-2026-11-14.md). A galaxy therefore carries a tier, one word, that says how far past
the science its makers want its rocky planets in the ecosphere pushed toward Earth's conditions:

| Tier | Target, Class M worlds per Sun-like system |
|---|---|
| `science` (the default) | whatever the laws give: nothing is moved |
| `sparse` | 1 in 1,000 |
| `scattered` | 1 in 100 |
| `plentiful` | 1 in 30 |
| `abundant` | 1 in 10 |

A tier turns three dials and nothing else (`tuning.ABUNDANCE_TIER_DIALS`, `Dials`):

1. `water_share`: the share of candidates (`candidate`) whose water class is drawn as Earth's (`low` of `physics.delivered_water`) and not from
   the Bern law at the distance the planet formed at.
2. `bare_share`: the share of candidates whose starting hydrogen envelope is taken as already lost (a bare core).
3. `request_share`: the share of systems that are given the request for a habitable world (`SystemConfig.HABITABLE_WORLD`), which raises the
   odds and never forces (`util.request_odds`).

No dial touches a class limit, names a class letter or sets a count: each one moves a draw, and the class of every body is still derived
from its finished conditions. The targets are what the dial values are calibrated to by census; the tier never checks or enforces them.

**Science changes nothing.** On `science` every dial is 0, every function here says so without a draw, and a caller must make no draw of its
own: a seeded system comes out byte for byte as it did before the tier existed.

The tier reaches a system through its config (`SystemConfig.ABUNDANCE_TIER`), which a galaxy sector sets from the galaxy's stored choice
(`brightStars.FillContext.apply_population`), as it sets the population and the galaxy's shape.
"""

from collections import namedtuple

from planetgen import tuning

TIERS = ("science", "sparse", "scattered", "plentiful", "abundant")
"""tuple[str]: Every tier, from no adjustment up: the choices of the `class_m` lever (`galaxy.dials`), in its order."""

DEFAULT = TIERS[0]
"""str: The tier of a galaxy that was told none, and of work outside any galaxy."""

TARGET_SYSTEMS_PER_CLASS_M = {"science": None, "sparse": 1000, "scattered": 100, "plentiful": 30, "abundant": 10}
"""dict: Tier -> the Sun-like systems per Class M world its dial values are calibrated to (Boss, 11:19 PM EST 2026-11-13); `None` for
`science`, which has no target."""

Dials = namedtuple("Dials", ["water_share", "bare_share", "request_share"])
"""The three dials of a tier, each a share from 0 (off) to 1."""

OFF = Dials(0.0, 0.0, 0.0)
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
    """The share of candidates given Earth's water class under the tier of `body`: 0 for a body that is not a `candidate`, and on `science`."""
    return dials_of(body).water_share if candidate(body) else 0.0


def bare_share(body):
    """The share of candidates whose starting envelope is taken as lost under the tier of `body`: 0 for a body that is not a `candidate`, and on `science`."""
    return dials_of(body).bare_share if candidate(body) else 0.0


def request_share(config):
    """The share of systems given the request for a habitable world under the tier of `config`; 0 on `science`."""
    return dials_of(config).request_share
