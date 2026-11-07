# planetgen/physics/envelope_share.py

"""
The envelope share of every body: the mass of a hydrogen-rich air over the mass of the body (GEN.1226)
=====================================================================================================

`envelope_mass_fraction` is the share of a body's mass that is a hydrogen-helium envelope. A giant stores the envelope its interior model
holds (`giant_draw`). A rocky planet or moon was drawn with none (`bulk_draw`), but a small world can keep a deep hydrogen air (the
class T gas dwarfs), and the body types (`planet_type`, GEN.1219) cut on this share, so a hydrogen-air dwarf must store the share its air
holds, not 0. Boss, 1:56 PM EDT 2026-10-31: store the true air share, the hydrogen-rich air's mass over the body's mass, not a class-based stand-in.

For a rocky body the air column is its surface pressure times the area over the gravity, `4 pi R^4 P / (G M)`, over the mass `M`:

    share = 4 pi R^4 P / (G M^2)      (the air is hydrogen-rich when `HYDROGEN_RICH_SHARE` of its pressure is hydrogen)

A body with an air that is not hydrogen-rich (nitrogen, carbon dioxide, a thin air of any kind) keeps 0, which is honest: its air is secondary
and not an envelope. The four solid fractions take the air's share out of themselves, so the five fractions still add up to 1
(`bulk_composition`); the radius, density and gravity were worked out from the solids and are not changed. `restore` puts the drawn
fractions back before the air is made again (`planets.refinish_rocky` solves the air once more at the star's final age), so the share never
compounds.

Use:

    envelope_share.restore(planet)     # before the air is solved: the drawn fractions
    envelope_share.fill(planet)        # after: the true air share, the solids scaled to match
"""

import math

from planetgen.physics import atmosphere, constants, gases
from planetgen.physics import unit_registry as units

HYDROGEN_RICH_SHARE = 0.3
"""float: STAND-IN ES1 (L0, the class T line of GEN.1142, registry CL1): an air with at least this share of its pressure as hydrogen is a hydrogen-rich air, whose
whole mass is an envelope; a smaller share is a secondary air with no envelope."""

SOLIDS = ("metal_mass_fraction", "rock_mass_fraction", "water_mass_fraction", "ice_mass_fraction")
"""tuple: The body attributes of the solid fractions that give way to the envelope (the ice species beside `ice_mass_fraction` are shares of the ice and do not change)."""

_km_to_m = units.converter("km", "m")

_FORMED = "_formation_fractions"


def air_share(body):
    """The hydrogen-rich air's mass over the body's mass; 0 for a body without a hydrogen-rich air (a giant, a body with no air, no mass or no radius)."""
    if getattr(body, "_captured_envelope", 0.0):
        return 0.0      # the envelope the planet captured (GEN.1323, `envelope_capture`) is its hydrogen air already; the share is not counted twice
    mass, radius_km, pressure = getattr(body, "mass", None), getattr(body, "radius", None), getattr(body, "atmospheric_pressure", None)
    if not mass or not radius_km or not pressure or pressure <= 0.0:
        return 0.0
    partials = atmosphere.gas_partials_pa(body)
    total = sum(partials.values())
    if total <= 0.0 or partials.get(gases.stored("h2").formula, 0.0) < HYDROGEN_RICH_SHARE * total:
        return 0.0
    radius = _km_to_m(radius_km)
    return min(4.0 * math.pi * radius ** 4 * pressure / (constants.G * mass * mass), 0.5)


def set_gas_detail(body):
    """Stores the helium and heavy-element shares of the envelope of a rocky body (`envelope_helium_fraction`, `envelope_z_fraction`): the disk's gas, so the star's
    (`bulk_draw.host_ratios`); `None` for both when the body has no envelope. A body with an envelope share stores both (the bulk rule of `tests/audit/sweep_bulk.py`)."""
    from planetgen.physics import bulk_draw      # imported here: bulk_draw is loaded by the planet before this module is used
    if not (getattr(body, "envelope_mass_fraction", None) or 0.0) > 0.0:
        body.envelope_helium_fraction = body.envelope_z_fraction = None
        return
    star = getattr(body, "star", None)
    if star is None:
        return              # a bare stand-in body of a test: no star, no ratios
    ratios = bulk_draw.host_ratios(star)
    body.envelope_helium_fraction = ratios.y_envelope * (1.0 - ratios.z)
    body.envelope_z_fraction = ratios.z


def restore(body):
    """Puts a rocky body's drawn fractions back (what `fill` scaled), before its air is solved again. Nothing to put back on a body not filled yet."""
    formed = getattr(body, _FORMED, None)
    if formed is None:
        return
    for name, value in formed.items():
        setattr(body, name, value)
    delattr(body, _FORMED)      # the next fill takes the fractions of the air solved again as the drawn ones


def fill(body):
    """Sets `envelope_mass_fraction` of a rocky planet or moon to its air's share and takes it out of the solid fractions. A giant keeps the envelope its interior holds."""
    if getattr(body, "body_type", None) != "t" or getattr(body, "envelope_mass_fraction", None) is None:
        return
    formed = getattr(body, _FORMED, None)
    if formed is None:
        formed = {name: getattr(body, name) for name in (*SOLIDS, "envelope_mass_fraction")}
        setattr(body, _FORMED, formed)
    share = air_share(body)
    body.envelope_mass_fraction = formed["envelope_mass_fraction"] + share
    set_gas_detail(body)
    keep = 1.0 - share
    for name in SOLIDS:
        if formed[name] is not None:
            setattr(body, name, formed[name] * keep)
