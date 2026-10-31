# planetgen/physics/body_colour.py

"""
The atmosphere colour a planet or moon stores (DB.96; docs/design/composition/atmosphere-colour.md section 14.1).

`atmosphere_colour` is the sRGB colour of the body's air as one integer, `red * 65536 + green * 256 + blue`, the hue the picture of the body
is painted with (MAP.244 reads it). `generate` works it out once, from what the body is made of by then: the air (`air_column`), its
star's temperature, the ozone column (`radiation.generate_dose`) and the haze, dust and cloud (`aerosols`, DB.68), over the ground
its stored surface mix makes (GEN.462, `ground_albedo`; the stand-in `atmosphere_colour.STAND_IN_GROUND_ALBEDO` where no mix is stored). A solid world gets the colour of the light its air
scatters down to the ground (GEN.418); a gas or ice giant the colour of its air at the reference level (GEN.419). A body with no air, no star
temperature or no mixture stores `None`.

The cloud cover is not drawn here: `aerosols.BODY_FIELDS` already stores it (`cloud_cover`, DB.68), and the picture reads that column.
"""

import numpy as np

from planetgen.physics import aerosols, air_column, atmosphere_colour, giant_colour, photometry, surface_composition

BODY_FIELDS = ("atmosphere_colour",)
"""tuple: What `generate` sets on a planet or moon, as stored in `planets` and `moons` (DB.96)."""

_CENTRES = atmosphere_colour.WAVELENGTHS_NM


def rgb_integer(rgb8):
    """The 24-bit integer of an 8-bit sRGB `(r, g, b)`."""
    red, green, blue = rgb8
    return (red << 16) | (green << 8) | blue


def ground_albedo(body):
    """
    The reflectance of a solid body's ground, one value per wavelength of the colour grid (GEN.462): the lab reflectance in the three visible
    bands (450, 550 and 650 nm centres; held flat beyond them) of the stored areal mix of ground types and ice, area-weighted. The part of the
    surface that has no measured visible spectrum (the open ocean and the metal, `surface_composition.band_albedo`) takes
    `atmosphere_colour.STAND_IN_GROUND_ALBEDO`, and so does the whole ground of a body whose mix is not stored. A surface under air is not
    space-weathered, so no darkening enters here.
    """
    stand_in = atmosphere_colour.STAND_IN_GROUND_ALBEDO
    mix = surface_composition.mix_of(body)
    if not mix:
        return stand_in
    composition = surface_composition.SurfaceComposition(
        ground_type=getattr(body, "ground_type", None) or max(mix, key=mix.get), mix=mix, components={},
        grain=getattr(body, "ground_grain", None) or surface_composition.POWDER)
    bands, covered = surface_composition.band_albedo(composition)
    values = [bands[band] for band in surface_composition.BANDS]
    if not covered > 0.0 or any(value is None for value in values):
        return stand_in
    blended = [covered * value + (1.0 - covered) * stand_in for value in values]
    return np.interp(_CENTRES, surface_composition.BAND_CENTRES_NM, blended)


def _solid_colour(body):
    temperature = getattr(getattr(body, "star", None), "temperature", None)
    if not temperature or not photometry.gives_visible_light(temperature):
        return None
    column = air_column.column_of(body)
    if column is None:
        return None
    stored = aerosols.of_body(body)
    haze, dust, cloud = aerosols.colour_layers(stored, _CENTRES)
    result = atmosphere_colour.atmosphere_colour(
        column, temperature_k=temperature, ozone_du=getattr(body, "ozone_column_du", None) or 0.0,
        ground_albedo=ground_albedo(body), haze=haze, dust=dust, cloud=cloud, cloud_cover=stored.cloud_cover,
        light=atmosphere_colour.body_light(body.star))
    return result.colour.rgb8


def colour_of(body):
    """The `(r, g, b)` of a planet or moon's atmosphere colour, or `None` for a body that has none (no air, no star temperature, no mixture)."""
    try:
        if getattr(body, "body_type", "t") == "g":
            found = giant_colour.body_giant_colour(body)
            return None if found is None else found.atmosphere.colour.rgb8
        return _solid_colour(body)
    except photometry.NoVisibleLight:  # a black hole's or a dead remnant's light holds none between 380 and 780 nm: nothing to colour
        return None


def generate(body):
    """Sets `BODY_FIELDS` on a planet or moon (call after its ozone column and aerosols). A body with no colour stores `None`."""
    rgb8 = colour_of(body)
    body.atmosphere_colour = None if rgb8 is None else rgb_integer(rgb8)
