# planetgen/physics/body_colour.py

"""
The atmosphere colour a planet or moon stores (DB.96; docs/design/composition/atmosphere-colour.md section 14.1).

`atmosphere_colour` is the sRGB colour of the body's air as one integer, `red * 65536 + green * 256 + blue`, the hue the picture of the body
is painted with (MAP.244 reads it). `generate` works it out once, from what the body is made of by then: the air (`air_column`), its
star's temperature, the ozone column (`radiation.generate_dose`) and the haze, dust and cloud (`aerosols`, DB.68), over the stand-in ground
(`atmosphere_colour.STAND_IN_GROUND_ALBEDO`) until the ground's own reflectance is stored. A solid world gets the colour of the light its air
scatters down to the ground (GEN.418); a gas or ice giant the colour of its air at the reference level (GEN.419). A body with no air, no star
temperature or no mixture stores `None`.

The cloud cover is not drawn here: `aerosols.BODY_FIELDS` already stores it (`cloud_cover`, DB.68), and the picture reads that column.
"""

from planetgen.physics import aerosols, air_column, atmosphere_colour, giant_colour

BODY_FIELDS = ("atmosphere_colour",)
"""tuple: What `generate` sets on a planet or moon, as stored in `planets` and `moons` (DB.96)."""

_CENTRES = atmosphere_colour.WAVELENGTHS_NM


def rgb_integer(rgb8):
    """The 24-bit integer of an 8-bit sRGB `(r, g, b)`."""
    red, green, blue = rgb8
    return (red << 16) | (green << 8) | blue


def rgb8_of(integer):
    """The `(r, g, b)` of a stored colour integer."""
    return ((integer >> 16) & 0xFF, (integer >> 8) & 0xFF, integer & 0xFF)


def _solid_colour(body):
    temperature = getattr(getattr(body, "star", None), "temperature", None)
    if not temperature or temperature <= 0.0:
        return None
    column = air_column.column_of(body)
    if column is None:
        return None
    stored = aerosols.of_body(body)
    haze, dust, cloud = aerosols.colour_layers(stored, _CENTRES)
    result = atmosphere_colour.atmosphere_colour(
        column, temperature_k=temperature, ozone_du=getattr(body, "ozone_column_du", None) or 0.0,
        ground_albedo=atmosphere_colour.STAND_IN_GROUND_ALBEDO, haze=haze, dust=dust, cloud=cloud, cloud_cover=stored.cloud_cover)
    return result.colour.rgb8


def colour_of(body):
    """The `(r, g, b)` of a planet or moon's atmosphere colour, or `None` for a body that has none (no air, no star temperature, no mixture)."""
    if getattr(body, "body_type", "t") == "g":
        found = giant_colour.body_giant_colour(body)
        return None if found is None else found.atmosphere.colour.rgb8
    return _solid_colour(body)


def generate(body):
    """Sets `BODY_FIELDS` on a planet or moon (call after its ozone column and aerosols). A body with no colour stores `None`."""
    rgb8 = colour_of(body)
    body.atmosphere_colour = None if rgb8 is None else rgb_integer(rgb8)
