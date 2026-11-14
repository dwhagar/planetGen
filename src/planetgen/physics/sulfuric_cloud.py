# planetgen/physics/sulfuric_cloud.py

"""
The sulfuric acid cloud of an air (GEN.1019, docs/design/composition/line-opacity.md; part of the cloud of `physics/aerosols.py`)
==================================================================================================================================

The condensate clouds of `aerosols` form where a gas reaches its saturation pressure. Sulfuric acid does not: its vapor is made in the
air from sulfur dioxide, an oxidiser and water, and the droplets stand where the air is cool enough, so the cloud has a supply, a base
and a mass limit but no saturation crossing.

* **Supply.** The gates of `acid_medium` (sulfur dioxide, an oxidiser, a trace of water): an air without sulfur gets no acid cloud.
* **Base.** The warmest level going up the dry adiabat (`aerosols.profile`) where sulfuric acid condenses out of the air's own acid vapor and
  water (`sulfur_cloud_set.acid_cloud_of_column`, the GEN.1424 curves of `sulfuric_vapor`), so the base moves with the air's sulfur and water (a drier
  air condenses cooler, a sulfur-richer one warmer); the ground where the ground already is. The same window gives the acid's strength. An air whose
  acid never condenses has no cloud.
* **Limit.** All the sulfur dioxide above the base made into acid, `m = x_SO2 P_base M_acid / (M g)` kg per m2. The depth of a cloud of
  that mass in droplets of radius `r` is `3 Q m / (4 rho r)` (`aerosols.CLOUD_Q_EXT`).
* **Draw.** The depth is uniform on the published Venus total (20 to 30), the droplet radius log-uniform on the published mode radii
  (0.3 to 3.85 micron), the depth held to the limit. The cloud covers the sky (`cover`). It is grey, as the other clouds are.

The numbers are in `physics/data/sulfuric_cloud.json` with their tier and tweak level; the registry rows are SC1 to SC5 of the gap plan.
"""

import math
import os


from planetgen.physics import acid_medium, aerosols, air_column, elements, sulfur_cloud_set

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "sulfuric_cloud.json")
ACID_MOLAR_MASS_KG_MOL = 0.098079
"""float: The molar mass of sulfuric acid, kg/mol (H2SO4: 2 x 1.008 + 32.06 + 4 x 15.999)."""


class SulfuricCloudError(ValueError):
    """A rule file that is refused."""


def _load():
    """The rule file's document, checked (`elements.load_rule_numbers`): every number has a value, a tier, a tweak and a basis."""
    return elements.load_rule_numbers(DATA_FILE, SulfuricCloudError)


DOCUMENT = _load()
NUMBERS = DOCUMENT["numbers"]
VERSION = DOCUMENT["version"]


def number(name):
    """The value of the rule file's number `name`."""
    return NUMBERS[name]["value"]


def supplied(column):
    """Whether the air has the sulfur dioxide, the oxidiser and the water to make the acid (the gates of `acid_medium`)."""
    return column.fraction("SO2") >= acid_medium.number("so2_min_share") \
        and column.fraction("O2") + column.fraction("CO2") >= acid_medium.number("oxidant_min_share") \
        and column.fraction("H2O") >= acid_medium.number("water_min_share")


def skin_k(energy_flux_w_m2):
    """The temperature the air turns isothermal at for an absorbed `energy_flux_w_m2` (`acid_medium.skin_k`), or `None` for no flux. A greenhouse
    world's ground is far warmer than the grey skin factor of `aerosols.SKIN_FACTOR` would put the cold top, and the cloud stands near the top."""
    return acid_medium.skin_k({"energy_flux_w_m2": energy_flux_w_m2}, None)


def window(column, temperature_k, skin=None):
    """The `sulfur_cloud_set.Cloud` (base temperature, base and top pressure, strength) of the air's acid cloud, or `None` when the air has no acid
    cloud: it lacks the sulfur, an oxidiser or the water (`supplied`), or its acid does not condense anywhere in its profile (GEN.1424)."""
    if not supplied(column):
        return None
    return sulfur_cloud_set.acid_cloud_of_column(column, temperature_k, skin)


def base_pressure_pa(column, temperature_k, skin=None):
    """The pressure of the cloud base, Pa, or `None` when the air has none. `skin` is the temperature the air turns isothermal at (`skin_k`), else
    the grey factor of the ground temperature."""
    found = window(column, temperature_k, skin)
    return None if found is None else found.base_pressure_pa


def acid_column_limit_kg_m2(column, temperature_k, skin=None):
    """The most acid the sulfur above the cloud base can make, kg per m2 of ground; 0 where there is no cloud."""
    base = base_pressure_pa(column, temperature_k, skin)
    if base is None:
        return 0.0
    return column.fraction("SO2") * base / column.gravity_m_s2 * ACID_MOLAR_MASS_KG_MOL / column.mean_molar_mass_kg_mol


def optical_depth(mass_kg_m2, radius_um):
    """The in-cloud depth at 550 nm of `mass_kg_m2` of acid in droplets of `radius_um`."""
    return 3.0 * aerosols.CLOUD_Q_EXT * mass_kg_m2 / (4.0 * number("density_kg_m3") * radius_um * 1.0e-6)


GOLDEN = sulfur_cloud_set.GOLDEN_CENTER
"""float: The multiplier that makes the depth draw from the one stored `aerosol_draw` (fractional part of `aerosol_draw x GOLDEN`), so depth and
radius are not the same number."""


def radius_um(radius_draw):
    """The droplet effective radius, micron, for a draw in [0, 1): log-uniform on the published mode radii."""
    low, high = number("radius_low_um"), number("radius_high_um")
    return math.exp(math.log(low) + radius_draw * (math.log(high) - math.log(low)))


def draw(column, temperature_k, tau_draw, radius_draw, skin=None):
    """
    `(tau_cloud, cloud_cover)` of the acid cloud for two draws in [0, 1): the depth within the published bracket and the droplet radius within
    the published modes, the depth held to what the sulfur can make. `(0, 0)` where there is no acid cloud.
    """
    limit = acid_column_limit_kg_m2(column, temperature_k, skin)
    if limit <= 0.0:
        return 0.0, 0.0
    low, high = number("optical_depth_low"), number("optical_depth_high")
    tau = min(low + tau_draw * (high - low), optical_depth(limit, radius_um(radius_draw)))
    return tau, number("cover")


def of_draw(column, temperature_k, aerosol_draw, skin=None):
    """
    `(tau_550, cloud_cover, radius_um)` of the acid cloud of a body from its stored `aerosol_draw` (the depth from its fractional part times
    `GOLDEN`, the radius from the draw itself), so the stored aerosol values repeat the cloud with no other column. `(0, 0, radius)` without a cloud.
    """
    radius = radius_um(aerosol_draw)
    tau, cover = draw(column, temperature_k, (aerosol_draw * GOLDEN) % 1.0, aerosol_draw, skin)
    return tau, cover, radius


def of_body(body):
    """`(tau_550, cloud_cover, radius_um)` of the acid cloud of a planet or moon from its stored values (its air, surface temperature, greenhouse
    factor and `aerosol_draw`); `(0, 0, 0)` for a body with no air or no stored aerosols."""
    temperature = getattr(body, "surface_temperature", None)
    draw_value = getattr(body, "aerosol_draw", None)
    column = air_column.column_of(body)
    if column is None or not temperature or draw_value is None:
        return 0.0, 0.0, 0.0
    return of_draw(column, temperature, draw_value, aerosols.skin_of(body))


def infrared_cloud(column, temperature_k, aerosol_draw, skin=None):
    """
    `(tau_550, radius_um, base_temperature_k, weight_percent)` for `radiative_convective.solve(cloud=...)`: the acid cloud's depth and droplet radius
    from the stored `aerosol_draw` (`of_draw`) and the base temperature and strength its air's acid fixes (`window`); `None` without an acid cloud.
    """
    found = window(column, temperature_k, skin)
    if found is None:
        return None
    tau, _cover, radius = of_draw(column, temperature_k, aerosol_draw, skin)
    if tau <= 0.0:
        return None
    return tau, radius, found.base_temperature_k, found.weight_percent
