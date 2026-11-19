# planetgen/physics/aerosol_material.py

"""
What the cloud and the haze of a body are made of (MAP.418)
===========================================================

A body stores the optical depths of its cloud, haze and dust (`tau_cloud`, `tau_haze`, `tau_dust`, `cloud_cover`; DB.68) but not what the
cloud and haze are made of. The material is not a new stored value: it follows from the air and the temperature the body stores, by the same
rules that made the depths, so the painter can colour each layer by its material (water white, sulfuric acid pale yellow, carbon dioxide ice
grey-white, ammonia ice and methane ice, dust tan, tholin orange) without a column of its own.

* **Cloud.** `tau_cloud` is the condensate cloud (`aerosols.draw_cloud`, the gases of `aerosols.CONDENSATE_DENSITY_KG_M3` that reach saturation in
  the air's profile) plus the sulfuric acid cloud (`sulfuric_cloud.of_body`, from the stored `aerosol_draw`). The acid share is recomputed from the
  stored draw; what is left of `tau_cloud` is the condensate. The condensate that makes most of it is the gas with the largest mass over its
  particle density (the depth of a cloud of one droplet radius goes as `3 Q m / (4 rho r)`, so for the one radius the cloud draw uses, the order is
  the order of `m / rho`). Where the acid's depth is over the rest it is the cloud's material.
* **Haze.** `tau_haze` is the hydrocarbon haze of the air's methane (`aerosols.haze_in`, tholin); `tau_dust` is the dust the wind lifts off a
  dry ground. The layer is the mix of the two by their depths.

Reuse: `aerosols.condensate_columns_kg_m2`, `aerosols.CONDENSATE_DENSITY_KG_M3`, `aerosols.skin_of`, `sulfuric_cloud.of_body`, `air_column.column_of`.
"""

from planetgen.physics import aerosols, air_column, sulfuric_cloud

ACID = "H2SO4"
"""str: The material name of the sulfuric acid cloud (the condensates are named by their gas formula)."""

THOLIN = "tholin"
"""str: The material name of the haze of methane's photochemistry (`aerosols.haze_in`)."""

DUST = "dust"
"""str: The material name of the wind-lifted dust of a dry ground (`aerosols.draw_dust`)."""


def cloud_material(body):
    """
    The material that makes most of the cloud of `body`: a gas formula of `aerosols.CONDENSATE_DENSITY_KG_M3` (H2O, NH3, CH4, CO2) or `ACID`. `None`
    for a body with no cloud stored, no air or no surface temperature.
    """
    depth = getattr(body, "tau_cloud", None)
    temperature = getattr(body, "surface_temperature", None)
    if not depth or depth <= 0.0 or not temperature:
        return None
    column = air_column.column_of(body)
    if column is None:
        return None
    acid = sulfuric_cloud.of_body(body)[0]
    rest = max(depth - acid, 0.0)
    if acid > rest:
        return ACID
    weights = {formula: mass / aerosols.CONDENSATE_DENSITY_KG_M3[formula] for formula, mass in aerosols.condensate_columns_kg_m2(column, temperature).items()}
    if not weights:
        return ACID if acid > 0.0 else None
    return max(weights, key=weights.get)


def haze_mix(body):
    """`{material: optical depth}` of the haze layer of `body` (tholin from `tau_haze`, dust from `tau_dust`), the layers that are over 0 only."""
    mix = {}
    for material, field in ((THOLIN, "tau_haze"), (DUST, "tau_dust")):
        value = getattr(body, field, None)
        if value is not None and value > 0.0:
            mix[material] = value
    return mix
