# planetgen/physics/cloud_infrared.py

"""
The infrared opacity of the sulfuric acid cloud (GEN.1020; the cloud of `physics/sulfuric_cloud.py`, the infrared model of `physics/line_opacity.py`)
=======================================================================================================================================

The acid cloud is stated by its depth at 550 nm and its droplet radius. This module turns that into an optical depth on each of the 130
infrared bands and puts it in the layers of an air column.

* **Optics per band** (`physics/data/acid_cloud_optics.json`, made by `scripts/build_acid_cloud_optics.py`): Mie spheres (`physics/mie.py`) of
  76 and 81 percent sulfuric acid at 298 K (the HITRAN aerosol indices of Myhre and others 2003), in a gamma size distribution of effective
  variance 0.2, for effective radii of 0.3 to 3.85 micron. Stored per band: extinction over the 550 nm extinction of the same droplets, the
  single-scattering albedo and the asymmetry factor.
* **Scattering in a non-scattering solver.** The column solver absorbs and emits only, so the scattered part is scaled by the similarity
  relation: the depth the solver sees is `tau_ext (1 - omega g)`, the absorbed part plus the scattering that is not forward.
* **Where it sits.** The cloud's mass is spread evenly in pressure over the layers between its base and its top. The base is the layer
  where the air is cooler than `sulfuric_cloud.base_temperature_k` going up; the top is `top_over_base_pressure` of the base pressure.

Limits: the indices are for 298 K droplets (the cloud is 250 to 370 K; the temperature dependence of the acid's infrared bands is small
compared with the strength spread of 76 to 81 percent); the indices are held at their 400 and 7500 cm-1 values outside that range; the scattering is
scaled, not solved. Registry rows CF1 to CF5.
"""

import json
import os

import numpy as np

from planetgen.physics import elements, sulfuric_cloud

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "acid_cloud_optics.json")

TOP_OVER_BASE_PRESSURE = 0.0267 / 1.391
"""float: The pressure of the cloud top over that of its base: [R0157: Table 1], the upper cloud's 70 km level
(0.0267 bar) over the lower cloud's 47.5 km level (1.391 bar). The same ratio is used for any air (row CF3)."""
DEFAULT_WEIGHT_PERCENT = 80.0
"""float: The acid strength the optics are read at, between the 76 and 81 percent tables: Gillmann and others 2024 give 80 percent at the equator."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        return json.load(handle)


DOCUMENT = _load()
RADII_UM = np.array(DOCUMENT["radii_um"], dtype=float)
PERCENTS = [int(p) for p in DOCUMENT["weight_percent"]]


def _table(percent, key):
    return np.array([row[key] for row in DOCUMENT["optics"][str(percent)]], dtype=float)


def optics(radius_um, weight_percent=DEFAULT_WEIGHT_PERCENT):
    """
    The per-band optics of the droplets: `(extinction over visible extinction, single-scattering albedo, asymmetry factor)`, three arrays of 130,
    interpolated in ln radius (held at the ends of the radii) and linearly in acid strength between the two tables.
    """
    position = np.interp(np.log(radius_um), np.log(RADII_UM), np.arange(len(RADII_UM)))
    low, high = int(np.floor(position)), int(np.ceil(position))
    share = position - low
    strength = np.clip((weight_percent - PERCENTS[0]) / (PERCENTS[1] - PERCENTS[0]), 0.0, 1.0)
    result = []
    for key in ("ext_over_visible", "omega", "g"):
        by_strength = []
        for percent in PERCENTS:
            table = _table(percent, key)
            by_strength.append((1.0 - share) * table[low] + share * table[high])
        result.append((1.0 - strength) * by_strength[0] + strength * by_strength[1])
    return tuple(result)


def effective_depth(tau_550, radius_um, weight_percent=DEFAULT_WEIGHT_PERCENT):
    """The depth on each band that the non-scattering solver sees for a cloud of depth `tau_550` at 550 nm: `tau (ext/visible) (1 - omega g)`."""
    ratio, omega, g = optics(radius_um, weight_percent)
    return tau_550 * ratio * (1.0 - omega * g)


def layer_depth(layer_pressure_pa, layer_mass_kg_m2, layer_temperature_k, tau_550, radius_um, weight_percent=DEFAULT_WEIGHT_PERCENT):
    """
    The cloud's depth in each layer of an air column, shape (layers, bands): zero without a cloud (`tau_550` of 0, or no layer that cool).

    Args:
        layer_pressure_pa, layer_mass_kg_m2, layer_temperature_k (ndarray): Each layer's pressure, mass per area and temperature (top first).
    """
    pressure = np.asarray(layer_pressure_pa, dtype=float)
    mass = np.asarray(layer_mass_kg_m2, dtype=float)
    temperature = np.asarray(layer_temperature_k, dtype=float)
    depth = np.zeros((len(pressure), DOCUMENT["bands"]))
    if tau_550 <= 0.0:
        return depth
    warm = temperature > sulfuric_cloud.number("base_temperature_k")
    first_warm = int(np.argmax(warm)) if warm.any() else len(pressure)
    if first_warm == 0:
        return depth
    base_pressure = pressure[first_warm - 1]
    inside = np.zeros(len(pressure), dtype=bool)
    inside[:first_warm] = pressure[:first_warm] >= base_pressure * TOP_OVER_BASE_PRESSURE
    share = mass * inside
    share = share / share.sum()
    return share[:, None] * effective_depth(tau_550, radius_um, weight_percent)[None, :]
