# planetgen/galaxy/profile_report.py

"""
A Sector Profile as plain data (API.27, UX.113): numbers, strings, lists and dicts that serialize to JSON as they
stand, so `GET /api/galaxy/profile` and the sector page show the same fields. Reads a `SectorProfile`
(`galaxy/profile.py`) and changes nothing.
"""

from planetgen.util import detmath


def _bubble(bubble, position_pc):
    if bubble is None:
        return None
    return {"radius_pc": bubble.radius_pc, "ambient_density_cm3": bubble.ambient_density_cm3,
            "temperature_k": bubble.temperature_k, "center_pc": list(bubble.center_pc),
            "centre_distance_pc": detmath.dist(bubble.center_pc, position_pc)}


def profile_report(profile):
    """
    The profile of a place as a dict.

    Returns:
        dict: `shell_index`, `band_index`, `wedge_index` (`null` for a profile of a point), `center_pc`, `edge_pc`,
            `volume_pc3`, `region`, `density`, `expected_systems`, `mix` (population -> density), `populations`
            (population -> `density`, `share`, `age_gy` as `[youngest, oldest]`, `fe_h` as a list of `[weight, mean,
            sigma]`, `dispersion_kms`, `rotation`), `gas` (`cold_fraction`, `ambient`, `ambient_density_cm3`,
            `mass_solar`, `superbubble` as `null` or the cavity's `radius_pc`, `ambient_density_cm3`,
            `temperature_k`, `center_pc` and `centre_distance_pc`), `dark_matter` (`density_solar_pc3`, `mass_solar`)
            and `kind_factors` (kind -> factor).
    """
    shell, band, wedge = profile.address if profile.address is not None else (None, None, None)
    gas = profile.gas
    return {
        "shell_index": shell, "band_index": band, "wedge_index": wedge,
        "center_pc": list(profile.position_pc), "edge_pc": profile.edge_pc, "volume_pc3": profile.volume_pc3,
        "region": profile.region, "density": profile.density, "expected_systems": profile.expected_systems,
        "mix": dict(profile.mix),
        "populations": {
            name: {"density": pop.density, "share": pop.share, "age_gy": list(pop.age_gy),
                   "fe_h": [list(peak) for peak in pop.fe_h], "dispersion_kms": pop.dispersion_kms,
                   "rotation": pop.rotation}
            for name, pop in profile.populations.items()},
        "gas": {"cold_fraction": gas.cold_fraction, "ambient": gas.ambient,
                "ambient_density_cm3": gas.ambient_density_cm3, "mass_solar": gas.mass_solar,
                "superbubble": _bubble(gas.superbubble, profile.position_pc)},
        "dark_matter": {"density_solar_pc3": profile.dark_matter.density_solar_pc3,
                        "mass_solar": profile.dark_matter.mass_solar},
        "kind_factors": dict(profile.kind_factors),
    }
