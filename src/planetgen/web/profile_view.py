# planetgen/web/profile_view.py

"""
The Sector Profile as the sector pages show it (UX.113): `profile_view(profile)` turns the answer of
`GET /api/galaxy/profile` (`galaxy/profile_report.py`) into the words and numbers of `partials/profile.html`.
"""

from planetgen import tuning
from planetgen.generation.run_galaxy import LARGE_RUN_WARNING_THRESHOLD
from planetgen.generation.stats import format_bytes, format_duration
from planetgen.physics.quantity import Quantity
from planetgen.web.lib.fmt import format_number, format_speed_kms

HEAVY_CELL_SECONDS = 120.0
"""float: A cell whose fill is expected to take this long on one worker gets the heavy-cell warning before it is asked."""

FAR_RUN_SECTORS = LARGE_RUN_WARNING_THRESHOLD
"""int: A neighbourhood of more sectors than this is a very large run, which the sector page says before it is confirmed."""

SHOWN_SHARE = 0.0005
"""float: A population holding less than this share of a cell is left out of its table."""

def _age_label(age_gy):
    """An age in Gy as the label of an open cluster's population names it: Myr below a Gy, three significant figures."""
    if age_gy >= 1.0:
        return f"{age_gy:.3g} Gy"
    return f"{Quantity('duration', age_gy, 'Gyr').Myr:.3g} Myr"


REGION_LABELS = {
    "nuclear": "Galactic nucleus", "bulge": "Bulge", "disk": "Disk", "halo": "Halo", "cluster": "Star cluster",
    "satellite": "Satellite galaxy",
}
POPULATION_LABELS = {
    "young": "Young disk", "intermediate": "Intermediate disk", "old": "Old disk", "bulge": "Bulge", "halo": "Halo",
    "cluster": "Cluster", "satellite": "Old satellite stars", "nuclear_cluster": "Nuclear cluster",
    "nuclear_disc": "Nuclear disc", "satellite_young": "Young satellite stars",
    "satellite_intermediate": "Intermediate satellite stars",
}
POPULATION_LABELS.update({
    name: "Open cluster stars, {low} to {high}".format(
        low=_age_label(tuning.STELLAR_POPULATION_AGE_RANGES_GY[name][0]),
        high=_age_label(tuning.STELLAR_POPULATION_AGE_RANGES_GY[name][1]))
    for name in tuning.OPEN_CLUSTER_POPULATIONS})
ROTATION_LABELS = {"disk": "Follows the disk", "bar": "Follows the bar", "orbit": "Orbits the galaxy", "none": "No net motion"}
AMBIENT_LABELS = {
    "disk": "Disk gas", "superbubble": "Superbubble cavity (hot, ionised gas)", "hot-halo": "Hot halo gas",
    "intergalactic": "Intergalactic medium",
}


def _sig(value):
    """Three significant figures below a thousand, a whole number to six digits, then scientific (UX.20)."""
    return format_number(value) if abs(value) >= 1000 else format_number(value, ".3g")


def _fe_h_text(peaks):
    """`"-0.10 ± 0.20"`, with `" + "` between the peaks of a place two components reach."""
    return " + ".join(f"{mean:+.2f} ± {sigma:.2f}" for _weight, mean, sigma in peaks).replace("-", "−")


def _age_text(age_gy):
    low, high = age_gy
    return f"{low:g} to {high:g}"


def profile_view(profile, charted=False):
    """
    The words and numbers for `partials/profile.html`.

    Args:
        profile (dict): `GET /api/galaxy/profile`'s answer.
        charted (bool): The cell has its sector already, so there is nothing left to warn about before running it.

    Returns:
        dict: `region`, `expected_systems`, `density`, `rows` (one per population above `SHOWN_SHARE`, most of the
            cell first: `label`, `share_percent`, `age`, `fe_h`, `dispersion`, `rotation`), `type_ia_percent` (the share of the place's supernova remnants that are Type Ia), `gas` (`cold_percent`,
            `ambient`, `ambient_density_cm3`, `mass_solar`, `superbubble`), `dark_matter` (`density_solar_pc3`,
            `mass_solar`), `heavy` (`None`, or the warning's figures: `systems`, `time`, `size`, `measured`) and
            `unresolved`.
    """
    rows = sorted(
        ({"label": POPULATION_LABELS.get(name, name), "share_percent": format_number(100.0 * pop["share"], ",.1f"),
          "age": _age_text(pop["age_gy"]), "fe_h": _fe_h_text(pop["fe_h"]), "dispersion": format_speed_kms(pop["dispersion_kms"]),
          "rotation": ROTATION_LABELS.get(pop["rotation"], pop["rotation"])}
         for name, pop in profile["populations"].items() if pop["share"] >= SHOWN_SHARE),
        key=lambda row: -float(row["share_percent"].replace(",", "")))
    gas = profile["gas"]
    estimate = profile.get("estimate") or {}
    heavy = None
    if not charted and (profile.get("unresolved") or estimate.get("seconds", 0.0) >= HEAVY_CELL_SECONDS):
        heavy = {"systems": estimate.get("systems"), "time": format_duration(estimate.get("seconds", 0.0)),
                 "size": format_bytes(estimate.get("bytes", 0)), "measured": estimate.get("measured", False),
                 "slow": estimate.get("seconds", 0.0) >= HEAVY_CELL_SECONDS}
    dark = profile["dark_matter"]
    return {
        "region": REGION_LABELS.get(profile["region"], profile["region"]),
        "expected_systems": _sig(profile["expected_systems"]), "density": _sig(profile["density"]), "rows": rows,
        "gas": {"cold_percent": format_number(100.0 * gas["cold_fraction"]),
                "ambient": AMBIENT_LABELS.get(gas["ambient"], gas["ambient"]),
                "ambient_density_cm3": _sig(gas["ambient_density_cm3"]), "mass_solar": _sig(gas["mass_solar"]),
                "superbubble": gas["superbubble"]},
        "dark_matter": {"density_solar_pc3": _sig(dark["density_solar_pc3"]), "mass_solar": _sig(dark["mass_solar"])},
        "type_ia_percent": format_number(100.0 * profile.get("type_ia_share", 0.0)),
        "heavy": heavy, "unresolved": profile.get("unresolved"),
    }
