# planetgen/physics/habitability_explain.py

"""
The composition behind a body's radiation rating (GEN.445, docs/design/composition/radiation-model.md section 9)
=============================================================================================================

The Radiation row of the habitability rating is a number and a verdict. This names the terms it is made of and the composition
that sets each one, from the body's stored record, so a reader can see why a world is Blue or Red and what would move it:

* cosmic rays (`dose_gcr_msv_yr`, times the heliosphere's multiplier): the shielding column `P / g` of the air and the
  magnetic cutoff of the dipole (`radiation.cutoff_rigidity_gv`);
* stellar particles (`dose_sep_msv_yr`): the star's flares, cut by the same column and field;
* trapped particles (`dose_belt_msv_yr`): a moon in the radiation belts of its parent planet's field (GEN.983);
* the ground (`dose_ground_msv_yr`): the gamma rays of the potassium, thorium and uranium the body stores (`radiogenic`) and the
  radon of its uranium;
* the ultraviolet (`uv_surface_index`, a multiple of Earth's DNA-damaging dose): the ozone column of its oxygen (`ozone_column_du`)
  and the rest of its air, under the star's own light.

`limit` says which of the dose and the ultraviolet sets the Radiation rating (the smaller of the two scores, as
`habitability.radiation_score` takes it), `tier` the rating's colour. Nothing is drawn or stored here.

Use:

    terms = habitability_explain.radiation_terms(planet)
    terms.gcr_msv_yr, terms.column_g_cm2, terms.air, terms.cutoff_gv, terms.potassium_ppm, terms.uv_index, terms.limit
"""

from dataclasses import dataclass

from planetgen.physics import (
    atmosphere,
    habitability,
    habitability_world,
    radiation,
    radiogenic,
)

AIR_SHARE_MIN = 0.01
"""float: A gas below this share of the air is left out of the list that names what shields (it moves the column by less than that)."""


@dataclass(frozen=True)
class RadiationTerms:
    """
    The radiation terms of one body.

    Attributes:
        total_msv_yr (float): The stored surface dose, mSv/yr.
        gcr_msv_yr (float): Cosmic rays with the heliosphere's multiplier.
        sep_msv_yr (float): Stellar particles.
        ground_msv_yr (float): The ground's gamma and radon.
        belt_msv_yr (float): The trapped particles of the parent planet's belts (a moon only; 0 for the rest).
        helio_mult (float): The compressed-heliosphere multiplier on the cosmic rays (1 in open space).
        column_g_cm2 (float): The air above the ground, `P / g`, g/cm2 (0 for an airless body).
        air (tuple): `(name, share)` of the gases of at least `AIR_SHARE_MIN` of the pressure, largest first, named as the
            gas table names them ("nitrogen").
        cutoff_gv (float): The dipole's equatorial cutoff rigidity, GV (0 without a field).
        potassium_ppm, thorium_ppm, uranium_ppm (float or None): The stored radiogenic content, ppm of the bulk silicate at
            formation; `None` for a body saved without it.
        uv_index (float): The ground DNA-damaging dose as a multiple of Earth's.
        ozone_du (float or None): The ozone column, Dobson units; `None` when none is stored.
        limit (str): `"dose"` or `"ultraviolet"`: which score is the smaller.
        tier (str): The Radiation rating's colour.
    """

    total_msv_yr: float
    gcr_msv_yr: float
    sep_msv_yr: float
    ground_msv_yr: float
    belt_msv_yr: float
    helio_mult: float
    column_g_cm2: float
    air: tuple
    cutoff_gv: float
    potassium_ppm: float | None
    thorium_ppm: float | None
    uranium_ppm: float | None
    uv_index: float
    ozone_du: float | None
    limit: str
    tier: str


def air_shares(body):
    """`((name, share), ...)` of the gases in the air of `body` that make up at least `AIR_SHARE_MIN` of its pressure,
    largest first."""
    pressures = atmosphere.gas_partials_pa(body)
    total = sum(pressures.values())
    if total <= 0.0:
        return ()
    shares = [(atmosphere.gas_name(formula), value / total) for formula, value in pressures.items()]
    return tuple(sorted(((name, share) for name, share in shares if share >= AIR_SHARE_MIN), key=lambda item: -item[1]))


def radiation_terms(body):
    """
    The `RadiationTerms` of a scored rocky body.

    Raises:
        ValueError: The body is a gas giant or has no stored dose (it is not scored).
    """
    values = habitability_world.values_from_body(body)
    if values.get("body_type") == "g" or values.get("surface_dose_msv_yr") is None:
        raise ValueError("a body with no stored dose has no radiation terms")
    world = habitability_world.world_from(values)
    helio = values.get("dose_helio_mult") or 1.0
    uv = habitability.uv_score(world)
    dose = habitability.band_score(world.dose_sv_yr, habitability.HUMAN.dose, log=True)
    content = radiogenic.content_of(body)
    return RadiationTerms(
        total_msv_yr=values["surface_dose_msv_yr"], gcr_msv_yr=(values.get("dose_gcr_msv_yr") or 0.0) * helio,
        sep_msv_yr=values.get("dose_sep_msv_yr") or 0.0, ground_msv_yr=values.get("dose_ground_msv_yr") or 0.0,
        belt_msv_yr=values.get("dose_belt_msv_yr") or 0.0, helio_mult=helio,
        column_g_cm2=radiation.column_g_cm2(values.get("atmospheric_pressure_pa") or 0.0, values.get("gravity_g")),
        air=air_shares(body), cutoff_gv=habitability_world.cutoff_gv(values),
        potassium_ppm=None if content is None else content.k0_ppm, thorium_ppm=None if content is None else content.th0_ppm,
        uranium_ppm=None if content is None else content.u0_ppm, uv_index=values.get("uv_surface_index") or 0.0,
        ozone_du=getattr(body, "ozone_column_du", None),
        limit="ultraviolet" if uv is not None and uv < dose else "dose",
        tier=habitability.tier(habitability.radiation_score(world)))
