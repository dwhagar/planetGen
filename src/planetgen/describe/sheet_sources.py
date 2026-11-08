# planetgen/describe/sheet_sources.py

"""
The one-line source of a physics row of the property sheet (UX.187): how the value was worked out,
so a reader of a page or of `GET /api/objects/<ref>/sheet` sees where a number comes from beside the number. Every line is
one sentence, written from the docstring of the module that makes the value, and a row with no line here has none: nothing
is invented for a column whose rule is not written down. `tests/test_sheet_sources.py` holds the table to the sheet's rows.
"""

_LID = ("a stagnant lid (one rigid plate) or a mobile one (plate tectonics, which needs a liquid surface ocean), drawn once with "
        "the interior")

_STEPS = ("Five hard steps, each an exponential waiting time that starts when the one before is done, faster on a warm world "
          "[R4642]")
_AGES = _STEPS + "."

_GROUND = ("The areal share of each ground type: what the hydrosphere covers, then the dominant crust by the heat flow and the mantle "
           "[R2024], then the modifiers the body's own water, air, tectonics and tides allow.")

_GROUND_TYPES = ("ocean", "ice", "ultramafic", "anorthositic", "basaltic", "felsic", "clay", "oxidized", "metal", "carbon", "evaporite",
                 "sulfur", "intermediate")

_BODY = {
    "eccentricity": "A planet draws from the Kepler multi-planet law [R1307] and is circularised by its tide; a regular "
                    "moon keeps a small forced value, a captured moon a thermal one.",
    "flux_variation": "The swing of the starlight between periapsis and apoapsis, peak to peak over the mean, from the eccentricity.",
    "tidal_flux_w_m2": "[R3464] tidal power from the stored eccentricity and the body's rheology.",
    "internal_heat_flux_w_m2": "Radiogenic plus primordial heat in one-box parameterised mantle convection ([R4339]; "
                               "[R1716]).",
    "urey_ratio": "The radiogenic share of the interior heat, worked out from the mantle and not a draw.",
    "lid_regime": "One of " + _LID + ".",
    "c_th0_ppm": "Follows the host star's thorium (its Eu relation, truncated at plus or minus 0.5 dex), not a draw of its own.",
    "c_u0_ppm": "Worked out from the thorium by the r-process production ratio and the decay of the galaxy's earlier ejecta.",
    "rotation_period_hours": "Birth spin worn down by the tide of what the body orbits, to a free, locked or resonant state "
                             "[R1584].",
    "axial_tilt_deg": "Isotropic for a rocky planet, Rayleigh about 15 degrees for a gas giant or an untidied moon, and on the orbit "
                      "normal for a locked body.",
    "magnetic_moment_a_m2": "The dipole of the dynamo that runs in the body: a rocky body's liquid core (its thermal history) or a giant's metallic or "
                            "ionic layer (the heat it radiates), by the buoyancy-flux scaling and the spin.",
    "dipole_tilt_deg": "A dipole-dominated dynamo keeps its axis near the spin axis (Rayleigh, scale 6 degrees; Earth 11, Jupiter 10); a multipolar one or an ice giant's thin "
                       "shell has it at any angle, isotropic (Uranus 59, Neptune 47) [R0784, R4196].",
    "dipole_offset_class": "Offset (over 0.15 radii from the center) for an ice giant, a multipolar dynamo and the iron-rich small body; centered otherwise [R4196].",
    "core_light_fraction": "A drawn 3 to 12 percent of the metal, split by the rock's FeO and the pressure at the base of the mantle: a reduced body puts Si and C in its core, an oxidised one S, O and H [R4608, R4833].",
    "ring_outer_km": "Ends at the Roche limit, 2.456 R (planet density over particle density)^(1/3), and inside the innermost large moon.",
    "ring_age_yr": "The years since the moon that fell inside the Roche limit, or the comet that passed inside it, made the ring; the ring thins by 1/e every half billion years.",
    "ring_optical_depth": "An output of the drawn mass, the extent and the particle size, 3 M / (4 density r area).",
    "ground_type": _GROUND,
    "ground_grain": "A surface under 0.1 Gyr old or molten is not yet regolith; an older one is a slab, crushed or powder by its grain.",
    **{f"ground_{name}_fraction": _GROUND for name in _GROUND_TYPES},
    "life_stage": "The latest stage whose age is not past the age of the star.",
    "life_age_abiogenesis_gyr": _AGES,
    "life_age_photosynthesis_gyr": _AGES,
    "life_age_complex_cells_gyr": _STEPS + "; the wait for complex cells starts only when the air holds free oxygen.",
    "life_age_multicellularity_gyr": _AGES,
    "life_age_technological_civilization_gyr": _AGES,
}
"""dict: A planet's or moon's column key -> its source line."""

_PLANET = {
    **_BODY,
    "orbital_inclination_deg": "A Rayleigh law about the star's equator, its dispersion from the Kepler multiples [R1311].",
}
"""dict: A planet's column key -> its source line (a moon's inclination has no line written)."""

_STAR = {
    "rotation_period_hours": "A white dwarf log-normal about a day, a giant a slow 1 to 10 km/s, a cool dwarf from "
                             "gyrochronology [R1429], a hot star from a drawn speed; a scattered star keeps its own draw.",
    "axial_tilt_deg": "Rayleigh about 15 degrees from the star's reference normal.",
    "surface_dipole_field_gauss": "A cool dwarf follows the large-scale field against the Rossby number, B ~ Ro^-1.4 [R4544, R5144]; an Ap dwarf has the polar field of its kind.",
    "nuclear_arrival_age_gy": "The age MIST's tracks put the star on the main sequence [R0773]; the star counts as pre-main-sequence "
                              "only until its radius has settled, the earlier age BHAC15 gives [R0254].",
    "disk_kind": "A young star holds a class 0, I or II disk by its age; a main-sequence star has a cold debris disk when "
                 "the dust of its outer belt is bright enough for a survey to find (about 1e-5 of the star's light [R4081], the dust law of [R4785]).",
    "disk_inner_au": "A young star's disk starts at the dust sublimation radius; a debris disk starts at the inner edge of "
                     "the star's outer belt.",
    "disk_outer_au": "A young star's disk ends log-uniform between 10 and 300 AU; a debris disk ends at the outer edge of the "
                     "star's outer belt.",
    "chemistry_branch": "A star older than the thick disk's start is alpha-rich with [R0349]'s logistic probability of its "
                        "[Fe/H] (age gate, [R1863]); one draw of the star's own decides it.",
    "chemistry_alpha_fe": "The mean [X/Fe] of O, Mg, Si and Ca: the branch's hinge fit at the star's [Fe/H] [R0349, R3294] "
                          "plus the star's residuals [R0535].",
    "chemistry_metal_fraction": "The mass fraction of everything heavier than helium, from [Fe/H] and the star's [X/Fe] on the proto-solar "
                                "mix of [R0168]; the composition the star formed with.",
    "chemistry_dredge_up_helium": "[R0467]: the rise in helium at the surface of a red giant of this mass and metal fraction.",
}
"""dict: A star's column key (or the key of a row worked out from the stored values, UX.254) -> its source line."""

_MOON = {key: line for key, line in _BODY.items() if not key.startswith("ring_")}
"""dict: A moon's column key -> its source line (a moon holds no ring system of its own, GEN.732, and its inclination has no line)."""

_NEUTRON_STAR = {
    "magnetic_inclination_deg": "Isotropic at birth, aligned with the spin by the dipole's torque over about 7 million years [R5143]; a recycled millisecond pulsar is isotropic again.",
}
"""dict: A neutron star's column key -> its source line."""

SOURCES = {"planet": _PLANET, "moon": _MOON, "star": _STAR, "neutron_star": _NEUTRON_STAR}
"""dict: Each kind with sources -> `{column key: one-line source}`; a kind or key not here has no source line."""


def source_of(kind, key):
    """The one-line source of a row of a kind's sheet, or `None` when none is written."""
    return SOURCES.get(kind, {}).get(key)
