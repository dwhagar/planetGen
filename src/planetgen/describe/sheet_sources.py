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

_BODY = {
    "eccentricity": "A planet draws from the Kepler multi-planet law [R1307] and is circularised by its tide; a regular "
                    "moon keeps a small forced value, a captured moon a thermal one.",
    "flux_variation": "The swing of the starlight between periapsis and apoapsis, peak to peak over the mean, from the eccentricity.",
    "tidal_flux_w_m2": "[R3464] tidal power from the stored eccentricity and the body's rheology.",
    "internal_heat_flux_w_m2": "Radiogenic plus primordial heat in one-box parameterised mantle convection ([R4339]; "
                               "[R1716]).",
    "urey_ratio": "The radiogenic share of the interior heat, an output of the mantle model and not a draw.",
    "lid_regime": "One of " + _LID + ".",
    "c_th0_ppm": "Follows the host star's thorium (its Eu relation, truncated at plus or minus 0.5 dex), not a draw of its own.",
    "c_u0_ppm": "Worked out from the thorium by the r-process production ratio and the decay of the galaxy's earlier ejecta.",
    "rotation_period_hours": "Birth spin worn down by the tide of what the body orbits, to a free, locked or resonant state "
                             "[R1584].",
    "axial_tilt_deg": "Isotropic for a rocky planet, Rayleigh about 15 degrees for a gas giant or an untidied moon, and on the orbit "
                      "normal for a locked body.",
    "magnetic_moment_a_m2": "The dipole of the dynamo that runs in the body: a rocky body's liquid core (its thermal history) or a giant's metallic or "
                            "ionic layer (the heat it radiates), by the buoyancy-flux scaling and the spin.",
    "ring_outer_km": "Ends at the Roche limit, 2.456 R (planet density over particle density)^(1/3), and inside the innermost large moon.",
    "ring_age_yr": "The years since the moon that fell inside the Roche limit, or the comet that passed inside it, made the ring; the ring thins by 1/e every half billion years.",
    "ring_optical_depth": "An output of the drawn mass, the extent and the particle size, 3 M / (4 density r area).",
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
    "nuclear_arrival_age_gy": "The age MIST's tracks put the star on the main sequence [R0773]; the star counts as pre-main-sequence "
                              "only until its radius has settled, the earlier age BHAC15 gives [R0254].",
    "disk_kind": "A young star holds a class 0, I or II disk by its age; a main-sequence star has a cold debris belt at the share "
                 "the DEBRIS survey found for its kind, 2 to 24 percent [R4081].",
    "disk_inner_au": "A young star's disk starts at the dust sublimation radius; a debris belt starts at 0.7 of its radius, the "
                     "blackbody radius of its drawn dust temperature times 1.74 [R4081].",
    "disk_outer_au": "A young star's disk ends log-uniform between 10 and 300 AU; a debris belt ends at 1.5 of its radius, the "
                     "blackbody radius of its drawn dust temperature times 1.74 [R4081].",
    "chemistry_branch": "A star older than the thick disc's start is alpha-rich with [R0349]'s logistic probability of its "
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

SOURCES = {"planet": _PLANET, "moon": _MOON, "star": _STAR}
"""dict: Each kind with sources -> `{column key: one-line source}`; a kind or key not here has no source line."""


def source_of(kind, key):
    """The one-line source of a row of a kind's sheet, or `None` when none is written."""
    return SOURCES.get(kind, {}).get(key)
