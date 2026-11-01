# planetgen/physics/mixture_oxygen.py

"""
The most oxygen an air can hold (GEN.693's companion; Boss, 10-24: biotic worlds 35 percent, abiotic worlds 95 percent)
========================================================================================================================

No natural air is pure oxygen: nitrogen, argon (radiogenic argon alone is 0.93 percent of Earth's air after 4.5 Gyr) and carbon
dioxide are always there, and a living world is held near its fire limit (a Phanerozoic peak of about 30 to 35 percent; above 25 to
30 percent fires spread so easily that plants would struggle). An air that lost several oceans to photolysis and hydrogen escape can
reach 90 percent or more [R2769, R4771] but is capped. `cap_oxygen` applies a cap and a
floor of nitrogen plus argon to mole fractions; the abiotic cap is applied by `mixture_evolved` to the abiotic-oxygen family and the
biotic one is for the life stage (`BIOTIC_O2_MAX`), which sets the oxygen a living world has.

    cap_oxygen({"O2": 0.99, "CO2": 0.01}, ABIOTIC_O2_MAX, INERT_FLOOR)     # O2 0.945, N2 0.005, CO2 0.05
"""

BIOTIC_O2_MAX = 0.35
"""float: The most oxygen, by mole, of the air of a world with life (Boss 10-24, from the Phanerozoic peak and the fire limit)."""

ABIOTIC_O2_MAX = 0.95
"""float: The most oxygen of the air of a lifeless world that lost its water (Boss 10-24)."""

INERT_FLOOR = 0.005
"""float: The least nitrogen plus argon, by mole, of any air with oxygen above a few percent; a design floor (Earth has 78 and 0.9 percent)."""


def cap_oxygen(mole_fractions, maximum, inert_floor=INERT_FLOOR):
    """
    `mole_fractions` with the oxygen at most `maximum` and nitrogen plus argon at least `inert_floor`; the other gases keep their ratios.

    Args:
        mole_fractions (dict): Gas-table formula -> mole fraction, summing to 1.
        maximum (float): The oxygen cap.
        inert_floor (float): The least N2 + Ar; 0 for no floor.

    Returns:
        dict: The new mole fractions (a new dict, summing to 1, largest first); the input if it is already inside the limits.
    """
    oxygen = mole_fractions.get("O2", 0.0)
    inert = mole_fractions.get("N2", 0.0) + mole_fractions.get("Ar", 0.0)
    if oxygen <= maximum and (inert >= inert_floor or oxygen <= 0.0):
        return mole_fractions
    new_oxygen = min(oxygen, maximum)
    others = {formula: value for formula, value in mole_fractions.items() if formula != "O2"}
    others_sum = sum(others.values())
    if others_sum > 0.0:
        scale = (1.0 - new_oxygen) / others_sum
        result = {formula: value * scale for formula, value in others.items()}
    else:
        result = {"N2": 1.0 - new_oxygen}
    deficit = inert_floor - (result.get("N2", 0.0) + result.get("Ar", 0.0))
    if deficit > 0.0 and new_oxygen > 0.0:
        take = min(deficit, new_oxygen)
        result["N2"] = result.get("N2", 0.0) + take
        new_oxygen -= take
    if new_oxygen > 0.0:
        result["O2"] = new_oxygen
    total = sum(result.values())
    return {formula: value / total for formula, value in sorted(result.items(), key=lambda item: -item[1])}
