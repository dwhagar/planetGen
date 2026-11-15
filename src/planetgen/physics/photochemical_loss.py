# planetgen/physics/photochemical_loss.py

"""
The ammonia and methane the star's ultraviolet takes out of an air (GEN.1732)
=============================================================================

The air of a world holds the ammonia and methane its volcanoes and its comets gave it, and nothing in the model destroyed them: an air that held 1e-4 to 1e-2 of ammonia kept it for ever, where the star's
far ultraviolet breaks ammonia in hours to years (`photolysis`, GEN.496: 6.7e3 s unshielded at 1 AU from the quiet Sun) and methane in years to decades, and the Earth's air holds ammonia
at parts per billion. On the Abundant meter that kept worlds at 0.34 to 0.36 of Earth's sunlight warm (284 to 286 K) with 68 Pa to 6.5 kPa of ammonia in their air and no carbon dioxide
(`reports/gen1732-photochem-prototype-measure-2026-11-15.md`). This module is the law of the loss. It is not wired into the one-shot air solve; the history of the air (`air_history`, GEN.1745)
calls it at every step.

**The law.** A gas the star's light breaks stands against a source, the gas its volcanoes and comets put into the air, and a loss at the photolysis rate `k` (`photolysis.rate_per_s`
of the star's far ultraviolet at the body's distance) per molecule, so

    dN/dt = S - k N,       N(t + dt) = N e^(-k dt) + (S / k) (1 - e^(-k dt))       (`advance`)

and a source that has run for a time `age` leaves `(S / k) (1 - e^(-k age))`, which is `S / k` (`standing`: the amount at which the loss equals the source) once the age is over a few lives of the gas, and never more than the inventory the source has put into the air. The source of a
gas that is given as one inventory over the body's age is the inventory over the age (`standing_from_inventory`), which is what the history's steps replace with the body's own supply at each age.

**The light's limit.** The far ultraviolet is absorbed in the column of one optical depth (`photolysis.shielded_column_mol_m2`, Hu's critical-pressure rule [R2026]); a source above what the light can break in that
column (`photolysis.photon_limited_flux_mol_m2_s`) is not destroyed as it comes: the air gains the difference (`standing_from_inventory` and `advance` both take it).

**What the loss leaves.** The hydrogen of the gas leaves the air (`LOST_TO`): a mole of ammonia leaves half a mole of nitrogen, a mole of methane leaves a mole of carbon dioxide (the carbon stays, the oxygen is the water's).

Stand-ins: `docs/design/research/gap-filling-plan.md` section 6.299.
"""

import math

from planetgen.physics import gases, photolysis

PHOTOLYSED = ("NH3", "CH4")
"""tuple: The gases of an air the loss acts on (the two the model held for ever; GEN.1732)."""

LOST_TO = {"NH3": ("N2", 0.5), "CH4": ("CO2", 1.0)}
"""dict: Formula of a photolysed gas to (the gas its loss leaves, the moles of it left per mole lost). The hydrogen escapes and is not stored; the oxygen of the CO2 is the water's (T7, L1: GEN.725 PC4)."""

YEAR_S = 365.25 * 86400.0


def rate_per_s(formula, far_uv, star_temperature_k):
    """The photolysis rate, per second, of `formula` under `far_uv` (`aerosols.FarUv`: the far ultraviolet flux relative to the quiet Sun at 1 AU, and its chromospheric share) from a star of `star_temperature_k`; 0 with no far ultraviolet."""
    if far_uv is None or not far_uv.ratio or far_uv.ratio < 0.0 or not star_temperature_k:
        return 0.0
    return photolysis.rate_per_s(formula, far_uv.ratio, star_temperature_k, chromospheric_share=far_uv.chromospheric_share or 0.0)


def standing(source_mol_m2_s, rate):
    """The amount, mol per m2 of ground, at which a loss of `rate` per second equals a source of `source_mol_m2_s`; infinity when there is no loss."""
    if rate <= 0.0:
        return math.inf
    return max(source_mol_m2_s, 0.0) / rate


def advance(amount_mol_m2, source_mol_m2_s, rate, seconds, light_limit_mol_m2_s=math.inf):
    """
    The amount, mol per m2, of a gas after `seconds` with the source `source_mol_m2_s` and the loss `rate` per molecule per second (the exact step of `dN/dt = S - k N`, stable for any step). The light breaks no more than
    `light_limit_mol_m2_s` (`photolysis.photon_limited_flux_mol_m2_s`) in the column that absorbs it: where the loss would pass it the amount is advanced at that loss.
    """
    amount = max(amount_mol_m2, 0.0)
    if seconds <= 0.0:
        return amount
    source = max(source_mol_m2_s, 0.0)
    if rate <= 0.0:
        return amount + source * seconds
    decay = math.exp(-rate * seconds)
    exact = amount * decay + source / rate * (1.0 - decay)
    if light_limit_mol_m2_s >= math.inf or rate * max(amount, exact) <= light_limit_mol_m2_s:
        return exact
    return max(amount + (source - light_limit_mol_m2_s) * seconds, 0.0)


def standing_from_inventory(inventory_mol_m2, age_s, rate, light_limit_mol_m2_s=math.inf):
    """
    The amount, mol per m2, of a gas of which the body's volcanoes and comets put `inventory_mol_m2` into the air over `age_s` (a source of the inventory over the age, constant), after the star's light has worked on it all that
    time: `S tau (1 - exp(-age / tau))` for the source `S = inventory / age` and the life `tau = 1 / rate` of the gas, which is the inventory less what was lost when the life is long and the source's `S tau` when it is short. Where the
    source is over what the light can break (`light_limit_mol_m2_s`) the air first fills to the absorbing column (`light_limit / rate`) and then gains the difference of source and limit.
    """
    inventory = max(inventory_mol_m2, 0.0)
    if inventory <= 0.0 or rate <= 0.0 or age_s <= 0.0:
        return inventory
    source = inventory / age_s
    life = 1.0 / rate
    if source > light_limit_mol_m2_s:
        column = light_limit_mol_m2_s / rate
        reached = -life * math.log1p(-light_limit_mol_m2_s / source)        # when the amount stands at the absorbing column
        if reached < age_s:
            return column + (source - light_limit_mol_m2_s) * (age_s - reached)
    return source * life * -math.expm1(-age_s / life)


def column_mol_m2(partial_pa, formula, gravity_m_s2):
    """The column, mol per m2 of ground, of a gas of `partial_pa` in an air of surface gravity `gravity_m_s2` (the weight of its part of the air over the molar mass)."""
    return partial_pa / (gravity_m_s2 * gases.gas(formula).molar_mass)


def partial_pa_of(column, formula, gravity_m_s2):
    """The partial pressure, Pa, of a gas of `column` mol per m2."""
    return column * gravity_m_s2 * gases.gas(formula).molar_mass


def light_limit(formula, far_uv, star_temperature_k):
    """The most the light can break of `formula` under `far_uv`, mol per m2 of ground per second: the column that absorbs it times the rate (`photolysis.photon_limited_flux_mol_m2_s`); infinity with no far ultraviolet."""
    if far_uv is None or not far_uv.ratio or not star_temperature_k:
        return math.inf
    return photolysis.photon_limited_flux_mol_m2_s(formula, far_uv.ratio, star_temperature_k, chromospheric_share=far_uv.chromospheric_share or 0.0)


def lose(partials_pa, gravity_m_s2, age_s, far_uv, star_temperature_k):
    """
    The partial pressures, Pa, of an air (`partials_pa`: formula to Pa) after the star's light has worked on its ammonia and methane for `age_s`, the amounts of the two taken as the inventory of a constant source over the age
    (`standing_from_inventory`): what is lost leaves its products (`LOST_TO`). Moles are conserved but for the hydrogen; the other gases are returned as they were.
    """
    out = dict(partials_pa)
    for formula in PHOTOLYSED:
        have = out.get(formula, 0.0)
        if have <= 0.0:
            continue
        rate = rate_per_s(formula, far_uv, star_temperature_k)
        if rate <= 0.0:
            continue
        column = column_mol_m2(have, formula, gravity_m_s2)
        left = standing_from_inventory(column, age_s, rate, light_limit(formula, far_uv, star_temperature_k))
        lost = column - left
        out[formula] = partial_pa_of(left, formula, gravity_m_s2)
        product, per_mole = LOST_TO[formula]
        out[product] = out.get(product, 0.0) + partial_pa_of(lost * per_mole, product, gravity_m_s2)
    return out
