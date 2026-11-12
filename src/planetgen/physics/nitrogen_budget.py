# planetgen/physics/nitrogen_budget.py

"""
Where the nitrogen of a rocky body goes before its air: the core and the melt (GEN.1584, from the research of GEN.1583)
=====================================================================================================================

The air of a rocky body takes the nitrogen its mantle holds (`mixture_rocky.rocky_inventory`) times the share the volcanoes bring up
(`surface_pressure.degassed_share` x an efficiency). Both were wrong for nitrogen alone: the efficiency was 1 (all that degasses stays
as gas), carbon's share being the Earth's calibration, and nothing took nitrogen back to the metal. Two pieces fix that, both neutral for an Earth twin.

* *The metal takes nitrogen (A).* Nitrogen is siderophile: the metal/silicate partition coefficient `D_N` is 10 to 31 in the experiments (Suer and others 2023, their
  section 2.3.3), 20 +/- 10 measured and 1.8 thermodynamic in [R2171], nominal 10 with a range of 1 to 100 in Johansen and others 2023, and it
  rises with the oxygen fugacity: 1 (dIW -4) to 39 (dIW -1) at 20 GPa and 3000 K, 5 to 265 at 5 GPa and 2000 K (Suer and others 2023, Table 1). The mantle keeps
  `1 / (1 + D_N x)` of the nitrogen the body brought, `x` the core over the mantle mass. `D_N` here is the 20 GPa pair, log-linear in the redox (`metal_silicate_d`), held between 1
  and the 5 GPa oxidized value 265. The inventory's bracket is the Earth's, which already has the Earth's partition inside it, so only a body more oxidized than the reference
  redox (`REFERENCE_DIW`, the Earth-mass mean of the magma-ocean redox `a1`) loses more; a more reduced one keeps what its bracket gives (`mantle_retention`).
* *The melt keeps nitrogen (C).* Johnson and Goldblatt put 4 of 21 to 35 x 1e18 kg of the bulk silicate Earth's nitrogen in the air (11 to 19 percent; Table 4), where the
  carbon-calibrated degassed share of the Earth twin is 36 percent: the share of the degassed nitrogen that stays gas is 0.3 to 0.55 (`DEGAS_EFFICIENCY_N`, drawn). Their inventory
  (17 to 31e18 kg) is also 1.4 to 2.6 times the old middle of ours, so the nitrogen per carbon atom bracket (`mixture_rocky.INVENTORY`, 0.012 to 0.23) is re-centered on the delivered ratio (CI
  0.059, Venus 0.074; the bulk silicate Earth's 0.020 is the low end): the product of the two keeps the Earth twin at its 78 kPa.

Every number is a design constant with its tier and tweak level in `docs/design/research/gap-filling-plan.md` (section 6.281, NB1 to NB4).
"""

import math

REDUCED_DIW = -4.0
"""float: The redox of the reduced end of the partition table, log10 units above the iron-wuestite buffer (Suer and others 2023, Table 1)."""
OXIDIZED_DIW = -1.0
"""float: The oxidized end."""
D_N_REDUCED = 1.0
"""float: `D_N` at the reduced end, 20 GPa and 3000 K (Suer and others 2023, Table 1: 'embryo or intermediate magma ocean')."""
D_N_OXIDIZED = 39.0
"""float: `D_N` at the oxidized end, same conditions."""
D_N_MAX = 265.0
"""float: The most `D_N` takes (the oxidized value at 5 GPa and 2000 K of the same table)."""
CORE_FORMATION_OFFSET = -2.5
"""float: The redox the metal separated at, below the magma-ocean surface redox `a1` of the recipes: the Earth's surface magma ocean ended near IW+0.5 ([R4154] scale) while its core formed near IW-2 (the
`a1` row's own note: about IW-2 for a small body whose redox was set at core formation). Stand-in NB2 (tier T6, tweak level L2)."""
REFERENCE_DIW = 0.5
"""float: The magma-ocean redox (`a1`) at which the inventory's bracket already holds the Earth's partition: the Earth-mass mean."""

_SLOPE = math.log(D_N_OXIDIZED / D_N_REDUCED) / (OXIDIZED_DIW - REDUCED_DIW)


def metal_silicate_d(delta_iw):
    """`D_N`, the nitrogen concentration in the metal over that in the silicate, at the redox `delta_iw` (log-linear between the table's ends, held between 1 and `D_N_MAX`)."""
    return min(max(D_N_REDUCED * math.exp(_SLOPE * (delta_iw - REDUCED_DIW)), D_N_REDUCED), D_N_MAX)


def mantle_share(delta_iw, core_over_mantle):
    """The share of a body's nitrogen the silicate keeps against the metal: `1 / (1 + D_N x)`, with the metal's redox the magma ocean's `delta_iw` plus `CORE_FORMATION_OFFSET`."""
    return 1.0 / (1.0 + metal_silicate_d(delta_iw + CORE_FORMATION_OFFSET) * max(core_over_mantle, 0.0))


def mantle_retention(delta_iw, core_over_mantle):
    """
    The factor on the nitrogen of the inventory's bracket for a body of redox `delta_iw` and core over mantle ratio `core_over_mantle`: the mantle share against that at
    `REFERENCE_DIW`, never above 1 (a body more reduced than the reference keeps what the bracket gives; one more oxidized keeps less). 1 with no redox or no core.
    """
    if delta_iw is None or not core_over_mantle or core_over_mantle <= 0.0:
        return 1.0
    return min(1.0, mantle_share(delta_iw, core_over_mantle) / mantle_share(REFERENCE_DIW, core_over_mantle))
