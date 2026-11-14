# planetgen/physics/surface_liquid_world.py

"""
The liquid a world holds on its surface (GEN.1035 wiring; the table and the rules are `surface_liquids`)
=========================================================================================================

Sets `COLUMN_FIELDS` on a planet or moon from its surface temperature, pressure, air, ices, ground and hydrosphere (called by
`surface_composition.generate_surface`, after the hydrosphere and the ground). It draws nothing. A body has at most one surface liquid:

1. **Water.** A world whose hydrosphere is a liquid surface (`surface ocean`, `brine seep`, `hycean`) has `water`, covering `ocean_fraction`
   to the depth the hydrosphere stored. Every other state of water (dry, vapor, ice, ice-covered ocean) leaves the surface to the rules below: a
   frozen or buried water ocean is not a surface liquid.
2. **Lava.** A surface hot enough to be molten (at least the liquidus of the ground's rock: ultramafic ground melts to komatiite, felsic ground to
   rhyolite, the rest to basalt, and a surface above the peridotite liquidus a peridotitic magma ocean) is covered by it all over: the temperature is the
   body's mean, so a mean above the liquidus is a molten surface (`physics.mixture_lava` is the vapor over it). The depth of a magma ocean is not
   known (no rule yet). A lava surface holds in its own vapor, so its boiling window is checked at least at one bar.
3. **Seas of the volatile ices and the air.** Otherwise `surface_liquids.dominant` picks the liquid that is stable at the surface temperature and pressure
   and that the world supplies: a volatile of the air at its saturation pressure condenses, and the ices the body holds (`ice_<species>_mass_fraction`)
   supply their liquid. A methane sea of a cold world with nitrogen air is Titan's. The share of the ice that lies on the surface as a liquid is
   `SURFACE_SHARE_OF_ICE`, spread over the surface like water (`hydrosphere.ocean_fraction`'s hypsometry); the liquid that condenses out of an
   air over saturation lies at the mass the excess partial pressure weighs (`(p - p_sat) / g` per square metre).

Every number that is not a published value is a stand-in in the gap registry (docs/design/research/gap-filling-plan.md section 6.69, rows LQ11 to LQ14).
"""

import math

from planetgen.physics import acid_medium, constants, hydrosphere, ice_species, mixture, planet_type, saturation, strange_solvents, surface_liquids

LIQUID_FIELDS = ("surface_liquid", "surface_liquid_fraction", "surface_liquid_depth_km")
"""tuple: What a body stores of the liquid itself: its id in `surface_liquids` (`None` for none), the share of the surface it covers, and its mean depth over
that share, km (`None` when not known)."""

SOLUTE_FIELDS = ("surface_liquid_acid_wt", "surface_liquid_acid_layer", "surface_liquid_ammonia_wt")
"""tuple: What a body stores of what its liquids hold (GEN.1215, docs/design/composition/surface-liquid-composition.md): the strength of liquid sulfuric acid,
percent of its mass, and where it is (`pool` at the ground or `cloud`), and the ammonia of a water ocean, percent of the liquid's mass. `None` where the
body has no such liquid."""

MIX_FIELDS = ("surface_liquid_methane_mol", "surface_liquid_ethane_mol", "surface_liquid_nitrogen_mol")
"""tuple: What a body stores of the makeup of a methane, ethane or nitrogen sea (GEN.1218b): the mole fraction of each, by Raoult's law with the nitrogen's reduced solubility. `None`
for any other liquid."""

SALT_FIELD = "surface_liquid_salt_wt"
"""str: The dissolved salt of a water sea, percent of the solution's mass, as sodium chloride that lowers the water activity as much as the sea's salts do (GEN.1218b)."""

SECOND_FIELDS = ("surface_liquid_2", "surface_liquid_2_fraction")
"""tuple: The second liquid of a body that holds two (GEN.1218b, Boss's answer 1 of 11:56 AM EDT 2026-10-31): its id and the share of the surface it covers beside the main one."""

COLUMN_FIELDS = LIQUID_FIELDS + SOLUTE_FIELDS + MIX_FIELDS + (SALT_FIELD,) + SECOND_FIELDS
"""tuple: Everything a body stores about its surface liquids (columns of `planets` and `moons`)."""

HYDROCARBON_SEAS = {"methane": "CH4", "ethane": "C2H6", "nitrogen": "N2"}
"""dict: The liquids whose makeup is stored as mole fractions (`MIX_FIELDS`), and the formula of each."""

OSMOTIC_COEFFICIENT = 0.93
"""float: STAND-IN GN4 (T5, L1): the osmotic coefficient of sodium chloride in a sea of 0.1 to 1 mol/kg (0.93 to 0.90 at 25 C), used to turn a water activity into the salt that makes it:
`ln a_w = -nu M_w phi m` with `nu` 2 and `M_w` 0.018015 kg/mol."""
NACL_G_PER_MOL = 58.44
"""float: The molar mass of sodium chloride, g/mol."""
NACL_SATURATION_ACTIVITY = 0.753
"""float: The water activity of a sodium chloride solution saturated at 25 C (753 per mille). A sea with a lower activity holds other salts (chlorides of calcium and magnesium, perchlorates)
and has no sodium chloride equivalent, so it stores no salinity (`_salinity`)."""
WATER_MOLAR_MASS_KG = 0.018015
"""float: The molar mass of water, kg/mol."""

WATER_ID = "water"
WATER_STATES = ("surface ocean", "brine seep", "hycean")
"""tuple: The hydrosphere states whose water is a liquid at the surface."""

SURFACE_SHARE_OF_ICE = 1.3e-3
"""float: STAND-IN LQ11 (L2): the share of a body's ice of a volatile that lies on the surface as its liquid. Tuned so that a Titan-like body with 1 percent of its
mass as methane ice (a guess, the ices are drawn) covers 1.7 percent of its surface, against the 1.6 percent the lakes and seas of Titan cover [R1857]. Retuned from 3e-4 when GEN.1147's two-level hypsometry replaced the cover law: it covers less of the surface for the same layer."""

ICE_LIQUID = {"CH4": "methane", "N2": "nitrogen", "CO2": "carbon_dioxide", "CO": "carbon_monoxide", "HCN": "hydrogen_cyanide",
              "SO2": "sulfur_dioxide"}
"""dict: The liquid each stored ice species supplies when it melts. Ammonia is not here: ice with water makes the ammonia-water solution (`AMMONIA_WATER_LIQUID`)."""
ETHANE_LIQUID = "ethane"
"""str: The liquid the ethane among the ices with no column of their own supplies when it melts (GEN.1218a)."""
AMMONIA_WATER_LIQUID = "ammonia_water"
"""str: The liquid the ammonia ice makes with the water ice of a cold body (GEN.1218a): the eutectic of the table, liquid from 176 K."""
AMMONIA_WATER_MASS_PERCENT = 32.0
"""float: The ammonia in `AMMONIA_WATER_LIQUID`, percent of its mass ([R2248] eutectic, the table's 32 percent). Source tier T2, tweak level L0 (registry LQ15)."""

LAVA_OF_GROUND = (("ultramafic", "komatiite_melt"), ("felsic", "rhyolite_melt"), ("basaltic", "basalt_melt"))
"""tuple: The lava a ground type of the surface melts to; every other ground (and a ground not drawn) melts to basalt. STAND-IN LQ12."""
LAVA_DEFAULT = "basalt_melt"
MAGMA_OCEAN = "peridotite_magma"
LAVA_IDS = (*(liquid for _ground, liquid in LAVA_OF_GROUND), LAVA_DEFAULT, MAGMA_OCEAN)
"""tuple: Every liquid id that is a lava; the depth of one is not known."""
LAVA_MIN_PRESSURE_PA = 1.0e5
"""float: A lava surface holds in its own vapor, so its window is checked at least at this pressure."""


def _liquid_water(body):
    body.surface_liquid = WATER_ID
    body.surface_liquid_fraction = body.ocean_fraction
    body.surface_liquid_depth_km = body.ocean_depth_km


def _lava(body, temperature_k, pressure_pa):
    """The molten surface of a body, or `None` when it is cooler than the liquidus of its ground's rock (a magma ocean over the peridotite liquidus)."""
    pressure = max(pressure_pa, LAVA_MIN_PRESSURE_PA)
    if surface_liquids.is_stable(MAGMA_OCEAN, temperature_k, pressure):
        return MAGMA_OCEAN
    lava = max(LAVA_OF_GROUND, key=lambda pair: getattr(body, f"ground_{pair[0]}_fraction", None) or 0.0)
    liquid = lava[1] if (getattr(body, f"ground_{lava[0]}_fraction", None) or 0.0) > 0.0 else LAVA_DEFAULT
    return liquid if surface_liquids.is_stable(liquid, temperature_k, pressure) else None


def _reservoir(body):
    """`{liquid id: mass fraction of the body}` of the ices that supply a liquid."""
    out = {}
    for formula, liquid in ICE_LIQUID.items():
        amount = getattr(body, ice_species.FIELD_OF[formula], None) if formula in ice_species.FIELD_OF else None
        if amount and amount > 0.0:
            out[liquid] = amount
    other = getattr(body, ice_species.FIELD_OF[ice_species.OTHER], None) or 0.0
    ethane = other * ice_species.ethane_share_of_other()
    if ethane > 0.0:
        out[ETHANE_LIQUID] = ethane
    ammonia = getattr(body, ice_species.FIELD_OF["NH3"], None) or 0.0
    water_ice = getattr(body, "water_mass_fraction", None) or 0.0
    if ammonia > 0.0 and water_ice > 0.0:         # a solution needs both: the ammonia fills the eutectic's 32 percent, the water the rest
        out[AMMONIA_WATER_LIQUID] = min(ammonia / (AMMONIA_WATER_MASS_PERCENT / 100.0), water_ice)
    return out


def _seas(body, temperature_k, pressure_pa):
    """`[(liquid id, surface mass per square metre in kg), ...]` of the seas a cold or temperate body could hold, best first (`surface_liquids.candidates`)."""
    partials = mixture.partial_pressures_pa(body)
    reservoir = _reservoir(body)
    gravity = (body.gravity or 0.0) * constants.EARTH_GRAVITY
    area_m2 = 4.0 * math.pi * (body.radius * 1000.0) ** 2
    out = []
    for pick in surface_liquids.candidates(temperature_k, pressure_pa, partials, reservoir):
        if pick.basis == "reservoir":
            out.append((pick.liquid_id, SURFACE_SHARE_OF_ICE * reservoir[pick.liquid_id] * body.mass / area_m2))
        else:
            row = surface_liquids.liquid(pick.liquid_id)
            excess = partials[row.saturation] - saturation.pressure_pa(row.saturation, temperature_k)
            out.append((pick.liquid_id, max(excess, 0.0) / gravity if gravity > 0.0 else 0.0))
    return out


def _cover(liquid_id, mass_per_m2, temperature_k, gravity):
    """`(share of the surface, mean depth over that share in km)` of a sea of `mass_per_m2` kg per square metre, or `None` if it covers none."""
    layer_km = mass_per_m2 / surface_liquids.density_kg_m3(liquid_id, temperature_k) / 1000.0
    covered = hydrosphere.ocean_fraction(layer_km, gravity)
    return (covered, layer_km / covered) if covered else None


def generate(body):
    """Sets `COLUMN_FIELDS` on a planet or moon (call after its hydrosphere and ground): `None` for a giant, a world with no liquid at its surface and a
    body saved before the columns."""
    for field in COLUMN_FIELDS:
        setattr(body, field, None)
    if not planet_type.has_ground(body):      # GEN.1334: nor does a planet under a deep envelope hold a liquid at a surface
        return
    _main_liquid(body)
    _second_liquid(body)
    _solutes(body)
    _mix(body)
    _salinity(body)


def acid_level(values):
    """`acid_medium`'s level (`NONE`, `CLOUD`, `POOL`) of a body's column-named `values`: read from the stored layer when the body has one, else worked out
    from its air."""
    layer = values.get("surface_liquid_acid_layer")
    if layer is None:
        return acid_medium.level(values)
    return acid_medium.POOL if layer == POOL_LAYER else acid_medium.CLOUD


def ammonia_water_mass_fraction(values):
    """The ammonia of a body's water as a share of the water's mass: read from the stored percent of the liquid when the body has one, else worked out
    from its bulk composition (`strange_solvents`)."""
    percent = values.get("surface_liquid_ammonia_wt")
    if percent is None:
        return strange_solvents.ammonia_water_mass_fraction(values)
    return percent / (100.0 - percent)


POOL_LAYER, CLOUD_LAYER = "pool", "cloud"
"""str: The places `surface_liquid_acid_layer` says liquid sulfuric acid is."""

AMMONIA_OCEAN_STATES = ("surface ocean", "brine seep", "ice-covered ocean", "hycean")
"""tuple: The hydrospheres whose water is a liquid ocean (under a lid or not), the ones that hold the body's ammonia as a solute."""


def _solutes(body):
    """Sets `SOLUTE_FIELDS`: the strength and place of the liquid sulfuric acid the air makes (`acid_medium`, GEN.888), and the ammonia in a liquid water
    ocean (`strange_solvents`, GEN.894), as a percent of the liquid's mass."""
    from planetgen.physics import habitability_world        # imported here: it reads this module's chain
    found = acid_medium.medium(habitability_world.values_from_body(body))
    if found.level != acid_medium.NONE:
        body.surface_liquid_acid_wt = found.weight_percent
        body.surface_liquid_acid_layer = POOL_LAYER if found.level == acid_medium.POOL else CLOUD_LAYER
    if getattr(body, "surface_liquid", None) == AMMONIA_WATER_LIQUID:
        body.surface_liquid_ammonia_wt = AMMONIA_WATER_MASS_PERCENT
    elif getattr(body, "hydrosphere", None) in AMMONIA_OCEAN_STATES:
        share = strange_solvents.ammonia_water_mass_fraction({"ice_mass_fraction": getattr(body, "ice_mass_fraction", None),
                                                              "water_mass_fraction": getattr(body, "water_mass_fraction", None)})
        if share:
            share *= _liquid_concentration(body)
            body.surface_liquid_ammonia_wt = min(100.0 * share / (1.0 + share), AMMONIA_WATER_MASS_PERCENT)


def _liquid_concentration(body):
    """
    How many times the ammonia of a water layer is concentrated in its liquid (GEN.1563): an ammonia-water ocean under an ice lid keeps the ammonia of the water that freezes out of it,
    because that ice holds almost none (the residual liquid of a freezing ocean climbs toward the eutectic, the way Titan's and Enceladus's oceans are thought to reach tens of percent).
    The water that freezes out of the ocean from above is the ice shell (`ice_shell_km`); the high-pressure ice beneath (`hp_ice_km`) is taken to keep its own ammonia, in hydrates, as the
    deep ice of a differentiated body does. By thickness (an ice's density is within a tenth of the liquid's): the factor is `(liquid + shell) / liquid`, 1 for a layer with no shell or a body
    with no ocean stored. Ammonia beyond the eutectic does not stay in the liquid: the caller holds the result at `AMMONIA_WATER_MASS_PERCENT`. Registry LQ16: tier 6, tweak level 1 (the
    high-pressure ice counted in raised the ammonia of nearly every ice-covered ocean to the eutectic, 75 of 76 in 180 seeds).
    """
    liquid = getattr(body, "ocean_depth_km", None) or 0.0
    if liquid <= 0.0:
        return 1.0
    return 1.0 + (getattr(body, "ice_shell_km", None) or 0.0) / liquid


def _main_liquid(body):
    """Sets `LIQUID_FIELDS`: water, lava, or the sea of a volatile."""
    if getattr(body, "hydrosphere", None) in WATER_STATES:
        _liquid_water(body)
        return
    temperature = body.surface_temperature or 0.0
    pressure = body.atmospheric_pressure or 0.0
    if not temperature > 0.0:
        return
    lava = _lava(body, temperature, pressure)
    if lava is not None:
        body.surface_liquid, body.surface_liquid_fraction = lava, 1.0
        return
    if not pressure > 0.0:
        return                              # in vacuum a liquid boils away
    seas = _seas(body, temperature, pressure)
    if not seas:
        return
    liquid_id, mass_per_m2 = seas[0]
    liquid_id = _named_by_its_bulk(body, liquid_id, temperature, pressure)
    found = _cover(liquid_id, mass_per_m2, temperature, (body.gravity or 0.0) * constants.EARTH_GRAVITY)
    if found is not None:
        body.surface_liquid, body.surface_liquid_fraction = liquid_id, found[0]
        body.surface_liquid_depth_km = found[1]


LABEL_FLOOR = 0.5
"""float: The least mole fraction a hydrocarbon sea's label holds in its stored mix (GEN.1249). The supplied liquid is the bulk of the sea, so this holds unless the gases dissolved in it
outweigh it (a dense, cold nitrogen air), and then the sea is named for the gas that makes most of it. A rule, not a fitted number (the label is the majority of the sea)."""


def _named_by_its_bulk(body, liquid_id, temperature_k, pressure_pa):
    """
    The name of a methane, ethane or nitrogen sea (GEN.1249): the liquid the supply gives, which is the bulk of the sea and the solvent the air's gases dissolve in (`_mix`). When the
    dissolved gases outweigh it (it holds under `LABEL_FLOOR` of the sea) the sea is named for the stable liquid that makes most of it, so the label always holds the majority of the stored mix.
    """
    for _ in HYDROCARBON_SEAS:
        if liquid_id not in HYDROCARBON_SEAS:
            return liquid_id
        fractions = surface_liquids.ideal_sea(temperature_k, mixture.partial_pressures_pa(body), solvent=HYDROCARBON_SEAS[liquid_id])
        own = fractions.get(HYDROCARBON_SEAS[liquid_id], 0.0)
        if own >= LABEL_FLOOR or own <= 0.0:
            return liquid_id            # the bulk, or a sea of its own liquid alone (`_mix`)
        by_formula = {formula: name for name, formula in HYDROCARBON_SEAS.items()}
        larger = [name for formula, name in by_formula.items()
                  if fractions.get(formula, 0.0) > own and surface_liquids.is_stable(name, temperature_k, pressure_pa)]
        if not larger:
            return liquid_id
        liquid_id = max(larger, key=lambda name: fractions[HYDROCARBON_SEAS[name]])
    return liquid_id


def _second_liquid(body):
    """Sets `SECOND_FIELDS`: the best other sea the body supplies and has room for beside its main liquid (a lava surface is molten all over and has none)."""
    main = getattr(body, "surface_liquid", None)
    if main is None or main in LAVA_IDS:
        return
    room = 1.0 - (getattr(body, "surface_liquid_fraction", None) or 0.0)
    temperature, pressure = body.surface_temperature or 0.0, body.atmospheric_pressure or 0.0
    if room <= 0.0 or not temperature > 0.0 or not pressure > 0.0:
        return
    gravity = (body.gravity or 0.0) * constants.EARTH_GRAVITY
    for liquid_id, mass_per_m2 in _seas(body, temperature, pressure):
        if liquid_id == main:
            continue
        found = _cover(liquid_id, mass_per_m2, temperature, gravity)
        if found is not None:
            body.surface_liquid_2, body.surface_liquid_2_fraction = liquid_id, min(found[0], room)
            return


def _mix(body):
    """
    Sets `MIX_FIELDS` for a methane, ethane or nitrogen sea (GEN.1249): the liquid the supply gives is the bulk of the sea, and the other gases of the air dissolve in it by Raoult's
    law (`surface_liquids.ideal_sea` with that liquid as the solvent: each volatile `x = p / P_sat`, nitrogen's reduced solubility included, the supplied liquid the rest). The label
    and the mix agree: the labelled liquid is at least what the dissolved gases leave of the sea, whatever the air holds of its own gas. Ethane is in a sea only when it is the supplied
    liquid (nothing in the generator makes it from methane's photolysis). A sea whose dissolved gases overfill it, or that has no gas of the air to dissolve, is taken as that liquid alone.
    """
    main = getattr(body, "surface_liquid", None)
    if main not in HYDROCARBON_SEAS:
        return
    solvent = HYDROCARBON_SEAS[main]
    fractions = surface_liquids.ideal_sea(body.surface_temperature, mixture.partial_pressures_pa(body), solvent=solvent)
    if fractions.get(solvent, 0.0) <= 0.0:
        fractions = {solvent: 1.0}
    for field, formula in zip(MIX_FIELDS, HYDROCARBON_SEAS.values()):
        setattr(body, field, fractions.get(formula, 0.0))


def _salinity(body):
    """Sets `SALT_FIELD` for a water sea from its water activity: the sodium chloride that lowers it as far (`ln a_w = -2 M_w phi m`), percent of the solution's mass; none under the activity of a saturated sodium chloride solution."""
    activity = getattr(body, "water_activity", None)
    if getattr(body, "surface_liquid", None) != WATER_ID or not activity or not NACL_SATURATION_ACTIVITY <= activity <= 1.0:
        return
    molality = -math.log(activity) / (2.0 * WATER_MOLAR_MASS_KG * OSMOTIC_COEFFICIENT)
    grams = molality * NACL_G_PER_MOL
    body.surface_liquid_salt_wt = 100.0 * grams / (1000.0 + grams)


