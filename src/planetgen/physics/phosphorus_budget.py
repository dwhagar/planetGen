# physics/phosphorus_budget.py

"""
The phosphorus budget of a rocky body (GEN.1367): the star's phosphorus through the core, the crust and its minerals to the phosphate in the ocean that the life scores read.

Phosphorus is one chain, and each link is a number a body already stores or this module derives:

1. **The mantle's P2O5** (`rock_p2o5_fraction`, GEN.1362, `mantle_oxides`): the star's P, the share its formation zone kept (condensation T50 1229 K, a leading-term volatility cut) and the
   share the core took (94 percent for the Earth). A reduced body has already lost most of it to the metal, and a Mars-like one (more FeO, less metal) keeps ten times the Earth's.
2. **The crust's P2O5.** Phosphorus is incompatible: the melt that makes the volcanic crust carries `rock_types.batch_enrichment(F)` times the mantle's, with the melt fraction F of the body's
   heat flow (`surface_composition.melt_fraction`). The Earth's mantle (0.021 weight percent) gives its crust 0.14.
3. **The host of the crust's P.** Under a low oxygen fugacity phosphorus is a phosphide, and the crust's P sits in schreibersite, (Fe,Ni)3P: the share is a logistic of the mantle's redox
   (`mantle_delta_iw`) about `PHOSPHIDE_MIDPOINT_DIW`. The rest is phosphate: apatite where the body has liquid water (it needs the volatiles), merrillite where it is dry.
4. **The weathering supply.** Apatite dissolves, merrillite faster, schreibersite barely without oxygen and water together (meteoritic schreibersite corrodes to phosphite and phosphate, which
   is slow). The supply is the crust's P2O5 times the reactive share of its hosts, times the erosion of the lid (a mobile lid uplifts and strips, a stagnant one does not), times the area the
   rain reaches (the land of a surface ocean; only the sea floor of one under ice), over the volume of the ocean the phosphate dissolves in.
5. **The ocean's phosphate.** The Earth is the one calibration: its 2.3 umol/L comes from its 0.14 percent crust, its 29 percent of land and its 1.3e18 m3 of ocean, so a body's phosphate is
   2.3 times each of its ratios to the Earth's. A sealed ocean (liquid between ice shell and high-pressure ice) has no contact with rock and holds the floor of meteoritic and hydrothermal
   supply, `FLOOR_UMOL_L`; a soda lake has no free calcium to precipitate apatite [R4399], so its phosphate climbs until the ceiling of the ocean's solubility.

Left out, and said so: the burial of phosphorus with iron oxides in an oxygenated ocean (the Earth's calibration has it, another oxidation state would change it), and the recycling by
a biosphere (the life timeline of GEN.750 runs after this and does not feed back). The numbers that stand in are in the registry (GEN.725, `docs/design/research/gap-filling-plan.md` section 6.221).

`phosphorus` (high, limited, starved) is kept as the name of the band of the number (`phosphorus_class`), for a page that wants words.
"""

import math

from planetgen.physics import constants, interior_heat, rock_types, surface_composition

EARTH_CRUST_P2O5_WT = 0.14
"""float: The Earth's crust, weight percent P2O5: continental 0.13, MORB 0.15 (registry row PB1)."""
EARTH_OCEAN_UMOL_L = 2.3
"""float: The Earth's ocean phosphate, umol/L: the deep ocean's 2 to 3 (registry row PB2; also `habitability_world`'s old `limited` value)."""
EARTH_LAND_FRACTION = 0.29
"""float: The Earth's land share of its surface."""
EARTH_OCEAN_VOLUME_M3 = 1.335e18
"""float: The Earth's ocean, cubic metres."""
EARTH_AREA_M2 = 4.0 * math.pi * (constants.EARTH_RADIUS_KM * 1e3) ** 2
"""float: The Earth's surface, square metres."""

PHOSPHIDE_MIDPOINT_DIW = -3.0
PHOSPHIDE_SCALE_DIW = 0.6
"""The mantle redox (log10 above iron-wuestite) at which half of the crust's P is a phosphide and the e-fold width of the logistic (registry rows PB3, PB4)."""
APATITE_WET_SHARE = 1.0
APATITE_DRY_SHARE = 0.2
"""float: The share of the phosphate that is apatite in a crust with liquid water and in a dry one; merrillite is the rest (registry row PB5)."""
MERRILLITE_FACTOR = 3.0
SCHREIBERSITE_FACTOR = 0.05
"""float: How fast merrillite and schreibersite weather against apatite, 1 (registry rows PB6, PB7)."""
STAGNANT_EROSION = 0.15
"""float: The erosion of a stagnant lid against a mobile one (registry row PB8)."""
SEAFLOOR_AREA_FRACTION = 0.05
"""float: The share of the surface whose rock the water reaches with no land above it (the sea floor's weathering; registry row PB9)."""
SODA_RESIDENCE = 400.0
"""float: How much longer phosphate stays dissolved in a soda lake, which has no calcium to take it out (registry row PB10)."""
CEILING_UMOL_L = 1000.0
"""float: The most phosphate an ocean holds, umol/L: the soda lakes' (registry row PB11)."""
FLOOR_UMOL_L = 0.01
"""float: The least, umol/L: meteoritic and hydrothermal supply into a sealed ocean (registry row PB12)."""
STARVED_BELOW_UMOL_L = 0.1
HIGH_FROM_UMOL_L = 100.0
"""float: The bands of `phosphorus_class` (registry row PB13)."""

SURFACE_RAIN_STATES = ("surface ocean", "brine seep", "hycean")
"""tuple: The hydrospheres that rain on land; an ocean under ice is reached by its sea floor only."""
LIQUID_STATES = SURFACE_RAIN_STATES + ("ice-covered ocean",)
"""tuple: The hydrospheres with liquid water in the crust's reach, the apatite's."""

COLUMN_FIELDS = ("crust_p2o5_wt", "crust_apatite_share", "crust_merrillite_share", "crust_schreibersite_share", "ocean_phosphate_umol_l", "phosphorus")
"""tuple: What `generate` sets on a planet or moon: the crust's P2O5 (weight percent), the shares of its P in the three hosts (they add to 1), the ocean's phosphate (umol/L) and its band."""


def crust_p2o5_wt(rock_p2o5_fraction, melt_fraction):
    """The crust's P2O5 in weight percent: the mantle's (a mass fraction of the rock) times the batch melt's enrichment."""
    return 100.0 * rock_p2o5_fraction * rock_types.batch_enrichment(melt_fraction)


def host_shares(delta_iw, wet):
    """The shares of the crust's phosphorus in apatite, merrillite and schreibersite, from the mantle's redox (`None` is the Earth's, all phosphate) and whether it has liquid water."""
    if delta_iw is None:
        phide = 0.0
    else:
        phide = 1.0 / (1.0 + math.exp((delta_iw - PHOSPHIDE_MIDPOINT_DIW) / PHOSPHIDE_SCALE_DIW))
    apatite = APATITE_WET_SHARE if wet else APATITE_DRY_SHARE
    return (1.0 - phide) * apatite, (1.0 - phide) * (1.0 - apatite), phide


def reactive_share(shares):
    """What of a crust's P2O5 the water can take, against an all-apatite crust's 1: the hosts' shares times how fast each weathers."""
    apatite, merrillite, schreibersite = shares
    return apatite + MERRILLITE_FACTOR * merrillite + SCHREIBERSITE_FACTOR * schreibersite


def phosphorus_class(umol_l):
    """The band of an ocean's phosphate: `starved`, `limited` or `high`."""
    return "starved" if umol_l < STARVED_BELOW_UMOL_L else "high" if umol_l >= HIGH_FROM_UMOL_L else "limited"


def ocean_phosphate(crust_wt, shares, erosion, area_ratio, volume_ratio, kind):
    """
    The ocean's phosphate in umol/L: the Earth's times the ratios of the crust's P2O5, its reactive share, its erosion and the rock area the rain reaches, over the ocean's volume
    (`volume_ratio` against the Earth's), longer-lived in a soda lake and the floor for a sealed ocean; never above `CEILING_UMOL_L` nor below `FLOOR_UMOL_L`.
    """
    if kind == "ice-sealed":
        return FLOOR_UMOL_L
    supply = (crust_wt / EARTH_CRUST_P2O5_WT) * reactive_share(shares) * erosion * area_ratio
    value = EARTH_OCEAN_UMOL_L * supply / max(volume_ratio, 1e-12) * (SODA_RESIDENCE if kind == "soda" else 1.0)
    return min(CEILING_UMOL_L, max(FLOOR_UMOL_L, value))


def _flux_of(body):
    found = interior_heat.update_heat(body)
    return found.own_flux_w_m2 if found is not None else 0.0


def generate(body):
    """
    Sets `COLUMN_FIELDS` on a planet or moon from its stored rock, melt, redox, lid and ocean (call after the ground and the water; no draw). `None` for a giant, a carbon rock or a row with no
    stored P2O5 (from before GEN.1362); the ocean's phosphate and band are `None` for a body with no liquid ocean.
    """
    for name in COLUMN_FIELDS:
        setattr(body, name, None)
    p2o5 = getattr(body, "rock_p2o5_fraction", None)
    if getattr(body, "body_type", None) in ("g", "a") or p2o5 is None:
        return
    wet = getattr(body, "hydrosphere", None) in LIQUID_STATES
    melt = surface_composition.melt_fraction(_flux_of(body))
    crust = crust_p2o5_wt(p2o5, melt)
    shares = host_shares(getattr(body, "mantle_delta_iw", None), wet)
    body.crust_p2o5_wt = crust
    body.crust_apatite_share, body.crust_merrillite_share, body.crust_schreibersite_share = shares
    kind = getattr(body, "ocean_class", None)
    if kind is None:
        return
    erosion = 1.0 if interior_heat.lid_of(body) == interior_heat.MOBILE else STAGNANT_EROSION
    rained = body.hydrosphere in SURFACE_RAIN_STATES
    area_fraction = max(getattr(body, "land_fraction", None) or 0.0 if rained else 0.0, SEAFLOOR_AREA_FRACTION)
    area_m2 = 4.0 * math.pi * (body.radius * 1e3) ** 2
    area_ratio = area_fraction / EARTH_LAND_FRACTION * area_m2 / EARTH_AREA_M2
    volume = (body.ocean_depth_km or 0.0) * 1e3 * (getattr(body, "ocean_fraction", None) or 0.0) * area_m2
    body.ocean_phosphate_umol_l = ocean_phosphate(crust, shares, erosion, area_ratio, volume / EARTH_OCEAN_VOLUME_M3, kind)
    body.phosphorus = phosphorus_class(body.ocean_phosphate_umol_l)
