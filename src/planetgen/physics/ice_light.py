# planetgen/physics/ice_light.py

"""
Light through ice and snow (GEN.1726, docs/design/habitability/light-through-ice.md)
===================================================================================

A frozen ocean or lake is not dark. Photosynthetic light passes through the ice and snow on top of it as it passes through water: each layer takes a fixed share
of what reaches it for every metre, so what is left under a lid is the light that got in at the top (the surface reflects the rest) times one exponential for the ice
and one for the snow:

    T = (1 - albedo) * exp(-k_ice * h_ice) * exp(-k_snow * h_snow)

`T` is the share of the photosynthetically active light (400 to 700 nm) at the surface that reaches the water under the lid. It is the two-layer Beer-Lambert law of
[R5264] (equation 3.1, fitted to light under Antarctic pack ice) with the reflection at the top that [R5265] take off before they fit the same law to a perennially
ice-covered lake. Snow is the main absorber: about 30 times stronger per metre than sea ice and about 100 times stronger than clear lake ice, so a world's snow depth on
its ice decides most of what gets through.

The module is pure numbers: it reads no body and draws nothing. It does not say how thick the ice or the snow is (that is the climate's, and no law for it is built yet),
and the surface light is the caller's, taken from the world's own star and orbit, so a dim star moves every depth this law gives. The one use built is the light-limited
production of ocean under ice (GEN.1723): `light_below(surface, ...)` and the inverse `thickest_ice_m`.

Where the coefficients come from (the rows IL1 to IL7 of section 6.293 of the gap plan; `docs/design/habitability/light-through-ice.md` has the checks):

* `K_ICE_PER_M`: sea ice 1.26 per metre, measured 1.26 +- 0.13 and expected 0.8 to 1.5 ([R5264, R5267] median 0.8 at 600 nm); clear fresh lake ice 0.32 (Lake A, [R5265]) and as low as 0.2 for thick
  perennial Antarctic lake ice. The ice of a salty ocean takes the sea value (brine pockets and air bubbles scatter), the ice of a fresh lake the lake value.
* `K_SNOW_PER_M`: 30 per metre, measured 31.8 +- 0.7 on Antarctic pack ice and expected 10 to 100 [R5264]; a single wind-packed layer on a lake measured 2.1 [R5265], and temperate snow
  3 to 30 [R5268], which sets the low end.
* The albedos: bare sea ice 0.5 to 0.7, clear lake ice with its scattering top layer ("candled", 0.84), snow 0.8 to 0.97, melting ice 0.2 to 0.6 ([R5266]: 0.78, 0.60, 0.44 and 0.20 as the ice warms and melts).

What it leaves out. The coefficients are for the Sun's spectrum: ice absorbs least near 400 to 500 nm and more toward the red, so under a cool red star the light under a lid falls faster
than this law says, and under a blue-white star slower. The surface layer of sea ice (the top 10 to 20 cm) attenuates more than the interior; it is folded into the albedo and the bulk k.
Ice algae on the underside of the ice shade themselves (k 20 to 50, [R5266]) and melt ponds let more light through; neither is counted. Every constant is a stand-in in section 6.293 of
`docs/design/research/gap-filling-plan.md` with its tier and tweak level.

Use:

    t = transmission(ice_m=2.0, snow_m=0.0, kind=FRESH)          # a 2 m lid of lake ice, no snow: 0.08
    par_below = light_below(400.0, ice_m=1.0, snow_m=0.1)         # umol photons m-2 s-1 under 1 m of sea ice with 10 cm of snow
    h = thickest_ice_m(1.0e-4, kind=SEA)                          # how thick a lid still passes one part in ten thousand: about 6.6 m
"""

import math

SEA = "sea"
"""str: Ice of a salty ocean (brine pockets, air bubbles)."""
FRESH = "fresh"
"""str: Ice of a fresh lake or a fresh-water ocean (clear, few bubbles)."""
KINDS = (SEA, FRESH)

K_ICE_PER_M = {SEA: 1.26, FRESH: 0.32}
"""dict: By kind, the attenuation of photosynthetic light by the ice, per metre. Sea ice: measured 1.26 +- 0.13 under Antarctic pack ice [R5264]; fresh lake ice 0.32, Lake A [R5265]. [IL1, IL2, T1, L0.]"""
K_ICE_RANGE_PER_M = {SEA: (0.8, 1.5), FRESH: (0.2, 0.32)}
"""dict: By kind, the range of the ice's attenuation, per metre: sea ice 0.8 to 1.5 [R5264, R5267], thick clear lake ice 0.2 to 0.32 [R5265]. [IL1, IL2, T1, L0.]"""

K_SNOW_PER_M = 30.0
"""float: The attenuation of photosynthetic light by snow, per metre: measured 31.8 +- 0.7 on Antarctic pack ice [R5264]. [IL3, T1, L0.]"""
K_SNOW_RANGE_PER_M = (10.0, 100.0)
"""tuple: The range of the snow's attenuation, per metre: the expected 10 to 100 of [R5264] (grain size, density and impurity), whose low end [R5265, R5268] reach (2.1 for wind-packed snow, 3 to 30 temperate). [IL3, T1, L0.]"""

ALBEDO_BARE = {SEA: 0.6, FRESH: 0.84}
"""dict: By kind, the albedo of the bare ice surface to photosynthetic light: sea ice 0.5 to 0.7, center 0.6 [R5264, R5266]; clear lake ice with its scattering top layer 0.84 [R5265]. [IL4, T1, L0.]"""
ALBEDO_BARE_RANGE = {SEA: (0.5, 0.7), FRESH: (0.7, 0.9)}
"""dict: By kind, the range of the bare ice albedo. The fresh range is the stand-in's own (the one measurement, 0.84, with the width of the sea ice's). [IL4, T6, L1.]"""
ALBEDO_SNOW = 0.9
"""float: The albedo of a snow surface to photosynthetic light: 0.8 to 0.97, the middle ([R5264]; 0.97 on Lake A, [R5265]). [IL4, T1, L0.]"""
ALBEDO_SNOW_RANGE = (0.8, 0.97)
"""tuple: The range of the snow albedo. [IL4, T1, L0.]"""
ALBEDO_MELTING = 0.4
"""float: The albedo of melting bare ice, 0.2 to 0.6: it falls from 0.78 to 0.60 to 0.44 to 0.20 as the ice warms and the surface melts [R5266]. [IL5, T1, L0.]"""
ALBEDO_MELTING_RANGE = (0.2, 0.6)
"""tuple: The range of the melting ice's albedo. [IL5, T1, L0.]"""
SNOW_COVER_DEPTH_M = 0.03
"""float: The depth of snow at which it has covered the ice under it, as an e-folding length, metres: a layer of 10 cm is 96 percent snow-bright, which is where a snow layer takes
5 percent of the light under it ([R5268]: 10 cm of temperate snow let 5 percent of the visible light through) and the surface is the snow's. [IL6, T7, L1.]"""


def _check_kind(kind):
    if kind not in KINDS:
        raise ValueError("the kind of ice is 'sea' or 'fresh', not %r" % (kind,))


def _check_depth(value, name):
    if not (value >= 0.0 and math.isfinite(value)):
        raise ValueError("%s must be a finite length of 0 or more, not %r" % (name, value))


def albedo(snow_m=0.0, kind=SEA, melting=False, albedo_ice=None, albedo_snow=None):
    """
    The albedo of the lid's top surface to photosynthetic light.

    Bare ice has its own albedo (`ALBEDO_BARE`, or `ALBEDO_MELTING` when `melting`); snow raises it to the snow's. A thin layer of snow is part way: the snow's weight is
    `1 - exp(-snow_m / SNOW_COVER_DEPTH_M)`, so no thickness of snow makes the albedo jump (IL6).
    """
    _check_kind(kind)
    _check_depth(snow_m, "the snow depth")
    ice = albedo_ice if albedo_ice is not None else (ALBEDO_MELTING if melting else ALBEDO_BARE[kind])
    snow = albedo_snow if albedo_snow is not None else ALBEDO_SNOW
    weight = 1.0 - math.exp(-snow_m / SNOW_COVER_DEPTH_M)
    return ice + (snow - ice) * weight


def transmission(ice_m, snow_m=0.0, kind=SEA, melting=False, k_ice=None, k_snow=None, albedo_ice=None, albedo_snow=None):
    """
    The share of the photosynthetic light at the surface that reaches the water under `ice_m` metres of ice and `snow_m` metres of snow: `(1 - albedo) exp(-k_ice ice_m)
    exp(-k_snow snow_m)`. The coefficients default to the central values of the rows; pass any of them to use another (the ends of the ranges, or a world's own).
    A lid of no ice and no snow still reflects (it is the surface of the water, whose albedo is the caller's).
    """
    _check_kind(kind)
    _check_depth(ice_m, "the ice thickness")
    _check_depth(snow_m, "the snow depth")
    ki = k_ice if k_ice is not None else K_ICE_PER_M[kind]
    ks = k_snow if k_snow is not None else K_SNOW_PER_M
    top = albedo(snow_m, kind, melting, albedo_ice, albedo_snow)
    return (1.0 - top) * math.exp(-ki * ice_m - ks * snow_m)


def light_below(surface, ice_m, snow_m=0.0, kind=SEA, melting=False, **coefficients):
    """The light under the lid, in the unit of `surface` (a flux of photosynthetic photons, or watts): `surface` times `transmission`."""
    return surface * transmission(ice_m, snow_m, kind, melting, **coefficients)


def thickest_ice_m(fraction, snow_m=0.0, kind=SEA, melting=False, k_ice=None, k_snow=None, albedo_ice=None, albedo_snow=None):
    """
    The thickest ice, metres, under which `transmission` is still at least `fraction`, with `snow_m` of snow on top. 0 when even the bare water under that snow does not reach `fraction`.
    """
    _check_kind(kind)
    _check_depth(snow_m, "the snow depth")
    if not 0.0 < fraction <= 1.0:
        raise ValueError("the share of the light must be above 0 and at most 1, not %r" % (fraction,))
    ki = k_ice if k_ice is not None else K_ICE_PER_M[kind]
    open_water = transmission(0.0, snow_m, kind, melting, ki, k_snow, albedo_ice, albedo_snow)
    if open_water <= fraction:
        return 0.0
    return math.log(open_water / fraction) / ki


def deepest_snow_m(fraction, ice_m=0.0, kind=SEA, melting=False, k_ice=None, k_snow=None, albedo_ice=None, albedo_snow=None, limit_m=100.0):
    """
    The deepest snow, metres, on `ice_m` of ice under which `transmission` is still at least `fraction` (found by bisection, because the snow also sets the albedo); 0 when the bare ice does not
    reach `fraction`, `limit_m` when no depth up to it stops it.
    """
    _check_kind(kind)
    _check_depth(ice_m, "the ice thickness")
    if not 0.0 < fraction <= 1.0:
        raise ValueError("the share of the light must be above 0 and at most 1, not %r" % (fraction,))

    def t(snow):
        return transmission(ice_m, snow, kind, melting, k_ice, k_snow, albedo_ice, albedo_snow)

    if t(0.0) < fraction:
        return 0.0
    if t(limit_m) >= fraction:
        return limit_m
    low, high = 0.0, limit_m
    for _ in range(80):
        mid = 0.5 * (low + high)
        low, high = (mid, high) if t(mid) >= fraction else (low, mid)
    return low
