# src/planetgen/physics/climate_branches.py

"""
The branches of a climate with an ice-albedo feedback (GEN.1046, W20 of `docs/design/systems/weird-worlds.md`; design in `docs/design/composition/climate-branches.md`).

A world with an ocean cools, freezes, and so reflects more of the light, which cools it further [R0580, R4008]. Written as a fixed point
`T = response(albedo(T))` the feedback can leave three answers for one flux: a warm one, a snowball, and an unstable one between them. This module is pure
numbers: it holds the albedo law of the ice and finds the fixed points of any `response`, which is the surface temperature the rest of the physics gives for a
held Bond albedo (the surface-air fixed point, or a toy in the tests). It reads no body and draws nothing; choosing a branch is the caller's draw.

Stand-ins (registry rows CC15 and CC16 of `gap-filling-plan.md`): the two albedos (CC15) and the ramp: its ends and its exponent (CC16).
"""

from dataclasses import dataclass

from scipy.optimize import brentq

ICE_FREE_ALBEDO = 0.30
"""float: The Bond albedo of a world with no ice (Earth's, 0.30, a test target and a stand-in for ocean and cloud together; row CC15)."""
ICE_COVERED_ALBEDO = 0.62
"""float: The Bond albedo of a world covered in ice and snow, a ramp's cold end (sea ice and snow 0.5 to 0.8, clouds thinner than on a warm world; row CC15)."""
ICE_FREE_ABOVE_K = 305.0
"""float: The mean surface temperature above which no ice is left (a global mean carries the whole spread from the equator to the poles; row CC16)."""
ICE_COVERED_BELOW_K = 240.0
"""float: The mean surface temperature below which the ice reaches the equator (row CC16)."""
ICE_FRACTION_EXPONENT = 3.0
"""float: The ice cover grows as this power of the way down the ramp: slowly while the ice is at the poles, faster as the ice line comes to the tropics, which is what a
one-dimensional ice-line balance does [R3312]. Cubic puts the surface albedo feedback of the Earth at 288 K at 0.34 W/m2 per K, inside the very likely range of
0.10 to 0.60 of IPCC AR6 (2021) table 7.10 (central 0.35). Row CC16."""

SCAN_RANGE_K = (150.0, 700.0)
"""tuple: The temperatures `fixed_points` looks between."""
SCAN_STEP_K = 1.0
"""float: The spacing of the scan; two fixed points closer than this are not told apart."""
SLOPE_STEP_K = 0.05
"""float: The step the slope of the loop is taken over."""


@dataclass(frozen=True)
class IceAlbedoLaw:
    """The Bond albedo as a function of the mean surface temperature: the ice-free value above `free_above_k`, the covered one below `covered_below_k`, and between them
    the ice cover grows as `exponent` power of the way down (1 is a straight line)."""

    ice_free: float = ICE_FREE_ALBEDO
    ice_covered: float = ICE_COVERED_ALBEDO
    free_above_k: float = ICE_FREE_ABOVE_K
    covered_below_k: float = ICE_COVERED_BELOW_K
    exponent: float = ICE_FRACTION_EXPONENT

    def __post_init__(self):
        if not 0.0 <= self.ice_free <= self.ice_covered < 1.0:
            raise ValueError("the albedos must satisfy 0 <= ice-free <= ice-covered < 1")
        if not self.covered_below_k < self.free_above_k:
            raise ValueError("the ice must be covered below a temperature under the one it is free above")
        if not self.exponent > 0.0:
            raise ValueError("the exponent of the ice cover must be positive")

    def ice_fraction(self, temperature_k):
        """The share of the world the ice has covered, 0 (ice-free) to 1 (covered)."""
        share = (self.free_above_k - temperature_k) / (self.free_above_k - self.covered_below_k)
        return min(max(share, 0.0), 1.0) ** self.exponent

    def albedo(self, temperature_k):
        return self.ice_free + (self.ice_covered - self.ice_free) * self.ice_fraction(temperature_k)


@dataclass(frozen=True)
class Branch:
    """One fixed point: the temperature, the albedo it carries, the ice fraction, and whether a small push away from it dies out."""

    temperature_k: float
    albedo: float
    ice_fraction: float
    stable: bool


def fixed_points(response, law=None, *, scan_range_k=SCAN_RANGE_K, step_k=SCAN_STEP_K):
    """
    Every fixed point of `T = response(law.albedo(T))` in `scan_range_k`, coldest first.

    `response(albedo)` is the mean surface temperature, K, of the world when its Bond albedo is held at `albedo`; it must fall as the albedo rises. A fixed point is
    stable when the loop's slope there, `d response(law.albedo(T)) / dT`, is under 1 (the ice-albedo gain is under unity). The coldest and the warmest of a
    bistable climate are stable and the one between them is not.
    """
    law = law or IceAlbedoLaw()

    def residual(temperature_k):
        return temperature_k - response(law.albedo(temperature_k))

    low, high = scan_range_k
    points, previous_t, previous_r = [], low, residual(low)
    steps = int((high - low) / step_k)
    for i in range(1, steps + 1):
        t = low + i * step_k
        r = residual(t)
        if previous_r == 0.0:
            root = previous_t
        elif previous_r * r < 0.0:
            root = brentq(residual, previous_t, t, xtol=1.0e-9)
        else:
            previous_t, previous_r = t, r
            continue
        if not points or root - points[-1].temperature_k > step_k / 2:
            slope = (response(law.albedo(root + SLOPE_STEP_K)) - response(law.albedo(root - SLOPE_STEP_K))) / (2.0 * SLOPE_STEP_K)
            points.append(Branch(root, law.albedo(root), law.ice_fraction(root), slope < 1.0))
        previous_t, previous_r = t, r
    return tuple(points)


def stable_branches(response, law=None, **kwargs):
    """The stable fixed points of `fixed_points`, coldest first (a snowball first, a warm world last)."""
    return tuple(branch for branch in fixed_points(response, law, **kwargs) if branch.stable)


def is_bistable(response, law=None, **kwargs):
    """Whether the climate holds two stable branches at once."""
    return len(stable_branches(response, law, **kwargs)) >= 2


def bistable_window(response_at_flux, fluxes, law=None, **kwargs):
    """
    `(lowest, highest)` of the fluxes in `fluxes` at which the climate holds two stable branches, or `None` if it does at none.

    `response_at_flux(flux)` returns the `response` of a world given that flux (for the same body: the surface-air fixed point with its absorbed light scaled to
    the flux). The window's width, as a share of the middle flux, is what W20 asks to be 10 to 20 percent for an Earth twin.
    """
    inside = [flux for flux in fluxes if is_bistable(response_at_flux(flux), law, **kwargs)]
    return (min(inside), max(inside)) if inside else None
