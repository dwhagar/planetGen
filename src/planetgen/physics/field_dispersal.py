# planetgen/physics/field_dispersal.py

"""
How an interstellar asteroid field spreads, ages and ends (GEN.1100; docs/design/systems/interstellar-asteroid-fields.md section 5, docs/design/research/interstellar-field-sources.md)
=====================================================================================================================================

Pure functions of a seed and an age, no simulation. A field has no star to hold it: every body flies in a straight line from where the field was made, so everything about it is a
formula of its age.

* **Size.** The rms radius of a cloud made from a ball of rms radius `R0` and random velocities of rms speed `v` grows as `sqrt(R0^2 + (v t)^2)` (`rms_radius_m`; the design writes
  `R0 + v t`, the same once `v t` is much more than `R0`). Radius and age are one number: `age_yr_at_radius` is the inverse, so a field's radius can be a result of its age
  (births are uniform in time, `age_from_unit_yr`) and the step is reversible.
* **The end.** Past one light year the cloud is thinner than the interstellar background of 100 m bodies (IF11) and is no longer a field: it is retired (Boss, 4:54 PM EDT 2026-10-30),
  at `retirement_age_yr` (3 Myr for 0.1 km/s, 300 kyr for 1 km/s, 30 kyr for 10 km/s).
* **No grinding.** Collisions stop once the optical depth of the pieces of a size falls under 1, at `freeze_out_radius_m`: 6.5e-4 AU in 5 hours for 100 km pieces of one Earth mass, 0.2 AU in two months for
  metre pieces. Beyond that the size distribution is what the first hours made, so the mass is fixed: `column_density_kg_m2` falls as `1 / t^2` with the mass conserved.
* **Shear.** Galactic shear moves a body at an offset `x` along the orbit by `2 A x t` (`shear_displacement_m`, Oort `A` of `OORT_A_KMS_KPC`); it is a second-order step before 100 Myr or 10 pc,
  and `member_position_m` includes it, valid to `SHEAR_VALID_YEARS` (no epicycle).
* **Bodies.** A body is `(seed, rank)`: its velocity is drawn once from a stream of that pair (so a body keeps its velocity at every age and never depends on the order bodies are made), its position is
  `x0 + v t` plus shear. Two shapes: `ball` (a collision: from a point, isotropic velocities) and `annulus` (a belt unbound by a supernova: from a ring of radius `a` on a straight tangent at the
  speed it leaves with, `unbound_belt_speed_ms`, a thin pancake).

Stand-ins are IF15 to IF18 of `docs/design/research/gap-filling-plan.md` section 6.105. Nothing reads the module yet: the field columns, the host rule and the bring-to-epoch step (GEN.1103,
GEN.1104, GEN.1107) wire it.
"""

import math

from planetgen.physics import constants
from planetgen.util import detmath, draw

RETIREMENT_RADIUS_M = constants.LIGHTYEAR_M
"""float: A field whose rms radius passes this is retired: one light year (Boss, 4:54 PM EDT 2026-10-30; the top of the old radius range, `tuning.ASTEROID_FIELD_RADIUS_RANGE_LY`). IF15, tier T9, L1."""
OORT_A_KMS_KPC = 15.0
"""float: Oort's constant `A`, km/s per kpc, the shear of the disc near the Sun (the design uses 15; 15.3 +/- 0.4 in recent fits, recalled). IF16, tier T2, L0."""
SHEAR_VALID_YEARS = 2.6e7
"""float: The age to which the shearing-sheet drift without epicycles holds, `1 / kappa` with the epicycle frequency of the solar neighbourhood (37 km/s per kpc), years. IF16."""
DEFAULT_DENSITY_KG_M3 = 3000.0
"""float: The bulk density of the pieces, for the freeze-out radius (silicate rubble; the design's table of 5.3 follows from it). IF17, tier T6, L0."""
SHAPES = ("ball", "annulus")
"""tuple: The shape a field starts from."""

_KMS_PER_KPC_TO_PER_S = 1000.0 / (1000.0 * constants.PARSEC_M)
"""float: 1 km/s per kpc in 1/s."""


class FieldDispersalError(ValueError):
    """A radius, speed, mass, age or rank argument is not in range."""


def _positive(**values):
    for name, value in values.items():
        if not value > 0.0:
            raise FieldDispersalError(f"{name} {value!r} is not positive")


def _years_to_s(age_yr):
    if age_yr < 0.0:
        raise FieldDispersalError(f"the age {age_yr!r} is negative")
    return age_yr * constants.SECONDS_PER_YEAR


# --- size and age --------------------------------------------------------------------------------------------------------

def rms_radius_m(age_yr, rms_speed_ms, initial_radius_m=0.0):
    """The rms radius of a field `age_yr` after it was made, `sqrt(R0^2 + (v t)^2)`, m, for pieces of rms speed `rms_speed_ms` (`rogue_collision.unbound_rms_speed_ms`)."""
    _positive(rms_speed_ms=rms_speed_ms)
    if initial_radius_m < 0.0:
        raise FieldDispersalError("the initial radius is negative")
    return detmath.hypot(initial_radius_m, rms_speed_ms * _years_to_s(age_yr))


def age_yr_at_radius(radius_m, rms_speed_ms, initial_radius_m=0.0):
    """The age at which the rms radius is `radius_m`, years: the inverse of `rms_radius_m`. A radius under the initial one is age 0 (a field is not younger than its birth)."""
    _positive(rms_speed_ms=rms_speed_ms)
    if radius_m <= initial_radius_m:
        return 0.0
    return math.sqrt(radius_m ** 2 - initial_radius_m ** 2) / rms_speed_ms / constants.SECONDS_PER_YEAR


def retirement_age_yr(rms_speed_ms, initial_radius_m=0.0, retirement_radius_m=RETIREMENT_RADIUS_M):
    """The age at which a field passes the retirement radius and is retired, years (30 kyr at 10 km/s, 300 kyr at 1 km/s, 3 Myr at 0.1 km/s)."""
    return age_yr_at_radius(retirement_radius_m, rms_speed_ms, initial_radius_m)


def is_retired(age_yr, rms_speed_ms, initial_radius_m=0.0, retirement_radius_m=RETIREMENT_RADIUS_M):
    """Whether a field of this age has passed the retirement radius."""
    return rms_radius_m(age_yr, rms_speed_ms, initial_radius_m) > retirement_radius_m


def age_from_unit_yr(unit, rms_speed_ms, initial_radius_m=0.0, retirement_radius_m=RETIREMENT_RADIUS_M):
    """
    The age of a field drawn at random among those that stand: births are uniform in time, so the age is uniform from 0 to the retirement age, `unit` in `[0, 1)` times it. Its radius is then a
    result of the age (`rms_radius_m`), not a second draw.
    """
    if not 0.0 <= unit < 1.0:
        raise FieldDispersalError(f"the unit draw {unit!r} is outside 0 to 1")
    return unit * retirement_age_yr(rms_speed_ms, initial_radius_m, retirement_radius_m)


def column_density_kg_m2(mass_kg, radius_m):
    """The mass of the field over the disc of its rms radius, `M / (pi R^2)`, kg/m2; falls as `1 / t^2` with the mass fixed (no grinding, no loss)."""
    _positive(mass_kg=mass_kg, radius_m=radius_m)
    return mass_kg / (math.pi * radius_m ** 2)


# --- freeze-out ----------------------------------------------------------------------------------------------------------

def freeze_out_radius_m(total_mass_kg, diameter_m, density_kg_m3=DEFAULT_DENSITY_KG_M3):
    """
    The radius at which the pieces of `diameter_m` stop colliding (optical depth 1): `sqrt(3 M / (2 pi rho D))` for all the mass in pieces of that diameter, m. For one Earth mass
    of 3000 kg/m3: 6.5e-4 AU for 100 km pieces, 6.5e-3 AU for 1 km, 0.2 AU for 1 m, 6.5 AU for 1 mm.
    """
    _positive(total_mass_kg=total_mass_kg, diameter_m=diameter_m, density_kg_m3=density_kg_m3)
    return math.sqrt(3.0 * total_mass_kg / (2.0 * math.pi * density_kg_m3 * diameter_m))


def freeze_out_time_s(total_mass_kg, diameter_m, speed_ms, density_kg_m3=DEFAULT_DENSITY_KG_M3):
    """The time the cloud takes to reach `freeze_out_radius_m` at `speed_ms`, seconds (5 hours for 100 km pieces at 5 km/s, 2 days for 1 km, 2 months for 1 m, 6 years for 1 mm)."""
    _positive(speed_ms=speed_ms)
    return freeze_out_radius_m(total_mass_kg, diameter_m, density_kg_m3) / speed_ms


# --- shear and the belt --------------------------------------------------------------------------------------------------

def shear_displacement_m(offset_m, age_yr):
    """The drift along the orbit of a body held at the offset `offset_m` from the field's centre for `age_yr`, `2 A x t`, m (0.3 pc at 1 pc after 10 Myr)."""
    return 2.0 * OORT_A_KMS_KPC * _KMS_PER_KPC_TO_PER_S * offset_m * _years_to_s(age_yr)


def unbound_belt_speed_ms(progenitor_mass_kg, remnant_mass_kg, belt_radius_m):
    """
    The speed at infinity of a belt body on a circular orbit of radius `belt_radius_m` round a star that loses mass at once to a remnant, `sqrt(G (M0 - 2 M_rem) / a)`, m/s; 0 when the
    remnant holds more than half of the mass (the belt stays bound, IF10). A belt at 50 AU round 10 solar masses leaves a 1.4 solar-mass remnant at 11.3 km/s (13.3 km/s orbital).
    """
    _positive(progenitor_mass_kg=progenitor_mass_kg, remnant_mass_kg=remnant_mass_kg, belt_radius_m=belt_radius_m)
    excess = progenitor_mass_kg - 2.0 * remnant_mass_kg
    return math.sqrt(constants.G * excess / belt_radius_m) if excess > 0.0 else 0.0


# --- bodies --------------------------------------------------------------------------------------------------------------

def member_velocity_ms(seed, rank, rms_speed_ms, shape="ball"):
    """
    The velocity of body `rank` of the field `seed`, m/s, the same at every age and for every call: for a `ball` a normal draw of dispersion `rms_speed_ms / sqrt(3)` in each component, for an
    `annulus` the speed `rms_speed_ms` along the tangent at the body's azimuth (`member_azimuth_rad`).
    """
    _positive(rms_speed_ms=rms_speed_ms)
    if shape not in SHAPES:
        raise FieldDispersalError(f"the shape {shape!r} is none of {SHAPES}")
    if rank < 0:
        raise FieldDispersalError(f"the rank {rank!r} is negative")
    if shape == "annulus":
        phi = member_azimuth_rad(seed, rank)
        return (-rms_speed_ms * math.sin(phi), rms_speed_ms * math.cos(phi), 0.0)
    stream = draw.Stream(f"field-velocity:{seed}:{rank}")
    sigma = rms_speed_ms / math.sqrt(3.0)
    return (stream.gauss(0.0, sigma), stream.gauss(0.0, sigma), stream.gauss(0.0, sigma))


def member_azimuth_rad(seed, rank):
    """The azimuth of body `rank` on its ring, from a stream of `(seed, rank)`, radians in `[0, 2 pi)`."""
    return draw.Stream(f"field-azimuth:{seed}:{rank}").random() * 2.0 * math.pi


def member_position_m(seed, rank, age_yr, rms_speed_ms, shape="ball", belt_radius_m=0.0):
    """
    The position of body `rank` of the field `seed` at `age_yr`, m, from the field's centre: `x0 + v t` plus the shear drift along the orbit, `-2 A (x0 t + v_x t^2 / 2)` (the shift of a body at the offset `x0 + v_x t`,
    integrated). `x0` is the origin for a `ball` and the point of the ring of radius `belt_radius_m` at the body's azimuth for an `annulus`. A pure function of its arguments: not of the order bodies
    were made in, and linear in the age apart from the shear.

    Raises:
        FieldDispersalError: A negative age, rank or radius, an unknown shape, or an annulus with no radius.
    """
    seconds = _years_to_s(age_yr)
    if shape == "annulus":
        _positive(belt_radius_m=belt_radius_m)
        phi = member_azimuth_rad(seed, rank)
        start = (belt_radius_m * math.cos(phi), belt_radius_m * math.sin(phi), 0.0)
    else:
        start = (0.0, 0.0, 0.0)
    vx, vy, vz = member_velocity_ms(seed, rank, rms_speed_ms, shape)
    shear = OORT_A_KMS_KPC * _KMS_PER_KPC_TO_PER_S
    drift = shear * (2.0 * start[0] * seconds + vx * seconds ** 2)       # 2 A (x0 t + v_x t^2 / 2): the along-track shear of a body that is at x0 + v_x t
    return (start[0] + vx * seconds, start[1] + vy * seconds - drift, start[2] + vz * seconds)
