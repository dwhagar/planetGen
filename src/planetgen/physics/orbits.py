# planetgen/physics/orbits.py

"""
Orbits
======

Orbital mechanics for generated systems: the habitable zone, Hill spheres
and mutual Hill radii, the Holman-Wiegert stability limits for binaries,
circular orbital speed, orbital position, the parent's reflex offset and
the shortest orbit update worth storing.
"""

import math

from planetgen import tuning
from planetgen.physics import constants as physical_constants
from planetgen.util.checks import finite_domain


@finite_domain()
def calculate_habitable_zone(luminosity):
    """
    Calculates the inner and outer boundaries of the habitable zone for a star.

    The habitable zone is defined as the region around a star where liquid
    water could exist on a planet's surface. This calculation is based on the
    star's luminosity.

    Args:
        luminosity (float): The luminosity of the star in Watts.

    Returns:
        tuple: A tuple containing the inner and outer radii of the habitable
               zone in AU.
    """
    solar_lum = luminosity / physical_constants.SOLAR_LUMINOSITY
    inner_radius = math.sqrt(solar_lum / 1.1)
    outer_radius = math.sqrt(solar_lum / 0.53)
    return (inner_radius, outer_radius)


@finite_domain()
def calculate_hill_sphere(distance_m, body_mass_kg, central_mass_kg):
    """
    Calculates the Hill sphere radius for a celestial body.

    The Hill sphere is the region around a celestial body where its own
    gravity is the dominant force for attracting satellites. This function
    calculates the radius of this sphere.

    Args:
        distance_m (float): The distance (semi-major axis) between the smaller
                            body and the larger central body, in meters.
        body_mass_kg (float): The mass of the smaller body (e.g., a planet) in kilograms.
        central_mass_kg (float): The mass of the larger central body (e.g., a star) in kilograms.

    Returns:
        float: The radius of the Hill sphere in meters.
    """
    return distance_m * (body_mass_kg / (3 * central_mass_kg)) ** (1 / 3)


@finite_domain(clamped=("companion_mass_fraction", "eccentricity"))
def holman_wiegert_critical_semimajor_axis(binary_separation_au, companion_mass_fraction, eccentricity):
    """
    Holman, M. & Wiegert, P. (1999), AJ 117:621, "Long-Term Stability of
    Planets in Binary Systems" -- the empirical fit for an S-type (wide)
    binary's critical semi-major axis: the largest orbit around ONE star of
    the pair that remains long-term stable against the other star's
    periodic gravitational perturbation.

        a_crit / a_bin = 0.464 - 0.380*mu - 0.631*e + 0.586*mu*e
                         + 0.150*e^2 - 0.198*mu*e^2

    `mu` is the *perturbing companion's* mass fraction of the pair's total
    mass, `mu = M_companion / (M_this_star + M_companion)` -- this is
    evaluated once per star, using that star's own companion, so calling
    this twice for one pair (once from each star's perspective) generally
    yields two different `a_crit` values unless the two masses are equal.

    Valid over roughly `mu` in [0.1, 0.9] and `e` in [0.0, 0.8] (Holman &
    Wiegert's own numerical grid doesn't extend meaningfully further) --
    inputs are clamped to `physical_constants.HOLMAN_WIEGERT_MU_RANGE`/
    `HOLMAN_WIEGERT_ECCENTRICITY_RANGE` rather than extrapolated, since a
    saturated-at-the-boundary estimate is more useful than either an
    exception or a silently invalid extrapolation.

    Args:
        binary_separation_au (float): The binary pair's own semi-major
                                      axis (separation), in AU.
        companion_mass_fraction (float): `mu`, as defined above (0-1).
        eccentricity (float): The binary orbit's own eccentricity (0-1).

    Returns:
        float: The critical semi-major axis, in AU (same unit as
              `binary_separation_au`) -- this star's own maximum stable
              planetary orbit given the companion's influence.

    Example (equal-mass, circular pair -- a standard reference case):
        `holman_wiegert_critical_semimajor_axis(1.0, 0.5, 0.0)` gives
        `0.464 - 0.380*0.5 = 0.274` exactly, consistent with the commonly
        cited ~0.27-0.30 * a_bin figure for this case.
    """
    mu_min, mu_max = physical_constants.HOLMAN_WIEGERT_MU_RANGE
    e_min, e_max = physical_constants.HOLMAN_WIEGERT_ECCENTRICITY_RANGE
    mu = min(max(companion_mass_fraction, mu_min), mu_max)
    e = min(max(eccentricity, e_min), e_max)

    ratio = (0.464 - 0.380 * mu - 0.631 * e
             + 0.586 * mu * e + 0.150 * e ** 2 - 0.198 * mu * e ** 2)
    return ratio * binary_separation_au



@finite_domain(clamped=("secondary_mass_fraction", "eccentricity"))
def holman_wiegert_circumbinary_a_crit_au(binary_separation_au, secondary_mass_fraction, eccentricity=0.0):
    """
    Holman, M. & Wiegert, P. (1999), AJ 117:621 -- the empirical fit for a
    P-type (circumbinary) orbit's critical semi-major axis: the SMALLEST
    orbit around both stars of a close pair that stays long-term stable.

        a_crit / a_bin = 1.60 + 5.10*e - 2.22*e^2 + 4.12*mu - 4.27*e*mu
                         - 5.09*mu^2 + 4.61*e^2*mu^2

    `mu` is the lighter star's fraction of the pair's total mass. Inputs
    are clamped to `tuning.HOLMAN_WIEGERT_P_TYPE_MU_RANGE`/
    `HOLMAN_WIEGERT_P_TYPE_ECCENTRICITY_RANGE` (the fit's tested grid)
    rather than extrapolated. A circular equal-mass pair gives about
    2.39 * a_bin.

    Args:
        binary_separation_au (float): The pair's own separation, in AU.
        secondary_mass_fraction (float): `mu`, as defined above.
        eccentricity (float): The pair's orbital eccentricity.

    Returns:
        float: The innermost stable circumbinary orbit, in AU.
    """
    mu_min, mu_max = tuning.HOLMAN_WIEGERT_P_TYPE_MU_RANGE
    e_min, e_max = tuning.HOLMAN_WIEGERT_P_TYPE_ECCENTRICITY_RANGE
    mu = min(max(secondary_mass_fraction, mu_min), mu_max)
    e = min(max(eccentricity, e_min), e_max)

    ratio = (1.60 + 5.10 * e - 2.22 * e ** 2 + 4.12 * mu - 4.27 * e * mu
             - 5.09 * mu ** 2 + 4.61 * e ** 2 * mu ** 2)
    return ratio * binary_separation_au

@finite_domain()
def mutual_hill_radius_au(mass1_kg, mass2_kg, distance1_au, distance2_au, central_mass_kg):
    """
    Gladman (1993), Icarus 106:247, "Dynamical stability of the outer solar
    system and the delivery of comets" -- the mutual Hill radius of two
    orbiting bodies:

        R_H,mutual = ((m1 + m2) / (3 * M_central))^(1/3) * ((a1 + a2) / 2)

    Gladman's own derivation assumes both bodies orbit the SAME central
    mass -- used here (see `systemData.StarSystem._validate_cross_star_clearance`)
    across two planets that orbit *different* stars of a wide binary, this
    is a physically-motivated extension of the criterion's spirit rather
    than a literal textbook application; see that method's docstring for
    the specific choice of `central_mass_kg` and its justification.

    Args:
        mass1_kg (float): First body's own mass, in kg.
        mass2_kg (float): Second body's own mass, in kg.
        distance1_au (float): First body's distance from whatever it
                              orbits, in AU.
        distance2_au (float): Second body's distance from whatever it
                              orbits, in AU.
        central_mass_kg (float): The mass, in kg, both distances above are
                                 measured against (see docstring above for
                                 how this generator chooses it when the two
                                 bodies orbit different stars).

    Returns:
        float: The mutual Hill radius, in AU (same unit as
              `distance1_au`/`distance2_au`).
    """
    return ((mass1_kg + mass2_kg) / (3 * central_mass_kg)) ** (1 / 3) * ((distance1_au + distance2_au) / 2)


@finite_domain()
def mutual_hill_radius_m(distance1_m, distance2_m, mass1_kg, mass2_kg, central_mass_kg):
    """
    Calculates the *mutual* Hill radius of two bodies that orbit the same
    primary -- the length scale real orbital-dynamics stability criteria
    (Gladman 1993; Chambers, Wetherill & Boslough 1996; Smith & Lissauer
    1999/2009) use for how close two adjacent planets' orbits can safely
    be, as distinct from `calculate_hill_sphere`'s single-body sphere of
    gravitational dominance (the right tool for "how far can a satellite
    orbit *this* body," not for "how close can two planets orbit each
    other").

    R_H,mutual = ((a1 + a2) / 2) * ((m1 + m2) / (3 * M_central)) ** (1/3)

    -- i.e. the single-body formula generalized to use the *pair's*
    combined mass and average distance, rather than either body's own
    mass and distance alone. See
    `tuning.MUTUAL_HILL_RADII_SEPARATION` for how this
    generator turns this length into an actual minimum separation.

    Args:
        distance1_m (float): The first body's distance from the shared
                             primary, in meters.
        distance2_m (float): The second body's distance from the shared
                             primary, in meters.
        mass1_kg (float): The first body's mass, in kilograms.
        mass2_kg (float): The second body's mass, in kilograms.
        central_mass_kg (float): The shared primary's mass, in kilograms.

    Returns:
        float: The mutual Hill radius, in meters.
    """
    avg_distance_m = (distance1_m + distance2_m) / 2
    return avg_distance_m * ((mass1_kg + mass2_kg) / (3 * central_mass_kg)) ** (1 / 3)


@finite_domain()
def circular_orbital_speed_kms(distance_au, period_years):
    """
    Tangential speed of a circular orbit, given its radius and period:
    `v = 2*pi*r / T`. Exact (constant at every point of the orbit), since
    this generator only ever models circular orbits for planets and moons
    (see `planetPhysics.generate_orbital_motion_properties`) -- unlike
    `calculate_galactic_orbit`, which has to *assume* a rotation-curve
    model to get a speed at all, a planet's/moon's period is already known
    exactly from Kepler's third law (`planetPhysics.
    calculate_orbital_period_years`), so speed here is a direct
    geometric consequence of the two, not a separate physical model.

    Args:
        distance_au (float): Orbital radius (semi-major axis), in AU.
        period_years (float): Orbital period, in years.

    Returns:
        float: Orbital speed, in km/s.
    """
    circumference_km = 2 * math.pi * distance_au * physical_constants.AU_TO_KM
    period_seconds = period_years * physical_constants.SECONDS_PER_YEAR
    return circumference_km / period_seconds


@finite_domain()
def orbital_position_au(distance_au, inclination_deg, ascending_node_deg, phase_deg):
    """
    Converts a circular orbit's elements -- radius, inclination, ascending
    node, and current phase -- into a 3D Cartesian position relative to
    the orbit's primary (the body actually being orbited: a star/binary
    system center for a planet, a planet for a moon).

    Standard orbital-plane-to-reference-frame rotation, specialized for a
    circular orbit: `orbital_phase_deg` already plays the role of the
    argument of latitude `u = omega + true_anomaly` directly (no separate
    argument-of-periapsis term, since a circular orbit has no periapsis to
    measure one from -- see `planetPhysics.generate_orbital_motion_properties`'s
    docstring). `inclination_deg`/`ascending_node_deg` orient the orbital
    plane itself; `phase_deg` is where the body sits within it:

        u = radians(phase_deg), i = radians(inclination_deg), Om = radians(ascending_node_deg)
        x = r * (cos(Om)*cos(u) - sin(Om)*sin(u)*cos(i))
        y = r * (sin(Om)*cos(u) + cos(Om)*sin(u)*cos(i))
        z = r * sin(u)*sin(i)

    At `i = 0` (an uninclined orbit) this correctly collapses to a flat
    circle in the primary's own reference plane (`z = 0` always); `node`
    and `phase` become degenerate there (no inclined plane left for the
    ascending node to describe the crossing of), so only their sum
    matters: `(r*cos(node+u), r*sin(node+u), 0)` -- the standard,
    physically expected behavior at this degenerate case, same as in real
    orbital mechanics, not a bug.

    Args:
        distance_au (float): Orbital radius, in AU.
        inclination_deg (float): Orbital plane tilt, in degrees, relative
                                 to the primary's reference plane.
        ascending_node_deg (float): Longitude of the ascending node, in
                                    degrees.
        phase_deg (float): Current argument of latitude (position angle
                           around the orbit), in degrees.

    Returns:
        tuple: `(x_au, y_au, z_au)`, relative to the primary, in the same
              reference frame `inclination_deg`/`ascending_node_deg` are
              measured against.
    """
    u = math.radians(phase_deg)
    i = math.radians(inclination_deg)
    node = math.radians(ascending_node_deg)

    cos_u, sin_u = math.cos(u), math.sin(u)
    cos_i = math.cos(i)
    cos_node, sin_node = math.cos(node), math.sin(node)

    x = distance_au * (cos_node * cos_u - sin_node * sin_u * cos_i)
    y = distance_au * (sin_node * cos_u + cos_node * sin_u * cos_i)
    z = distance_au * sin_u * math.sin(i)

    return x, y, z


@finite_domain()
def circular_orbital_velocity_au_per_year(distance_au, inclination_deg, ascending_node_deg, phase_deg, period_years):
    """
    The velocity, relative to the primary, of a body on the circular orbit
    `orbital_position_au` places it on: the position's derivative with the
    phase, times the phase's rate of change (360 degrees a period, in the
    direction of increasing phase). Same axes as `orbital_position_au`.

    Args:
        distance_au (float): Orbital radius, in AU.
        inclination_deg (float): Orbital plane tilt, in degrees.
        ascending_node_deg (float): Longitude of the ascending node, in degrees.
        phase_deg (float): Current argument of latitude, in degrees.
        period_years (float): Orbital period, in years. Positive.

    Returns:
        tuple: `(vx, vy, vz)` in AU per year.
    """
    if period_years <= 0:
        raise ValueError(f"the orbital period must be positive, got {period_years!r}")
    u = math.radians(phase_deg)
    i = math.radians(inclination_deg)
    node = math.radians(ascending_node_deg)
    cos_u, sin_u = math.cos(u), math.sin(u)
    cos_i = math.cos(i)
    cos_node, sin_node = math.cos(node), math.sin(node)
    scale = distance_au * 2 * math.pi / period_years
    return (
        scale * (-cos_node * sin_u - sin_node * cos_u * cos_i),
        scale * (-sin_node * sin_u + cos_node * cos_u * cos_i),
        scale * cos_u * math.sin(i),
    )


def calculate_reflex_offset(parent_mass_kg, children):
    """
    A parent body's own displacement from its nominal fixed point, caused
    by the combined gravitational pull of everything orbiting it -- the
    "wobble"/reflex-motion half of a proper two-body (barycentric)
    treatment, mirrored on the SQL side by the correlated-subquery
    `UPDATE`s in `_db.advance_orbital_phases`.

    For a single child, this is the exact two-body barycentric formula:
    `offset = -(child_mass / (parent_mass + child_mass)) * relative_vector`,
    where `relative_vector` is the child's own already-stored position
    relative to the parent (e.g. `Planet.position_x/y/z`,
    `BinaryStarProxy.binary_mutual_position_x/y/z`) -- deliberately never
    recomputed or changed by this function, since a large amount of
    existing physics (insolation, Hill sphere, tidal locking) depends on
    that vector remaining the *true* separation, not a barycenter-reduced
    one; this only computes the *parent's* own small displacement.

    For multiple children (a star with several planets, a planet with
    several moons), each child's individual pairwise pull is summed --
    the standard linear-superposition approximation real radial-velocity
    work uses for multi-planet reflex motion, exact to leading order
    whenever the parent is much more massive than any single child (true
    for every relationship in this generator).

    Args:
        parent_mass_kg (float): The parent body's own mass, in kg.
        children (list): `(mass_kg, x_au, y_au, z_au)` tuples, one per
                         body orbiting the parent, each already in the
                         parent's own reference frame.

    Returns:
        tuple: `(x_au, y_au, z_au)`, the parent's own offset from its
              nominal fixed point -- `(0.0, 0.0, 0.0)` if `children` is
              empty.
    """
    offset_x = offset_y = offset_z = 0.0
    for child_mass_kg, x_au, y_au, z_au in children:
        mu_child = child_mass_kg / (parent_mass_kg + child_mass_kg)
        offset_x -= mu_child * x_au
        offset_y -= mu_child * y_au
        offset_z -= mu_child * z_au
    return offset_x, offset_y, offset_z


@finite_domain()
def minimum_update_interval_years(period_years):
    """
    The shortest `elapsed_years` worth advancing a body's
    `orbital_phase_deg` for at all -- below this, the phase delta added is
    smaller than `orbital_phase_deg`'s own floating-point resolution, so
    `MOD(orbital_phase_deg + delta, 360)` is guaranteed to round right
    back to the exact value already stored: a wasted write that changes
    nothing. Exists specifically to guard `_db.advance_orbital_phases`
    against that silent no-op, not as a narrative/display stat -- unlike
    this package's other derived quantities, there is no scale at which a
    human would want to read this number (it lands in the nanosecond
    range for any realistic orbital period, see below).

    Derivation: `orbital_phase_deg` ranges over `[0, 360)`, stored as an
    IEEE 754 double (MySQL `DOUBLE`, Python `float` -- identical
    representation). The coarsest (least precise) representable step
    anywhere in that range is the unit-in-the-last-place at magnitudes
    just under 360 -- `math.ulp(360.0)` -- used here as a single,
    domain-wide conservative bound rather than a per-row value that would
    depend on the body's current phase (tighter near 0, coarser near 360)
    and so would itself need updating every time phase does, for no real
    benefit. A phase delta at or above this many degrees is guaranteed to
    change the stored value, regardless of where in `[0, 360)` the
    current phase happens to sit; anything smaller might not.

    `elapsed_years / period_years * 360 >= ulp_deg`
    `elapsed_years >= period_years * ulp_deg / 360`

    Worked example: a 1-year period gives a floor around 5e-9 seconds --
    roughly 16 orders of magnitude below `planetgen.cli.orbits`'s own "once a
    month or so" real-world cadence (see that module's docstring), so
    this guard exists for correctness/defensiveness (a future caller
    advancing time in much smaller steps, e.g. a fast-forward simulation)
    rather than because today's actual usage pattern ever comes close to
    triggering it.

    Args:
        period_years (float): The body's own orbital period, in years.
                              Always positive and finite for a real
                              generated planet/moon (Kepler's third law on
                              a positive distance and mass).

    Returns:
        float: The minimum `elapsed_years` worth calling
              `_db.advance_orbital_phases` for, for a body with this
              period.
    """
    ulp_deg = math.ulp(360.0)
    return period_years * ulp_deg / 360
