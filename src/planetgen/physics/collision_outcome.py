# planetgen/physics/collision_outcome.py

"""
What a collision does (GEN.345; docs/design/orbits/collisions-and-mergers.md sections 4 and 5.1)
=========================================================================================

One pure function, `resolve`, from two bodies that have just come within `R1 + R2` of each other (found by
`physics/closest_approach.py`) to an `Outcome`: the code of the row of the section 5.1 table, the bodies that survive and
what they became, the debris field, the mass a gas giant lost, the numbers the choice was made from and the flags the
admin report carries. It reads no database, draws no random number and moves nothing; the writer (GEN.197, `record_event`)
and the products' generators (GEN.110) take the `Outcome` from here.

The pair is put in key order first, so `resolve(a, b, ...)` and `resolve(b, a, ...)` agree whatever order a broad phase
found them in, and the same input gives the same output. `Outcome.seed` is `collision:<key a>:<key b>:<window start>`: the
part of the draw stream's name that depends on the pair (`draw.Stream(f"{galaxy_seed}:{outcome.seed}")`), for the
products' generators to continue from. Nothing here needs a draw: every row of the table is a function of the numbers.

**The numbers** (section 4.1), all SI. `v_inf` is the relative speed of the stored velocities, `v_esc` the mutual escape
speed at contact, `v_imp = sqrt(v_inf^2 + v_esc^2)`; `Q_R = 0.5 mu v_imp^2 / M_tot` the specific impact energy and
`Q*_RD` the catastrophic disruption threshold of Leinhardt and Stewart (2012) for a strengthless planet (c* = 1.9,
mu_bar = 0.36), at the pair's mass ratio; `x = Q_R / Q*_RD`. A pass whose impact parameter is beyond `b_crit = R_target /
(R_target + R_projectile)` only overlaps a lens of the projectile, and `Q_R` and `Q*_RD` are then worked out for the
interacting mass, the projectile's lens fraction `(3 R_p l^2 - l^3) / (4 R_p^3)` of its mass, `l = R_t + R_p - b`.

**The rows** (section 5.1), by the kinds T (rocky or icy, below the gas giant threshold), G (gas giant), B (brown dwarf),
S (star) and X (neutron star or black hole):

- T + T: a merger (`b < b_crit` and `x < 0.1`, or `v_imp < 1.2 v_esc`), a graze-and-merge (`b >= b_crit` and `v_imp <=
  1.5 v_esc`), a hit-and-run (`b >= b_crit`, `v_imp > 1.5 v_esc`, `x < 1`: both survive, each stripped of half the eroded
  mass `x/2` of an interacting mass of twice the projectile's lens, as a field when that is at least 1 percent of the
  total, else only lost), a partial erosion (`b < b_crit`, `0.1 <= x < 1`: one body of `M_tot (1 - x/2)`, the rest a field)
  and a disruption (`x >= 1`: both destroyed, one field of `M_tot` with a hero body of the largest remnant). With
  `ALWAYS_DESTROY_TERRESTRIAL` every contact is the last, with no hero body.
- T or G or B + S or X (or the fluid and rigid Roche test on X): absorbed by the larger kind (S or X); a body whose pass
  is inside the Roche limit of a neutron star or black hole (1.26 rigid, 2.44 fluid, section 2.3) is a `tidal_disruption`;
  no debris either way.
- T + G or B: absorbed by the giant, which grows by the mass and keeps its momentum.
- G or B + G or B: a `giant_merger` that loses `f = clamp(Q_R / (2 Q*_RD), 0.01, 0.5)` of its mass as gas (no field),
  keeps the momentum of the rest, takes its radius from `giant_radius_km` and is flagged fast spinning. The result is
  then sorted by the one fusion function of metallicity (section 4.3): under the deuterium limit a gas giant, to the
  hydrogen-burning limit `hydrogen_burning_limit_sol(z)` a brown dwarf (`brown_dwarf_formed` when neither parent was
  one), at or above it a new star (`star_created`: no planets, age 0, "merger remnant"). A result within 3 percent of the
  hydrogen limit is flagged borderline, and one between it and the generator's 0.08 floor is flagged.
- S or X + S or X: a `stellar_merger` that conserves momentum and leaves the star's own properties to its generator.
- Not in the table, so a design extension: a G or B body with S or X is absorbed (as a T is).

Radii: T bodies add by volume; a giant's follows `planets.giant_radius_km` (so does a brown dwarf's). Stars and compact
remnants have `radius_m = None` in the outcome: the generator that makes the product works it out.

Conservation: `conservation_error` gives the relative mass and momentum error of an outcome, the gas a giant lost counted
at the pair's centre-of-mass velocity.
"""

import dataclasses
import math

from planetgen import tuning
from planetgen.physics import constants
from planetgen.physics import planets
from planetgen.physics import relations
from planetgen.physics import stellar_abundances, unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")
_mjup_to_kg = units.converter("Mjup", "kg")

KINDS = ("T", "G", "B", "S", "X")
"""tuple: The kinds of body section 5.1 names."""

CODES = ("merger", "graze_and_merge", "hit_and_run", "erosion", "disruption", "absorbed", "tidal_disruption",
         "giant_merger", "stellar_merger")
"""tuple: The rows of the table a collision can end in (`Outcome.code`)."""

RIGID_ROCHE_FACTOR = 1.26
FLUID_ROCHE_FACTOR = 2.44
"""float: The Roche coefficients (section 2.3): rocky and icy bodies, and gas giants, brown dwarfs and stars."""

RHO_1 = 1000.0
"""float: The reference density, kg/m^3, of the catastrophic disruption threshold (section 4.1)."""

CRITICAL_COEFFICIENT, CRITICAL_MU_BAR = 1.9, 0.36
"""float: c* and mu_bar for a strengthless planet. The small-body set (5 and 0.37) is not used: no kind here is one."""

HERO_FADE_X = 1.8
"""float: The x past which the largest remnant is `0.1 M_tot (x / 1.8)^-1.5` (section 4.1); between 1 and this it runs
straight from 0.5 to 0.1 (a stand-in, the source gives no value there)."""


@dataclasses.dataclass(frozen=True)
class Body:
    """
    One body at the impact. `position_m` and `velocity_ms` are in the galactic axes (any common origin; only differences
    and the centre of mass use them), `composition` is `stony`, `carbonaceous` or `icy`, `metallicity_z` the mass share of
    metals the fusion limits read, and `key` an orderable name (the object ID).
    """

    key: str
    kind: str
    mass_kg: float
    radius_m: float
    position_m: tuple = (0.0, 0.0, 0.0)
    velocity_ms: tuple = (0.0, 0.0, 0.0)
    composition: str = "stony"
    metallicity_z: float = stellar_abundances.PROTO_Z


@dataclasses.dataclass(frozen=True)
class Survivor:
    """A body that is still there: `key` is the larger parent's when two made one, `radius_m` is `None` where the
    product's own generator sets it (a star or a compact remnant)."""

    key: str
    kind: str
    mass_kg: float
    radius_m: float | None
    position_m: tuple
    velocity_ms: tuple
    composition: str
    flags: tuple = ()


@dataclasses.dataclass(frozen=True)
class Debris:
    """The asteroid field an event makes: its mass, centre-of-mass position and velocity, starting radius (`R1 + R2`),
    composition (`mixed` when the parents differ) and the mass of its largest fragment (`hero_mass_kg`, 0 for none)."""

    mass_kg: float
    position_m: tuple
    velocity_ms: tuple
    radius_m: float
    composition: str
    hero_mass_kg: float


@dataclasses.dataclass(frozen=True)
class Outcome:
    """
    What the collision did. `code` is one of `CODES`; `records` the event records the table names (the code, then
    `brown_dwarf_formed` or `star_created` for a giant merger that crossed a fusion line); `destroyed` the keys that are
    gone; `survivors` the bodies left (in key order); `debris` the field, or `None`; `lost_mass_kg` the mass that left
    with no field (a giant's gas, a small hit-and-run strip, mass accreted by a compact object counts in its survivor);
    `numbers` the quantities the choice was made from; `flags` what the report carries.
    """

    code: str
    records: tuple
    destroyed: tuple
    survivors: tuple
    debris: Debris | None
    lost_mass_kg: float
    numbers: dict
    flags: tuple
    seed: str


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _weighted(m1, v1, m2, v2):
    total = m1 + m2
    return tuple((m1 * x + m2 * y) / total for x, y in zip(v1, v2))


def _norm(v):
    return math.sqrt(sum(c * c for c in v))


def escape_speed(m_total, r_sum):
    """The mutual escape speed at contact, m/s: `sqrt(2 G (M1 + M2) / (R1 + R2))`."""
    return relations.escape_speed(m_total, r_sum)


def catastrophic_energy(m_target, m_projectile):
    """
    `Q*_RD`, J/kg, the specific impact energy at which the largest remnant is half the total mass (Leinhardt and Stewart
    2012): `[(gamma + 1)^2 / (4 gamma)]^(2 / (3 mu_bar) - 1) c* (4/5) pi rho_1 G R_C1^2` for `gamma = M_p / M_t` and `R_C1`
    the radius of the total mass at `RHO_1`. Earth + Earth: 6.4e7.
    """
    if not (m_target > 0.0 and m_projectile > 0.0) or m_projectile > m_target:
        raise ValueError("the target is the larger mass and both are positive")
    gamma = m_projectile / m_target
    r_c1 = (3.0 * (m_target + m_projectile) / (4.0 * math.pi * RHO_1)) ** (1.0 / 3.0)
    base = CRITICAL_COEFFICIENT * 0.8 * math.pi * RHO_1 * constants.G * r_c1 ** 2
    return ((gamma + 1.0) ** 2 / (4.0 * gamma)) ** (2.0 / (3.0 * CRITICAL_MU_BAR) - 1.0) * base


def lens_fraction(r_projectile, depth):
    """The share of a projectile's mass in the lens `depth` deep that overlaps the target: `(3 R l^2 - l^3) / (4 R^3)`,
    0 for no overlap and 1 for the whole body."""
    depth = min(max(depth, 0.0), 2.0 * r_projectile)
    return (3.0 * r_projectile * depth ** 2 - depth ** 3) / (4.0 * r_projectile ** 3)


def largest_remnant_fraction(x):
    """`M_lr / M_tot` for `x = Q_R / Q*_RD`: `1 - x/2` below 1, straight to 0.1 at `HERO_FADE_X`, then
    `0.1 (x / 1.8)^-1.5` (section 4.1)."""
    if x < 1.0:
        return 1.0 - x / 2.0
    if x < HERO_FADE_X:
        return 0.5 - 0.4 * (x - 1.0) / (HERO_FADE_X - 1.0)
    return 0.1 * (x / HERO_FADE_X) ** -1.5


def _volume_radius(*radii):
    return sum(r ** 3 for r in radii) ** (1.0 / 3.0)


def _giant_radius_m(mass_kg):
    return planets.giant_radius_km(mass_kg) * 1000.0


def _field_composition(a, b):
    return a.composition if a.composition == b.composition else "mixed"


def _centre(a, b):
    return _weighted(a.mass_kg, a.position_m, b.mass_kg, b.position_m)


def _velocity(a, b):
    return _weighted(a.mass_kg, a.velocity_ms, b.mass_kg, b.velocity_ms)


def _survivor(body, mass_kg, radius_m, position_m, velocity_ms, flags=(), kind=None, composition=None):
    return Survivor(body.key, kind or body.kind, mass_kg, radius_m, position_m, velocity_ms,
                    composition or body.composition, tuple(flags))


def _fusion(survivor, z, parents, flags, records):
    mass_sol = _kg_to_msun(survivor.mass_kg)
    deuterium_kg = _mjup_to_kg(stellar_abundances.deuterium_burning_limit_mj(z))
    hydrogen_sol = stellar_abundances.hydrogen_burning_limit_sol(z)
    if abs(mass_sol / hydrogen_sol - 1.0) <= tuning.COLLISION_BORDERLINE_SHARE:
        flags.append("borderline")
    if mass_sol >= hydrogen_sol:
        records.append("star_created")
        flags.extend(("merger remnant", "no planets"))
        if mass_sol < tuning.IMF_BREAKS_SOL[0]:
            flags.append("below the generator's floor")
        return dataclasses.replace(survivor, kind="S", radius_m=None, composition=survivor.composition)
    if survivor.mass_kg >= deuterium_kg:
        if all(p.kind != "B" for p in parents):
            records.append("brown_dwarf_formed")
        return dataclasses.replace(survivor, kind="B")
    return dataclasses.replace(survivor, kind="G")


def _roche(compact, body):
    factor = RIGID_ROCHE_FACTOR if body.kind == "T" else FLUID_ROCHE_FACTOR
    rho_compact = relations.mean_density(compact.mass_kg, compact.radius_m)
    rho_body = relations.mean_density(body.mass_kg, body.radius_m)
    return factor * compact.radius_m * (rho_compact / rho_body) ** (1.0 / 3.0)


def resolve(a, b, impact_parameter_m, window_start_unix, always_destroy_terrestrial=None):
    """
    The outcome of bodies `a` and `b` having come within `R1 + R2`.

    Args:
        a, b (Body): The two bodies at the impact, in either order.
        impact_parameter_m (float): Their closest approach along the straight relative path, metres, `0 <= b <= R1 + R2`
            (`closest_approach.sweep`'s miss distance).
        window_start_unix (float): The start of the window the hit was found in; part of the seed.
        always_destroy_terrestrial (bool, optional): Boss's literal rule for T + T (default `tuning.ALWAYS_DESTROY_TERRESTRIAL`).

    Returns:
        Outcome: See the module docstring.

    Raises:
        ValueError: For an unknown kind, a mass or radius that is not positive and finite, equal keys, or an impact
            parameter outside `[0, R1 + R2]`.
    """
    for body in (a, b):
        if body.kind not in KINDS:
            raise ValueError(f"unknown kind {body.kind!r}")
        if not (math.isfinite(body.mass_kg) and body.mass_kg > 0.0 and math.isfinite(body.radius_m) and body.radius_m > 0.0):
            raise ValueError(f"{body.key}: mass and radius must be positive and finite")
    if a.key == b.key:
        raise ValueError("a body does not collide with itself")
    if a.key > b.key:
        a, b = b, a
    r_sum = a.radius_m + b.radius_m
    if not (math.isfinite(impact_parameter_m) and 0.0 <= impact_parameter_m <= r_sum * (1.0 + 1e-12)):
        raise ValueError(f"the impact parameter {impact_parameter_m!r} is outside 0 to R1 + R2 = {r_sum!r}")
    destroy = tuning.ALWAYS_DESTROY_TERRESTRIAL if always_destroy_terrestrial is None else always_destroy_terrestrial
    seed = f"collision:{a.key}:{b.key}:{window_start_unix}"

    target, projectile = (a, b) if (a.mass_kg, b.key) >= (b.mass_kg, a.key) else (b, a)
    m_total = a.mass_kg + b.mass_kg
    v_inf = _norm(tuple(x - y for x, y in zip(a.velocity_ms, b.velocity_ms)))
    v_esc = escape_speed(m_total, r_sum)
    v_imp = math.sqrt(v_inf ** 2 + v_esc ** 2)
    sine = min(impact_parameter_m / r_sum, 1.0)
    b_crit = target.radius_m / r_sum
    grazing = sine >= b_crit
    depth = r_sum - impact_parameter_m
    lens = lens_fraction(projectile.radius_m, depth) if grazing else 1.0
    m_proj = projectile.mass_kg * lens
    if m_proj > 0.0:
        m_tot_eff = target.mass_kg + m_proj
        q_r = 0.5 * (target.mass_kg * m_proj / m_tot_eff) * v_imp ** 2 / m_tot_eff
        q_star = catastrophic_energy(target.mass_kg, m_proj)
        x = q_r / q_star
    else:
        q_r = q_star = x = 0.0
    numbers = {"v_inf": v_inf, "v_esc": v_esc, "v_imp": v_imp, "q_r": q_r, "q_star": q_star, "x": x, "b": sine,
               "b_crit": b_crit, "lens_fraction": lens, "total_mass_kg": m_total}
    centre, v_com = _centre(a, b), _velocity(a, b)
    flags, records = [], []

    def finish(code, survivors, destroyed, debris=None, lost=0.0, extra=()):
        flags.extend(extra)
        return Outcome(code, (code, *records), tuple(destroyed), tuple(survivors), debris, lost, numbers, tuple(flags), seed)

    def merged_into_larger(mass_kg, radius_m, kind=None, composition=None, extra=()):
        keeper = target
        gone = projectile.key
        return _survivor(keeper, mass_kg, radius_m, centre, v_com, extra, kind, composition), gone

    kinds = {a.kind, b.kind}
    compact = [body for body in (a, b) if body.kind == "X"]
    lighter = [body for body in (a, b) if body.kind != "X"]

    # --- A compact remnant against anything small enough to be torn up first -----------------------------------------
    if len(compact) == 1 and lighter[0].kind in ("T", "G", "B"):
        x_body, other = compact[0], lighter[0]
        inside = impact_parameter_m <= _roche(x_body, other)
        survivor = _survivor(x_body, m_total, None, centre, v_com, ("accreted",), composition=x_body.composition)
        return finish("tidal_disruption" if inside else "absorbed", [survivor], [other.key])

    if kinds & {"S", "X"}:
        if kinds <= {"S", "X"}:
            survivor, gone = merged_into_larger(m_total, None, extra=("outside the table: merged conserving momentum",))
            return finish("stellar_merger", [survivor], [gone])
        big = target if target.kind in ("S", "X") else projectile
        small = projectile if big is target else target
        extra = ("transient",) if big.kind == "S" else ()
        survivor = _survivor(big, m_total, None, centre, v_com, extra)
        return finish("absorbed", [survivor], [small.key])

    # --- Gas giants and brown dwarfs ---------------------------------------------------------------------------------
    if "G" in kinds or "B" in kinds:
        if "T" in kinds:
            giant = a if a.kind != "T" else b
            small = b if giant is a else a
            new = _survivor(giant, m_total, _giant_radius_m(m_total), centre, v_com)
            new = _fusion(new, giant.metallicity_z, (giant,), flags, records)
            return finish("absorbed", [new], [small.key])
        loss = min(max(q_r / (2.0 * q_star), tuning.GIANT_MERGER_LOSS_RANGE[0]), tuning.GIANT_MERGER_LOSS_RANGE[1])
        new_mass = (1.0 - loss) * m_total
        z = (a.mass_kg * a.metallicity_z + b.mass_kg * b.metallicity_z) / m_total
        new = _survivor(target, new_mass, _giant_radius_m(new_mass), centre, v_com,
                        ("rapidly rotating", "inflated", "hot"))
        new = _fusion(new, z, (a, b), flags, records)
        return finish("giant_merger", [new], [projectile.key], lost=loss * m_total)

    # --- Two terrestrial bodies -------------------------------------------------------------------------------------
    composition = _field_composition(a, b)
    merged_radius = _volume_radius(a.radius_m, b.radius_m)
    if destroy:
        debris = Debris(m_total, centre, v_com, r_sum, composition, 0.0)
        return finish("disruption", [], [a.key, b.key], debris)
    if (not grazing and x < tuning.COLLISION_EROSION_X) or v_imp < tuning.COLLISION_MERGE_SPEED_RATIO * v_esc:
        survivor, gone = merged_into_larger(m_total, merged_radius, composition=composition)
        return finish("merger", [survivor], [gone])
    if grazing and v_imp <= tuning.COLLISION_HIT_AND_RUN_SPEED_RATIO * v_esc and x < 1.0:
        survivor, gone = merged_into_larger(m_total, merged_radius, composition=composition)
        return finish("graze_and_merge", [survivor], [gone])
    if x >= 1.0:
        remnant = largest_remnant_fraction(x) * m_total
        debris = Debris(m_total, centre, v_com, r_sum, composition, remnant)
        return finish("disruption", [], [a.key, b.key], debris)
    if grazing:
        eroded = x / 2.0 * 2.0 * m_proj
        share = eroded / m_total
        kept = []
        taken = [0.0, 0.0]
        for index, body in enumerate((a, b)):
            taken[index] = eroded / 2.0
            mass = body.mass_kg - taken[index]
            kept.append(_survivor(body, mass, body.radius_m * (mass / body.mass_kg) ** (1.0 / 3.0), body.position_m,
                                  body.velocity_ms))
        if share >= tuning.COLLISION_FIELD_MIN_SHARE:
            momentum = tuple(taken[0] * u + taken[1] * v for u, v in zip(a.velocity_ms, b.velocity_ms))
            velocity = tuple(c / eroded for c in momentum)
            position = _weighted(taken[0], a.position_m, taken[1], b.position_m)
            return finish("hit_and_run", kept, [], Debris(eroded, position, velocity, r_sum, composition, 0.0))
        return finish("hit_and_run", kept, [], None, lost=eroded)
    remnant = largest_remnant_fraction(x) * m_total
    radius = merged_radius * (remnant / m_total) ** (1.0 / 3.0)
    survivor, gone = merged_into_larger(remnant, radius, composition=composition)
    return finish("erosion", [survivor], [gone], Debris(m_total - remnant, centre, v_com, r_sum, composition, 0.0))


def conservation_error(a, b, outcome):
    """
    The relative error of an outcome's mass and momentum: `(mass, momentum)`, the survivors, the debris and the gas a
    giant lost (carried at the pair's centre-of-mass velocity, as the hit-and-run strip is at the bodies' own) against the
    two parents. A hit-and-run's small unfielded strip is counted in `lost_mass_kg` at the mean of its parents' velocity.
    """
    parents = (a, b)
    mass_in = sum(p.mass_kg for p in parents)
    momentum_in = tuple(sum(p.mass_kg * p.velocity_ms[i] for p in parents) for i in range(3))
    mass_out = sum(s.mass_kg for s in outcome.survivors) + outcome.lost_mass_kg
    momentum_out = [sum(s.mass_kg * s.velocity_ms[i] for s in outcome.survivors) for i in range(3)]
    if outcome.debris is not None:
        mass_out += outcome.debris.mass_kg
        for i in range(3):
            momentum_out[i] += outcome.debris.mass_kg * outcome.debris.velocity_ms[i]
    if outcome.lost_mass_kg:
        carrier = _velocity(a, b) if outcome.code != "hit_and_run" else _weighted(1.0, a.velocity_ms, 1.0, b.velocity_ms)
        for i in range(3):
            momentum_out[i] += outcome.lost_mass_kg * carrier[i]
    scale = sum(p.mass_kg * _norm(p.velocity_ms) for p in parents) or 1.0
    return (abs(mass_out - mass_in) / mass_in,
            _norm(tuple(o - i for o, i in zip(momentum_out, momentum_in))) / scale)
