# planetgen/physics/comet_coma.py

"""
The coma of an active comet (GEN.455, part A: docs/design/composition/atmosphere-composition.md finding 17 and 27, part 04 of its parts and G1-comets.md)
=============================================================================================================================================

What a comet shows when its ice sublimates: how much water it loses a second, what else leaves with it, and how far each gas goes before the star's far
ultraviolet breaks it. Nothing here reads a Solar System comet: the nucleus and its ices are the comet's own (`draw_comet`, GEN.442), the light is its host's.

1. **The water.** The sublimation is limited by the light: the nucleus takes the light of its projected circle, `pi R^2 (1 - A) S / r^2` at the lit distance `r` (the
   distance that is lit like `r` AU round the Sun, `comet_activity.solar_equivalent_au`), and the ice that lies bare on a share `f` of its surface spends all of it
   on the molecules it frees, each at the latent heat of sublimation. So

       Q(H2O) at 1 AU lit = f * pi R^2 * (1 - A) * S / (h_sub / N_A)       molecules a second

   The share `f` (`ACTIVE_AREA_FRACTION_RANGE`, CZ1) and the nucleus Bond albedo `A` (`NUCLEUS_ALBEDO_RANGE`, CZ3) are drawn per comet from their observed ranges;
   `h_sub` is the shared materials table's heat of vaporisation plus that of fusion of water ice (derived, not a stand-in). Farther in, a comet makes more as the
   light grows as `1 / r^2` (`ACTIVITY_INDEX`, physics: the energy limit), and the water turns off at `comet_activity`'s limit, 3 AU of lit distance
   (`tuning.COMET_ACTIVITY_PERIHELION_THRESHOLD_AU`).
2. **The species.** Each other gas leaves at the nucleus's own ratio to water (`makeup.detail["volatiles"]`, percent by number): `Q(X) = Q(H2O) * ratio`. A gas the
   nucleus does not hold is not in the coma. The gases whose ice evaporates colder than water's (`comet_ice_budget.ICE_LINES_K`: CO2 and CO) are not switched off at
   the water's limit, but keep the same index farther out.
3. **The extent.** A gas flows out at `v(r) = v1 / sqrt(r)` (`OUTFLOW_KMS_1AU_RANGE`, CZ2: the expansion of a gas that the light heats, a speed of about a kilometre a second
   at 1 AU) and lives `tau(r) = tau1 * r^2` against the host's far ultraviolet (`photolysis.lifetime_s`, GEN.496, the host's temperature and the quiet Sun's lifetime at
   1 AU lit), so its parent scale length is `L(r) = v tau = L1 * r^1.5`, with `L1` the length at 1 AU lit. A gas with no photolysis row has no length (stored as none, never
   guessed).

A comet stores its numbers at 1 AU lit distance, so the rate and the length at any distance, and for the comet's own orbit as it moves, are computed when asked
(`production_at`, `scale_length_at`). Test targets, never inputs: the observed water production of comets, 5e26 to 1.4e31 molecules a second near 1 AU
(the SWAN sample, finding 17 of the atmosphere composition note), and the nucleus albedos of 0.02 to 0.10. [Stand-ins 6.327 CZ1 to CZ3.]
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import comet_activity, comet_ice_budget, constants, materials, photolysis

ACTIVE_AREA_FRACTION_RANGE = (0.003, 1.0)
"""tuple: The share of a nucleus's surface that lies bare in water ice and sublimates, drawn log-uniform: the widest range, between a nucleus whose few vents cover a third of a
percent and a body that is all ice. [Stand-in 6.327 CZ1, tier T7, tweak L2.]"""

OUTFLOW_KMS_1AU_RANGE = (0.4, 1.0)
"""tuple: The speed, km/s, at which a coma's gas leaves the nucleus at 1 AU lit distance, drawn log-uniform; it falls as the inverse square root of the lit distance. [Stand-in 6.327 CZ2, tier T6, tweak L2.]"""

NUCLEUS_ALBEDO_RANGE = (0.02, 0.10)
"""tuple: The Bond albedo of a comet nucleus, drawn uniform: a dark body, the darkest of the small bodies. [Stand-in 6.327 CZ3, tier T6, tweak L1.]"""

ACTIVITY_INDEX = 2.0
"""float: The power of the lit distance the production falls with, `Q ~ r^-2`: the light the nucleus takes (physics, the energy limit; not a stand-in)."""

SPECIES_FORMULAS = {"H2CO": "CH2O"}
"""dict: A volatile's name in the comet population -> its formula in the photolysis table where they differ."""

WATER = "H2O"
"""str: The species the coma's others are ratios to."""


FIELDS = ("water_production_1au_per_s", "coma_active_area_fraction", "coma_outflow_1au_kms", "nucleus_albedo")
"""tuple: The attributes of a comet that hold its coma's numbers (and the columns of its row): the water, the share of the surface that sublimates, the outflow speed and the albedo, all at 1 AU lit."""

SPECIES_FIELD = "coma_species"
"""str: The attribute that holds the species rows, `(species, production_1au_per_s, scale_length_1au_km or None)`, water first (the rows of `comet_coma`); empty for a comet with no coma."""

UNITS = {"water_production_1au_per_s": "1/s", "coma_active_area_fraction": "ratio", "coma_outflow_1au_kms": "km_s", "nucleus_albedo": "ratio"}
"""dict: The unit each of `FIELDS` is stored in."""


@dataclass(frozen=True)
class Coma:
    """
    A comet's coma at 1 AU lit distance.

    Attributes:
        water_production_1au_s (float): Molecules of water a second at 1 AU lit distance.
        active_area_fraction (float): The share of the surface that sublimates (CZ1).
        outflow_kms_1au (float): The gas speed at 1 AU lit (CZ2).
        nucleus_albedo (float): The Bond albedo (CZ3).
        species (tuple): `(species, production_1au_s, scale_length_1au_km or None)` for water and each volatile the nucleus holds, water first.
    """

    water_production_1au_s: float
    active_area_fraction: float
    outflow_kms_1au: float
    nucleus_albedo: float
    species: tuple


def sublimation_j_per_molecule():
    """The energy a molecule of water ice takes to leave, J: the shared table's heat of vaporisation plus that of fusion over Avogadro's number."""
    props = materials.material("water ice").props
    kilojoule_per_mol = props["h_vap_kj_mol"].value + props["h_fus_kj_mol"].value
    return kilojoule_per_mol * 1000.0 / constants.AVOGADRO


def water_production_1au(nucleus_diameter_km, active_area_fraction, nucleus_albedo):
    """The molecules of water a second a nucleus loses at 1 AU lit distance: `f pi R^2 (1 - A) S / (h_sub / N_A)`, with `S` the starlight at 1 AU round the Sun."""
    if not nucleus_diameter_km > 0.0:
        raise ValueError("a nucleus has a diameter")
    radius_m = nucleus_diameter_km * 500.0
    absorbed_w = active_area_fraction * math.pi * radius_m ** 2 * (1.0 - nucleus_albedo) * constants.SOLAR_CONSTANT_W_M2
    return absorbed_w / sublimation_j_per_molecule()


def scale_length_1au_km(species, outflow_kms_1au, teff_k):
    """The parent scale length, km, of `species` at 1 AU lit distance: the outflow speed times the lifetime against the host's far ultraviolet; `None` for a species the photolysis table lacks."""
    formula = SPECIES_FORMULAS.get(species, species)
    if formula not in photolysis.formulas():
        return None
    lifetime = photolysis.lifetime_s(formula, 1.0, teff_k)
    return None if not math.isfinite(lifetime) else outflow_kms_1au * lifetime


def draw_coma(rng, nucleus_diameter_km, volatiles_percent, teff_k=None):
    """
    The coma of an active comet.

    Args:
        rng: The stream drawn from (`uniform`), a keyed stream of the comet so no other draw moves.
        nucleus_diameter_km (float): The nucleus.
        volatiles_percent (dict): The comet's ratios to water, percent by number (`SmallBodyMakeup.detail["volatiles"]`).
        teff_k (float or None): The host's effective temperature, K; the Sun's for none (a comet with no host).

    Returns:
        Coma: The coma at 1 AU lit distance.
    """
    teff = photolysis.TEMPLATE_TEMPERATURES_K[1] if not teff_k else float(teff_k)
    area = math.exp(rng.uniform(math.log(ACTIVE_AREA_FRACTION_RANGE[0]), math.log(ACTIVE_AREA_FRACTION_RANGE[1])))
    speed = math.exp(rng.uniform(math.log(OUTFLOW_KMS_1AU_RANGE[0]), math.log(OUTFLOW_KMS_1AU_RANGE[1])))
    albedo = rng.uniform(*NUCLEUS_ALBEDO_RANGE)
    water = water_production_1au(nucleus_diameter_km, area, albedo)
    rows = [(WATER, water, scale_length_1au_km(WATER, speed, teff))]
    for species, percent in volatiles_percent.items():
        rows.append((species, water * percent / 100.0, scale_length_1au_km(species, speed, teff)))
    return Coma(water, area, speed, albedo, tuple(rows))


def production_at(rate_1au_s, species, solar_equivalent_au):
    """
    A species' production, molecules a second, at a lit distance: `rate * r^-2`; water and the gases that ride it are none beyond the water's limit
    (`tuning.COMET_ACTIVITY_PERIHELION_THRESHOLD_AU`), the gases whose ice evaporates colder than water's go on.
    """
    if not solar_equivalent_au > 0.0:
        raise ValueError("a lit distance must be positive")
    if solar_equivalent_au > tuning.COMET_ACTIVITY_PERIHELION_THRESHOLD_AU and not _keeps_going(species):
        return 0.0
    return rate_1au_s * solar_equivalent_au ** -ACTIVITY_INDEX


def scale_length_at(length_1au_km, solar_equivalent_au):
    """The parent scale length, km, at a lit distance: `L1 * r^1.5` (the speed falls as `r^-1/2`, the lifetime grows as `r^2`)."""
    if not solar_equivalent_au > 0.0:
        raise ValueError("a lit distance must be positive")
    return length_1au_km * solar_equivalent_au ** 1.5


def _keeps_going(species):
    """True for a gas whose ice evaporates colder than water's (CO2, CO): its coma goes on where the water's has stopped."""
    lines = comet_ice_budget.ICE_LINES_K
    return species in lines and lines[species] < lines[WATER]


def solar_equivalent_au(lit_distance_au, luminosity_w):
    """The distance round the Sun that is lit like `lit_distance_au` round a star of `luminosity_w` (`comet_activity.solar_equivalent_au`)."""
    return comet_activity.solar_equivalent_au(lit_distance_au, luminosity_w)


def attach(body, coma):
    """Sets the coma on a comet or an interstellar comet: `FIELDS` and `SPECIES_FIELD`; `None` (an inactive comet) leaves each of `FIELDS` none and no species."""
    values = (None,) * len(FIELDS) if coma is None else (coma.water_production_1au_s, coma.active_area_fraction, coma.outflow_kms_1au, coma.nucleus_albedo)
    for name, value in zip(FIELDS, values):
        setattr(body, name, value)
    setattr(body, SPECIES_FIELD, () if coma is None else coma.species)


def generate(body, star=None):
    """
    Draws and attaches the coma of an active comet (and the nothing of an inactive one) from its nucleus and ices, on a stream of its own keyed by what makes the
    nucleus (so no other number of the run moves).

    Args:
        body: A `Comet` or an `InterstellarComet` with `makeup`, `nucleus_diameter_km` and `is_active`.
        star: The host, for its temperature; none for a comet that has none.
    """
    if not body.is_active:
        attach(body, None)
        return
    from planetgen.util import draw
    key = f"comet-coma:{float(body.nucleus_diameter_km)!r}:{float(body.mass_kg)!r}:{float(getattr(body, 'perihelion_distance_au', 0.0))!r}:{float(getattr(body, 'velocity_kms', 0.0))!r}"
    attach(body, draw_coma(draw.Stream(key), body.nucleus_diameter_km, body.makeup.detail.get("volatiles") or {}, getattr(star, "temperature", None)))


def of(body):
    """The `Coma` a comet holds (what `attach` set), or `None` for one with none."""
    if getattr(body, FIELDS[0], None) is None:
        return None
    return Coma(body.water_production_1au_per_s, body.coma_active_area_fraction, body.coma_outflow_1au_kms, body.nucleus_albedo, tuple(tuple(row) for row in getattr(body, SPECIES_FIELD, ())))
