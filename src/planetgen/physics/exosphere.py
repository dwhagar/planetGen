# planetgen/physics/exosphere.py

"""
The exosphere of an airless body (GEN.455, part B: part 04 of docs/design/composition/atmosphere-composition-parts, section 5)
==============================================================================================================================

A body with no air at all still has a few atoms above its ground: the ones its surface lets go of, for as long as the light and the wind keep freeing them and
the star's light does not ionise them first. Nothing here reads a Solar System body: the rock, the ices, the star and the distance are the body's own, and the
Moon, Mercury, Io, Europa, Callisto and the comets are test points the numbers are held against, never inputs.

A body gets an exosphere only if it has a solid ground and its stored pressure is zero (a trace of air, however small, is the air's business: `generate` does nothing
for a body with any). The species are not a list by class: a species is in the exosphere exactly where a source supplies it, and the sources read the surface.

1. **The four sources**, each a flux of atoms per square centimetre of the whole sphere, `Phi`:

   * *Thermal release* of the ices the surface exposes: `x_i P_sat,i(T) / sqrt(2 pi m_i k T)` (`saturation.pressure_pa`), `x_i` the molecule share of species `i` in the
     ground ice (the area the ice covers times the share the stored ices hold).
   * *Photon-stimulated desorption* of Na and K held in the regolith: `Q N_S x_i F / 4`, the cross section `Q` (3 and 2 times 1e-20 cm2 for 250 to 400 nm light, the
     first drawn within its 1e-20 spread), `N_S` the 7.5e14 atoms per cm2 of the layer the photons reach, `x_i` the atom share of the element in the rock, and
     `F` the near-ultraviolet photon flux the star puts on the body (the reference flux at the Moon scaled by the star's light and its share in 250 to 400 nm). A quarter,
     because only the lit half is hit and a sphere is four times its disk.
   * *Micrometeoroid impact vaporisation*: the dust rains at `MIV_FLUX_1AU_G_CM2_S` at 1 AU and falls off as `d^-1.5`, each gram vaporises `VAPOR_PER_IMPACTOR` grams
     of the rock, and every element leaves in the rock's own proportions, at the plume's 3,500 K.
   * *Ion sputtering*: the wind's ions at the body's distance, `F_wind Y x_i / 4`, the yield `Y` being the H+ and He+ yields in the wind's 4 percent helium, cut to
     the share of a field's cusps for a body whose own field holds the wind off.

   The three sources of the rock work on the ground the ice does not cover (a body wholly under ice shows none of them).

2. **The loss.** Only the atoms that stay count: a share `f_esc` of each source leaves at once (the Weibull tail of the desorbed, the effusing tail of a Maxwell
   gas, the Thompson tail of the sputtered). The rest fly a ballistic arc of `2 v / g`; during it an atom is ionized by the star's light after `tau` (the Na
   photoionization rate by host temperature, the photolysis table's lifetime for a molecule; none read for the others, which are then lost to the ground only), so its mean time
   up is `tau (1 - exp(-t/tau))`. A returning atom sticks. The column, atoms per cm2, is the flux times that time: `N = Phi (1 - f_esc) tau (1 - exp(-t/tau))`.

3. **The stored record.** A list of `{"species", "column_cm2", "source", "scale_height_km"}`, the largest column first, at most `MAX_SPECIES`, none under `COLUMN_FLOOR_CM2`;
   `None` for a body with air, without a ground or without a star's light. `scale_height_km` is `k T / (m g)` at the source's launch temperature, and a species whose
   scale height is over the body's radius is not recorded (its atoms leave, no exosphere is held), nor is a thermal column over `COLLISIONAL_COLUMN_CM2` (that gas is an atmosphere).

The constants of the sources are read from part 04 section 5 [R2546, R3392, R4783]; the ones that are not are registered stand-ins (section 6.328). Draws come from
a keyed stream of the body, so no other number of a run moves.
"""

import math

from planetgen.physics import constants, gases, life_classes, magnetism, photolysis, radiation_field, saturation
from planetgen.util import draw

FIELD = "exosphere"
"""str: The attribute of a body that holds its exosphere (a list of dicts, or `None`)."""

BOLTZMANN = constants.BOLTZMANN
"""float: The Boltzmann constant, J/K."""

SOURCES = ("thermal release", "photon-stimulated desorption", "micrometeoroid impact vaporisation", "ion sputtering")
"""tuple: The four sources of part 04 section 5.1, the names a record carries."""

MAX_SPECIES = 8
"""int: The most species a record holds (the largest columns)."""

COLUMN_FLOOR_CM2 = 1.0e4
"""float: A column under this, atoms per cm2, is not recorded (a floor of presentation: nothing above a body's ground is measured under it) [Stand-in 6.328 EX8, tier T7, tweak L0.]"""

COLLISIONAL_COLUMN_CM2 = 1.0e15
"""float: The most a thermally released gas can hold above the ground and stay an exosphere, atoms per cm2: where a molecule meets another within a scale height (a gas-kinetic cross section of 1e-15 cm2) the gas is an atmosphere, which the air makes
[Stand-in 6.328 EX9, tier T6, tweak L0.]"""

PSD_CROSS_SECTION_CM2 = {"Na": 3.0e-20, "K": 2.0e-20}
"""dict: The photon-stimulated desorption cross section of the alkali, cm2, for 250 to 400 nm light [R4783]."""

PSD_CROSS_SECTION_SPREAD = 1.0e-20
"""float: The spread of the sodium cross section, cm2, the same fraction (a third) taken for potassium [R4783]: the cross section of a body is drawn within it."""

PSD_SITE_DENSITY_CM2 = 7.5e14
"""float: The atoms per cm2 of the surface layer the photons reach [R4783]."""

PSD_PHOTON_FLUX_REFERENCE_CM2_S = 2.0e15
"""float: The 250 to 400 nm photon flux at the Moon's distance from the Sun, per cm2 and second [R4783]; another star or distance scales it by its light and its share in the band [Stand-in 6.328 EX1, tier T6, tweak L1.]"""

PSD_BAND_NM = (250.0, 400.0)
"""tuple: The band of the desorbing photons, nm."""

PSD_WEIBULL = (575.0, 1.7)
"""tuple: The Weibull scale (m/s) and shape of the speed of the desorbed atoms [R4783]."""

PSD_LAUNCH_K = 1200.0
"""float: The effective temperature, K, of a desorbed atom for the scale height (the Na of Mercury is 1,100 to 1,400 K) [R2683]."""

IMPACT_PLUME_K = 3500.0
"""float: The temperature of the vapor a micrometeoroid makes, K (2,500 to 5,000 K; part 04 uses 3,000 to 4,000 K) [R4783]."""

MIV_FLUX_1AU_RANGE_G_CM2_S = (5.1e-17, 8.4e-17)
"""tuple: The mass of micrometeoroids that hits a square centimetre a second at 1 AU, g, by the Moon's measured range [R4783]; drawn uniform."""

MIV_DISTANCE_EXPONENT = 1.5
"""float: The flux falls as `d^-exponent` from the star [Stand-in 6.328 EX2, tier T7, tweak L1.]: it takes Mercury's range (1.66 to 3.56)e-16 from the Moon's."""

VAPOR_PER_IMPACTOR_RANGE = (1.0, 4.0)
"""tuple: The grams of rock vaporised per gram of micrometeoroid at the speeds of the inner Solar System (14 to 20 km/s), drawn log-uniform [Stand-in 6.328 EX3, tier T6, tweak L1.]"""

SPUTTER_YIELD_H_RANGE = (0.01, 0.1)
"""tuple: The atoms a solar wind proton sputters, drawn log-uniform [R4783]."""

SPUTTER_YIELD_HE_RANGE = (0.1, 1.0)
"""tuple: The atoms a He+ ion sputters, drawn log-uniform [R4783]."""

WIND_HELIUM_SHARE = 0.04
"""float: The helium's share of the wind's ions by number (the solar wind's) [R2546]."""

SPUTTER_BINDING_EV_RANGE = (2.6, 8.4)
"""tuple: The surface binding energy, eV, of the sputtered atom, drawn log-uniform [R4783]."""

WIND_FLUX_1AU_CM2_S = 3.2e8
"""float: The solar wind's ions on a square centimetre a second at 1 AU [R2546]; another star scales it by its wind's pressure at the body [Stand-in 6.328 EX6, tier T7, tweak L1.]"""

SHIELDED_SPUTTER_SHARE = 0.05
"""float: The share of the wind that reaches the surface of a body whose own field stands the wind off (the cusps) [Stand-in 6.328 EX4, tier T7, tweak L1.]"""

RETAINED_SPEED_SHARE = 0.5
"""float: The speed, as a share of the escape speed, a sputtered atom that stays flies its arc at (the ones above leave) [Stand-in 6.328 EX5, tier T7, tweak L1.]"""

NA_IONISATION_RATE_1AU = ((4500.0, 9.5e-7), (5772.0, 5.92e-6), (6800.0, 1.3e-5))
"""tuple: `(host temperature K, photoionization rate of Na at 1 AU, 1/s)` for a K, a G and an F star [R3392]; the rate in between is read linearly in the logarithm of both, outside from the end pair
[Stand-in 6.328 EX7, tier T6, tweak L0.]"""

ICE_FORMULAS = {"water_ice": "H2O", "co2_ice": "CO2", "ch4_ice": "CH4", "nh3_ice": "NH3", "n2_ice": "N2", "so2_ice": "SO2", "co_ice": "CO", "hcn_ice": "HCN"}
"""dict: The ices of the ground the thermal release reads, by the stored ice and the formula of the gas."""

ELEMENTS = {
    "Na": ("Na2O", 2), "K": ("K2O", 2), "Ca": ("CaO", 1), "Mg": ("MgO", 1), "Fe": ("FeO", 1), "SiO": ("SiO2", 1),
}
"""dict: The exospheric species the rock supplies: `{species: (oxide, atoms of the species in a formula unit of the oxide)}`; silicon leaves as SiO."""

OXIDE_ATOMS = {"SiO2": 3, "MgO": 2, "CaO": 2, "Al2O3": 5, "FeO": 2, "Na2O": 3, "K2O": 3}
"""dict: The atoms in a formula unit of each oxide of the rock, for the atom share of an element."""

OXIDE_G_MOL = {"SiO2": 60.084, "MgO": 40.304, "CaO": 56.077, "Al2O3": 101.961, "FeO": 71.844, "Na2O": 61.979, "K2O": 94.196}
"""dict: The molar mass of each oxide of the rock, g/mol."""


def _log_uniform(rng, bounds):
    return math.exp(rng.uniform(math.log(bounds[0]), math.log(bounds[1])))


def na_ionization_rate_1au(teff_k):
    """The photoionization rate of sodium at 1 AU from a star of `teff_k`, 1/s (linear in the logarithm of the rate and of the temperature through the three of `NA_IONISATION_RATE_1AU`)."""
    points = NA_IONISATION_RATE_1AU
    if teff_k <= points[1][0]:
        (t0, r0), (t1, r1) = points[0], points[1]
    else:
        (t0, r0), (t1, r1) = points[1], points[2]
    slope = (math.log(r1) - math.log(r0)) / (math.log(t1) - math.log(t0))
    return math.exp(math.log(r0) + slope * (math.log(teff_k) - math.log(t0)))


def atom_shares(oxides_percent):
    """`{species: share of all atoms}` of the species `ELEMENTS` names, from the rock's oxides in weight percent (the mass of each oxide over its molar mass makes the moles)."""
    total_atoms = 0.0
    moles = {}
    for oxide, percent in oxides_percent.items():
        if oxide in OXIDE_G_MOL and percent > 0.0:
            moles[oxide] = percent / OXIDE_G_MOL[oxide]
            total_atoms += moles[oxide] * OXIDE_ATOMS[oxide]
    if total_atoms <= 0.0:
        return {}
    return {species: moles.get(oxide, 0.0) * count / total_atoms for species, (oxide, count) in ELEMENTS.items() if moles.get(oxide, 0.0) > 0.0}


def mass_per_atom_kg(species):
    """The mass of one atom or molecule of `species`, kg (the gas table's molar mass over Avogadro's number)."""
    return gases.gas(species).molar_mass / constants.AVOGADRO


def escape_speed_m_s(gravity_m_s2, radius_km):
    """`sqrt(2 g R)`, m/s."""
    return math.sqrt(2.0 * gravity_m_s2 * radius_km * 1000.0)


def maxwell_escape_share(mass_kg, temperature_k, escape_m_s):
    """The share of a Maxwell gas's effusing flux that leaves: `(1 + x^2) exp(-x^2)` with `x^2 = m v_esc^2 / (2 k T)`."""
    x2 = mass_kg * escape_m_s ** 2 / (2.0 * BOLTZMANN * temperature_k)
    return (1.0 + x2) * math.exp(-x2) if x2 < 700.0 else 0.0


def weibull_escape_share(escape_m_s):
    """The share of the desorbed that leave: the Weibull tail `exp(-(v_esc / v0)^kappa)` of the speed of PSD."""
    scale, shape = PSD_WEIBULL
    return math.exp(-((escape_m_s / scale) ** shape)) if escape_m_s < 50.0 * scale else 0.0


def sputter_escape_share(mass_kg, binding_ev, escape_m_s):
    """The share of the sputtered that leave: the Thompson law's tail over the escape energy `E0`, `2 Eb / (E0 + Eb) - (Eb / (E0 + Eb))^2`."""
    e0 = 0.5 * mass_kg * escape_m_s ** 2 / constants.ELECTRON_VOLT_J
    ratio = binding_ev / (e0 + binding_ev)
    return 2.0 * ratio - ratio * ratio


def time_up_s(flight_s, lifetime_s):
    """The mean time an atom stays above the ground in one arc of `flight_s`, with a photoionization lifetime `lifetime_s` (`None` for none read): `tau (1 - exp(-t / tau))`, or `t`."""
    if lifetime_s is None or not math.isfinite(lifetime_s) or lifetime_s <= 0.0:
        return flight_s
    return lifetime_s * -math.expm1(-flight_s / lifetime_s)


def flight_s(speed_m_s, gravity_m_s2):
    """The time of a ballistic arc launched straight up at `speed_m_s`, `2 v / g`, s."""
    return 2.0 * speed_m_s / gravity_m_s2


class Inputs:
    """
    What the sources and the losses read of a body, all numbers: a body is turned into this by `inputs_of`, a test builds one by hand.

    Attributes:
        gravity_m_s2, radius_km, surface_temperature_k: The body's.
        oxides_percent (dict): The rock's oxides, weight percent, `{"Na2O": 0.5, ...}`.
        ice_shares (dict): `{formula: share of the ground's molecules}` of the ices the ground exposes (adds to 1 or less).
        ice_cover (float): The share of the ground the ice covers.
        teff_k (float): The host's effective temperature.
        flux_ratio (float): The host's light at the body over the Sun's at 1 AU.
        distance_au (float): The body's distance from the star (its planet's, for a moon).
        wind_ratio (float): The wind's pressure at the body over the Earth's.
        shielded (bool): Whether the body's own field stands the wind off.
    """

    def __init__(self, gravity_m_s2, radius_km, surface_temperature_k, oxides_percent, ice_shares, ice_cover, teff_k, flux_ratio, distance_au, wind_ratio, shielded):
        self.gravity_m_s2 = gravity_m_s2
        self.radius_km = radius_km
        self.surface_temperature_k = surface_temperature_k
        self.oxides_percent = oxides_percent
        self.ice_shares = ice_shares
        self.ice_cover = ice_cover
        self.teff_k = teff_k
        self.flux_ratio = flux_ratio
        self.distance_au = distance_au
        self.wind_ratio = wind_ratio
        self.shielded = shielded


def _entry(species, source, flux_cm2_s, stay, flight, lifetime, temperature_k, gravity):
    """One source's contribution: `{species, source, column_cm2, scale_height_km}`; the thermal release is held under `COLLISIONAL_COLUMN_CM2`."""
    column = flux_cm2_s * stay * time_up_s(flight, lifetime)
    if source == SOURCES[0]:
        column = min(column, COLLISIONAL_COLUMN_CM2)
    height = BOLTZMANN * temperature_k / (mass_per_atom_kg(species) * gravity) / 1000.0
    return {"species": species, "source": source, "column_cm2": column, "scale_height_km": height}


def _ionization_lifetime_s(species, inputs):
    """The photoionization or photolysis lifetime of `species`, s, at the body's light; `None` for a species no rate is read for."""
    if species == "Na":
        return 1.0 / (na_ionization_rate_1au(inputs.teff_k) * inputs.flux_ratio)
    if species in photolysis.formulas():
        return photolysis.lifetime_s(species, inputs.flux_ratio, inputs.teff_k)
    return None


def columns(inputs, rng):
    """
    The columns of every species the four sources give a body, largest first.

    Args:
        inputs (Inputs): The body.
        rng: The stream drawn from (`uniform`), a keyed stream of the body.

    Returns:
        list: Dicts `{"species", "column_cm2", "source", "scale_height_km"}`, the species' sources added and the largest source named, none under `COLUMN_FLOOR_CM2`, at most `MAX_SPECIES`.
    """
    g = inputs.gravity_m_s2
    v_esc = escape_speed_m_s(g, inputs.radius_km)
    shares = atom_shares(inputs.oxides_percent)
    parts = []
    # draws, in a fixed order: the cross section scale, the MIV flux, the vapor per impactor, the yields, the binding energy
    cross_scale = rng.uniform(-1.0, 1.0)
    miv_flux = rng.uniform(*MIV_FLUX_1AU_RANGE_G_CM2_S)
    vapor = _log_uniform(rng, VAPOR_PER_IMPACTOR_RANGE)
    yield_h = _log_uniform(rng, SPUTTER_YIELD_H_RANGE)
    yield_he = _log_uniform(rng, SPUTTER_YIELD_HE_RANGE)
    binding = _log_uniform(rng, SPUTTER_BINDING_EV_RANGE)

    # thermal release of the ground's ices
    temperature = inputs.surface_temperature_k
    if inputs.ice_cover > 0.0 and temperature > 0.0:
        for formula, share in inputs.ice_shares.items():
            if share <= 0.0 or not gases.has(formula):
                continue
            mass = mass_per_atom_kg(formula)
            pressure = saturation.pressure_pa(formula, temperature)
            flux_m2 = inputs.ice_cover * share * pressure / math.sqrt(2.0 * math.pi * mass * BOLTZMANN * temperature)
            speed = min(math.sqrt(8.0 * BOLTZMANN * temperature / (math.pi * mass)), 0.7 * v_esc)
            stay = 1.0 - maxwell_escape_share(mass, temperature, v_esc)
            parts.append(_entry(formula, SOURCES[0], flux_m2 / 1.0e4, stay, flight_s(speed, g), _ionization_lifetime_s(formula, inputs), temperature, g))

    bare = 1.0 - min(max(inputs.ice_cover, 0.0), 1.0)          # the rock sources work on the ground the ice does not cover
    # photon-stimulated desorption of the alkali
    sun_share = life_classes.band_fraction(photolysis.TEMPLATE_TEMPERATURES_K[1], PSD_BAND_NM)
    band_share = life_classes.band_fraction(inputs.teff_k, PSD_BAND_NM)
    photons = PSD_PHOTON_FLUX_REFERENCE_CM2_S * inputs.flux_ratio * (band_share / sun_share if sun_share > 0.0 else 0.0)
    scale = PSD_WEIBULL[0]
    stay_psd = 1.0 - weibull_escape_share(v_esc)
    speed_psd = min(scale * math.gamma(1.0 + 1.0 / PSD_WEIBULL[1]), 0.7 * v_esc)
    for species, cross in PSD_CROSS_SECTION_CM2.items():
        if species not in shares:
            continue
        cross_cm2 = cross * (1.0 + cross_scale * PSD_CROSS_SECTION_SPREAD / PSD_CROSS_SECTION_CM2["Na"])
        flux = bare * cross_cm2 * PSD_SITE_DENSITY_CM2 * shares[species] * photons / 4.0
        parts.append(_entry(species, SOURCES[1], flux, stay_psd, flight_s(speed_psd, g), _ionization_lifetime_s(species, inputs), PSD_LAUNCH_K, g))

    # micrometeoroid impact vaporisation of the refractories
    mass_flux = miv_flux * inputs.distance_au ** -MIV_DISTANCE_EXPONENT * vapor       # g of rock vapor per cm2 per s
    grams_oxides = {oxide: percent for oxide, percent in inputs.oxides_percent.items() if oxide in OXIDE_G_MOL and percent > 0.0}
    total_percent = math.fsum(grams_oxides.values())
    for species, (oxide, count) in ELEMENTS.items():
        if oxide not in grams_oxides or total_percent <= 0.0:
            continue
        atoms_per_g = grams_oxides[oxide] / total_percent / OXIDE_G_MOL[oxide] * count * constants.AVOGADRO
        mass = mass_per_atom_kg(species)
        speed = min(math.sqrt(8.0 * BOLTZMANN * IMPACT_PLUME_K / (math.pi * mass)), RETAINED_SPEED_SHARE * v_esc)
        stay = 1.0 - maxwell_escape_share(mass, IMPACT_PLUME_K, v_esc)
        parts.append(_entry(species, SOURCES[2], bare * mass_flux * atoms_per_g, stay, flight_s(speed, g), _ionization_lifetime_s(species, inputs), IMPACT_PLUME_K, g))

    # ion sputtering of the refractories
    if inputs.wind_ratio > 0.0:
        wind = WIND_FLUX_1AU_CM2_S * inputs.wind_ratio * (SHIELDED_SPUTTER_SHARE if inputs.shielded else 1.0)
        yield_ion = (1.0 - WIND_HELIUM_SHARE) * yield_h + WIND_HELIUM_SHARE * yield_he
        for species, share in shares.items():
            mass = mass_per_atom_kg(species)
            stay = 1.0 - sputter_escape_share(mass, binding, v_esc)
            speed = RETAINED_SPEED_SHARE * v_esc
            parts.append(_entry(species, SOURCES[3], bare * wind * yield_ion * share / 4.0, stay, flight_s(speed, g), _ionization_lifetime_s(species, inputs), IMPACT_PLUME_K, g))

    found = {}
    for entry in parts:
        row = found.setdefault(entry["species"], {"species": entry["species"], "column_cm2": 0.0, "source": entry["source"], "scale_height_km": entry["scale_height_km"], "_best": 0.0})
        row["column_cm2"] += entry["column_cm2"]
        if entry["column_cm2"] > row["_best"]:
            row["_best"], row["source"], row["scale_height_km"] = entry["column_cm2"], entry["source"], entry["scale_height_km"]
    bound = [row for row in found.values() if row["column_cm2"] >= COLUMN_FLOOR_CM2 and row["scale_height_km"] < inputs.radius_km]       # a scale height over the radius is no exosphere: the atoms leave
    ranked = sorted(bound, key=lambda row: -row["column_cm2"])
    return [{key: value for key, value in row.items() if key != "_best"} for row in ranked[:MAX_SPECIES]]


def has_exosphere(body):
    """Whether a planet or moon takes an exosphere: a solid ground and no air at all (a stored pressure of exactly zero)."""
    if getattr(body, "body_type", None) in ("g", "a"):
        return False
    pressure = getattr(body, "atmospheric_pressure", None)
    return pressure is not None and pressure <= 0.0 and getattr(body, "ground_ice_fraction", None) is not None


def _ice_shares(body):
    """`({formula: share of the ground's ice molecules}, cover)` of a body: water by `water_mass_fraction` and the others by the stored ice species, by moles, over the area the ground ice covers."""
    cover = getattr(body, "ground_ice_fraction", None) or 0.0
    masses = {"H2O": getattr(body, "water_mass_fraction", None) or 0.0}
    from planetgen.physics import ice_species
    for key, field in ice_species.FIELD_OF.items():
        if key != ice_species.OTHER:
            masses[key] = getattr(body, field, None) or 0.0
    moles = {formula: mass / (gases.gas(formula).molar_mass * 1000.0) for formula, mass in masses.items() if mass > 0.0 and gases.has(formula)}
    total = math.fsum(moles.values())
    return ({formula: value / total for formula, value in moles.items()} if total > 0.0 else {}), cover


def inputs_of(body):
    """The `Inputs` of a planet or moon, or `None` when it has no star's light or no rock to read."""
    star = getattr(body, "star", None)
    if star is None or not getattr(star, "temperature", None) or not getattr(star, "luminosity", None):
        return None
    distance = getattr(getattr(body, "build", None), "star_distance_au", None)
    if not distance:
        return None
    from planetgen.physics import surface_composition
    oxides = surface_composition._rock_oxides_of(body)
    if oxides is None:
        return None
    lit = radiation_field.lit_distance_au(star, distance)
    flux_ratio = star.luminosity / constants.SOLAR_LUMINOSITY / (lit * lit)       # the star's luminosity is in watts
    shares, cover = _ice_shares(body)
    shielded = (getattr(body, "magnetopause_rp", None) or 0.0) > 1.0
    wind_ratio = magnetism.wind_pressure_at(star, distance) / magnetism.EARTH_WIND_PA
    gravity = (body.gravity or 0.0) * constants.EARTH_GRAVITY
    if not (gravity > 0.0 and body.radius and body.surface_temperature):
        return None
    return Inputs(gravity, body.radius, body.surface_temperature, oxides, shares, cover, star.temperature, flux_ratio, distance, wind_ratio, shielded)


def _seed_of(body):
    return f"exosphere:{body.is_moon}:{body.mass!r}:{body.radius!r}:{body.distance!r}:{getattr(body.star, 'age', None)!r}"


def generate(body):
    """Sets `body.exosphere` (see `FIELD`): the list `columns` makes for a body that takes one (`has_exosphere`), `None` otherwise. The draws come from a keyed stream of the body."""
    setattr(body, FIELD, None)
    if not has_exosphere(body):
        return
    inputs = inputs_of(body)
    if inputs is None:
        return
    found = columns(inputs, draw.Stream(_seed_of(body)))
    setattr(body, FIELD, found or None)
