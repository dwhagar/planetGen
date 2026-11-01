# planetgen/generation/phenomena/system_types.py

"""
The system-type registry (GEN.268): `SYSTEM_TYPE_RATES`, one entry per row of
the catalogue of multi-star and star-like-object system types
(docs/design/systems/multi-star-systems.md, sections 13.1 to 13.6; 57 types), plus
the tuning knobs of section 13.7. It sits beside the phenomena registry
(`registry.py`) and holds data only.

Every entry carries the reference class its fraction applies to, the
default, the low and high of its range, the evidence tag, the populations
that may host it (the sector profile's `mix`), the generator that draws it,
the number of the GEN item that builds that generator (`None` when it exists
already; a number, not a code string, so no TODO code is shown anywhere), whether it is drawn at scatter time
(S), at fill (F) or both, and whether it is theoretical (no observation).
The member-tree generator and every later item of the multi-star train read
their rates from here and hold none of their own; changing a default is an
edit to this file, with no change to the generator code. The audit counts of
the catalogue (section 15) read the same entries.

A type whose rate is not a number of its own (it comes out of the Roche
classifier, or it is another entry's share) has `default=None` and says where
the number comes from in `derived`.
"""

from dataclasses import dataclass

GALAXY_SYSTEMS = 2.953e11
"""float: Systems in the default galaxy (the note's [C] figure); a `count`
entry is that many objects in it, so its fraction of systems is
`count / GALAXY_SYSTEMS`."""

GENERATOR_MODULE = "planetgen.generation.multistar"
"""str: The module the generators are written in; `SystemType.generator` is a
function name in it. The item named in `SystemType.owner` adds the function
when it builds the type. A generator that already exists is a full
`module:attr` path with owner `None`. The name is checked as an identifier
until its owner has landed; the item that adds the function lists its own
number in `test_system_type_rates._LANDED_OWNERS`, and the test then fails for
a name the module lacks. The module exists since the member-tree generator."""

POPULATIONS = ("young", "intermediate", "old", "bulge", "halo", "cluster", "satellite", "nuclear_cluster",
               "nuclear_disc", "satellite_young", "satellite_intermediate", *(f"open_cluster_{k:02d}" for k in range(12)))
"""tuple: The stellar populations of the sector profile's `mix` a type may be
restricted to; `ANY` means all of them."""

AUDITS = (10, 11, 12, 13, 14)
"""tuple: The audit numbers of the note's section 15: 10 type counts, 11
region gates, 12 the N tail, 13 interaction consistency, 14 group decay."""

MODES = ("S", "F", "S+F")
"""tuple: The valid `SystemType.mode` values: scatter, fill, both."""

UNITS = ("fraction", "count")
"""tuple: `fraction` of the reference class, or `count` in the default galaxy."""

ANY = POPULATIONS


@dataclass(frozen=True, eq=False)
class SystemType:
    """One entry of `SYSTEM_TYPE_RATES` (see the module docstring)."""
    section: str
    row: str
    unit: str
    reference: str
    default: float | None
    low: float | None
    high: float | None
    tag: str
    populations: tuple
    generator: str
    owner: int | None
    mode: str
    theoretical: bool
    audit: tuple
    by_class: dict | None
    derived: str | None
    note: str


def _t(section, row, reference, default, tag, populations, generator, owner, mode, audit=(10,), *,
       low=None, high=None, unit="fraction", theoretical=False, by_class=None, derived=None, note=""):
    return SystemType(section, row, unit, reference, default, low, high, tag, tuple(populations), generator,
                      owner, mode, theoretical, tuple(audit), by_class, derived, note)


def _count_fraction(count):
    return count / GALAXY_SYSTEMS


_TREE = "draw_member_tree"

SYSTEM_TYPE_RATES = {
    # --- 13.1 Architecture of stellar systems ---
    "close-binary": _t(
        "13.1", "Spectroscopic / close detached binary, a under 0.26 AU", "binaries", 0.066, "C", ANY,
        _TREE, 261, "S+F", (10, 11),
        note="Log-normal period law; the close fraction rises to 40% to 53% at [Fe/H] -1 to -3 (halo, thick disk). "
             "Eclipsing is derived from the inclination."),
    "intermediate-binary": _t(
        "13.1", "Intermediate binary, 0.26 to 50 AU", "binaries", 0.433, "C", ANY, _TREE, 261, "S+F"),
    "wide-binary": _t(
        "13.1", "Wide binary, 50 to 10,000 AU", "binaries", 0.436, "C", ANY, _TREE, 261, "S+F", (10, 11),
        note="Thinned in the bulge, clusters and halo."),
    "ultra-wide-pair": _t(
        "13.1", "Ultra-wide pair, 10,000 to 20,000 AU (0.05 to 0.1 pc)", "binaries", 0.065, "C",
        ("young", "intermediate", "old"), _TREE, 261, "S+F", (10, 11),
        note="Flag `fragile`; the cap is 20,000 AU for the root, a directive raises it."),
    "triple-2plus1": _t(
        "13.1", "Triple 2+1 (a close pair with a wide third)", "FGK systems", 0.08, "S", ANY, _TREE, 261, "S+F",
        by_class={"B": 0.22, "O": 0.30}, note="50% and 70% for B and O in a recent summary."),
    "quadruple-3plus1": _t(
        "13.1", "Quadruple 3+1 (a triple with a wide fourth)", "quadruples", 0.60, "I", ANY, _TREE, 261, "S+F"),
    "quadruple-2plus2": _t(
        "13.1", "Quadruple 2+2 (two close pairs bound wide)", "quadruples", 0.40, "I", ANY, _TREE, 261, "S+F",
        by_class={"F and G systems": 0.04}),
    "quintuple-to-septuple": _t(
        "13.1", "Quintuple to septuple", "systems", 2.4e-3, "S,C", ANY, _TREE, 261, "S+F", (10, 12),
        by_class={"N5": 2.4e-3, "N6": 2.9e-4, "N7": 3.0e-5},
        note="Richest in young massive regions. The multi-star prevalence knob scales the N >= 5 columns."),
    "eight-to-twelve-stars": _t(
        "13.1", "Eight to twelve stars", "systems", _count_fraction(1.1e6), "C", ("young", "intermediate"),
        _TREE, 261, "S+F", (10, 12),
        by_class={"N8": _count_fraction(1.1e6), "N9": _count_fraction(1.1e5), "N10": _count_fraction(1.1e4),
                  "N11": _count_fraction(1.1e3), "N12": _count_fraction(1.1e2)},
        note="Ratio 0.12 for N = 8 and 0.1 per step after; cap 12 stars. None confirmed; ten candidates "
             "with multiplicity up to 12."),
    "young-group": _t(
        "13.1", "Young non-hierarchical group (trapezium)", "massive systems younger than 0.3 Myr", 0.05, "I",
        ("young",), "draw_young_group", 273, "S+F", (10, 11, 14),
        note="O and early B systems; stored as a `group` root node; drawn at fill in star-forming cores and, "
             "rarely, at scatter. Never drawn older than its decay time."),
    "decayed-trapezium-remnant": _t(
        "13.1", "Decayed trapezium remnant", "groups older than their decay age", 1.0, "S", ("young", "cluster"),
        "draw_young_group", 273, "S+F", (10, 14),
        note="A group older than its decay time is drawn as its result: an ejected lightest star (a runaway row) "
             "plus a tighter hierarchical remnant."),
    "multiple-with-disk": _t(
        "13.1", "Multiple with a circumbinary or circumstellar disk", "stars younger than 1 Myr", 0.8, "S",
        ("young",), "draw_disk_descriptor", 242, "F", (10, 11), low=0.1, high=0.8,
        note="Disk fraction falls from about 0.8 at 1 Myr to 0.1 at 5 to 6 Myr; a descriptor on the node, no bodies."),

    # --- 13.2 Interacting binaries and merger products ---
    "contact-binary": _t(
        "13.2", "Contact binary (W UMa, over-contact)", "FGK dwarf systems", 0.002, "S", ("intermediate", "old"),
        "classify_interaction", 269, "S+F", (10, 11, 13), low=0.002, high=0.0077,
        by_class={"tertiary": 0.60, "massive": 0.01},
        note="1 in 500 (older estimate 1 in 130; up to 1 in 30 among blue stragglers). Emergent: both Roche lobes "
             "fill at the age. Default tertiary share 60%; massive over-contact 1% of O and B binaries [I, low]."),
    "algol": _t(
        "13.2", "Semi-detached (Algol)", "disk stars", 1.5e-3, "S", ("young", "intermediate"),
        "classify_interaction", 269, "S+F", (10, 11, 13), low=1e-3, high=2e-3,
        note="A lower limit. Emergent: the donor fills its lobe; `mass_reversed` marks a donor that was heavier "
             "at birth."),
    "blue-straggler-field": _t(
        "13.2", "Blue straggler and merger product, field", "old and intermediate stars", 5e-3, "I",
        ("intermediate", "old", "halo", "cluster"), "draw_merger_product", 270, "F", (10, 11),
        note="Triple mergers explain 20% to 25% of blue stragglers; a single star flagged `merger-product` with "
             "a rejuvenated age. Cluster blue stragglers are drawn by the cluster rules."),
    "massive-merger-product": _t(
        "13.2", "Massive-star merger product and mass gainer", "massive main-sequence stars", 0.08, "R",
        ("young",), "draw_merger_product", 270, "S", (10, 11),
        by_class={"mass_gainer": 0.20, "exchange_mass_with_companion": 0.70},
        note="About 8% are mergers and about 20% have gained mass; 70% of O stars exchange mass with a companion."),
    "hot-subdwarf-stripped": _t(
        "13.2", "Hot subdwarf and stripped helium star", "post-common-envelope pairs", None, "R",
        ("old", "young"), "draw_stripped_helium", 271, "F", (10, 13),
        derived="a conditional state of a post-common-envelope pair (knob post_ce_fraction); companion mass 0.3 to "
                "0.5 solar masses, kind `stripped-helium`"),
    "be-star": _t(
        "13.2", "Be star", "B stars", 0.20, "R", ("young",), "draw_be_star", 270, "S+F", (10, 11),
        by_class={"stripped_companion": 0.30},
        note="Flag `be` with a decretion disk; a few have a neutron star companion (Be/X-ray, below)."),
    "wolf-rayet-binary": _t(
        "13.2", "Wolf-Rayet binary", "Wolf-Rayet stars", 0.52, "S", ("young",), "draw_wolf_rayet_companion",
        243, "S+F", (10, 11), by_class={"WN": 0.52, "WC": 0.90}, note="Dsilva et al."),
    "luminous-blue-variable": _t(
        "13.2", "Luminous blue variable and hypergiant", "Galaxy count", 30, "R", ("young",), "draw_lbv_hypergiant",
        246, "S+F", (10, 11), unit="count", note="Tens in the Galaxy; as massive stars (binary 0.9)."),
    "rcb-fk-comae": _t(
        "13.2", "R Coronae Borealis and FK Comae stars", "Galaxy count", 300, "R", ("old",), "draw_merger_product",
        270, "F", (10, 11), unit="count", by_class={"FK Comae": 36},
        note="A few hundred R CrB stars and a few dozen FK Comae stars; single, flagged `merger-product`."),
    "thorne-zytkow": _t(
        "13.2", "Thorne-Zytkow object (a neutron star inside a giant)", "Galaxy count", 100, "S",
        ("old",), "draw_exotic_star", 270, "F", (10, 11), unit="count", low=20, high=200, theoretical=True,
        note="None confirmed; 3.4e-10 of systems. A giant with a neutron-star core, no planets, no companion; "
             "new kind `thorne-zytkow`."),

    # --- 13.3 White-dwarf systems ---
    "white-dwarf-companion": _t(
        "13.3", "White dwarf with a companion", "white dwarfs", 0.30, "S", ("intermediate", "old", "halo", "cluster"),
        _TREE, 261, "F", (10, 11), low=0.25, high=0.32, by_class={"wide": 0.22},
        note="25% to 26% observed, 32% in the 20 pc sample."),
    "close-double-white-dwarf": _t(
        "13.3", "Close double white dwarf", "white dwarfs", 0.095, "S", ("old",), "draw_white_dwarf_binary",
        271, "F", (10, 11), note="Within 4 AU (Maoz et al.); axis log-uniform 0.001 to 4 AU; Peters rule."),
    "double-white-dwarf-merger": _t(
        "13.3", "Double white dwarf merger product", "white dwarfs", 0.10, "S", ("old",), "draw_white_dwarf_merger",
        271, "F", (10, 11), by_class={"wide DWD heavier star": 0.26},
        note="The high-mass bump, 0.8 to 1.1 solar masses; 26% to 54% of wide double white dwarfs have a "
             "merger-product heavier star."),
    "white-dwarf-post-ce": _t(
        "13.3", "White dwarf plus main-sequence star, post-common-envelope pair", "white-dwarf companions", 0.10,
        "I", ("old", "halo", "cluster"), "draw_white_dwarf_binary", 271, "F", (10, 13),
        derived="the knob post_ce_fraction", note="Axis 0.005 to 2 AU, log-uniform; 22% are wide and never interacted."),
    "cataclysmic-variable": _t(
        "13.3", "Cataclysmic variable (all)", "white dwarfs", 6.5e-4, "S,C", ("intermediate", "old"),
        "draw_white_dwarf_binary", 271, "F", (10, 11, 13), by_class={"magnetic": 0.36},
        note="4.8e-6 per pc^3 from a 150 pc Gaia sample, 77% complete; 3.4e-5 per system, 1.0e7 in the Galaxy. "
             "Orbital period 1.3 to 10 h."),
    "am-cvn": _t(
        "13.3", "AM Canum Venaticorum binary", "white dwarfs", 7.5e-5, "S,C", ("old",), "draw_white_dwarf_binary",
        271, "F", (10, 13),
        note="5.5e-7 per pc^3 (+-3.7); 3.9e-6 per system, 1.2e6 in the Galaxy; period 5 to 65 min."),
    "symbiotic-barium": _t(
        "13.3", "Nova, symbiotic, Type Ia progenitor", "giants", 5e-6, "I",
        ("old", "bulge", "nuclear_cluster", "nuclear_disc"),
        "draw_giant_companion", 271, "S+F", (10, 11), by_class={"barium and CH (G and K giants)": 0.01},
        note="Symbiotic 276 confirmed, predicted 1.2e3 to 4e5 (low confidence). A Type Ia explosion is an event "
             "and out of scope."),
    "white-dwarf-brown-dwarf": _t(
        "13.3", "White dwarf plus brown dwarf", "white dwarfs", None, "R", ("old",), "draw_substellar_companion",
        275, "F", (10,), derived="part of the brown-dwarf companion rates (13.5)"),
    "white-dwarf-pollution": _t(
        "13.3", "White dwarf with a planet or debris", "cool white dwarfs", 0.375, "S", ("old",),
        "draw_white_dwarf_pollution", 253, "F", (10,), low=0.25, high=0.50,
        note="25% to 50% of cool white dwarfs are polluted; the default is the middle of the range."),

    # --- 13.4 Compact-object systems ---
    "high-mass-xray-binary": _t(
        "13.4", "High-mass X-ray binary (supergiant or Be)", "compact companions of living stars of 14 solar "
        "masses or more", 0.40, "S", ("young",), "draw_accreting_compact", 272, "S+F", (10, 11, 13),
        by_class={"be_xray_share": 0.60, "compact_companion": 0.003},
        note="The compact companion is 0.3% of such stars (0.05% to 2%); 40% of them accrete."),
    "low-mass-xray-binary": _t(
        "13.4", "Low-mass X-ray binary", "neutron stars with a low-mass star", 5e-5, "I",
        ("bulge", "old", "nuclear_cluster", "nuclear_disc"),
        "draw_accreting_compact", 272, "F", (10, 11, 13), by_class={"black hole": 1e-3, "accreting": 0.01},
        note="Companion of a low-mass star: 5e-5 of neutron stars and 1e-3 of black holes, then an accreting "
             "fraction of 1%."),
    "microquasar": _t(
        "13.4", "Microquasar and ultraluminous X-ray source", "Galaxy count", 50, "R", ("young",),
        "draw_accreting_compact", 272, "F", (10, 13), unit="count", low=24, high=1000,
        note="Several dozen; a flag on an `accreting` pair; jets along the spin axis."),
    "millisecond-pulsar-companion": _t(
        "13.4", "Pulsar with a white dwarf (millisecond pulsar)", "millisecond pulsars", 0.78, "S",
        ("old", "bulge", "cluster", "nuclear_cluster", "nuclear_disc"), "draw_accreting_compact", 272, "F", (10, 11),
        by_class={"white-dwarf": 0.8474, "black-widow": 0.08, "redback": 0.04, "neutron-star": 0.0026, "black-hole": 3e-5,
                  "other": 0.03},
        note="Companion mass gap 0.07 to 0.1 solar masses is empty; spiders carry an ablation flag. A neutron star "
             "companion is 0.2 percent of recycled pulsars (the double neutron star rate allows about 0.04 percent of "
             "neutron stars; 5 percent would be ten times too many; stand-in row AC10)."),
    "pulsar-massive-companion": _t(
        "13.4", "Pulsar with a main-sequence or Be star", "young neutron stars", 0.02, "S", ("young",),
        _TREE, 261, "F", (10, 11), low=0.01, high=0.03,
        note="The 1% to 3% young neutron star companion share (Be or main-sequence star)."),
    "double-neutron-star": _t(
        "13.4", "Double neutron star", "neutron stars", 4e-4, "C", ANY, _TREE, 261, "S+F", (10,),
        low=1.4e-4, high=9.5e-4, note="Promoted to the scatter."),
    "neutron-star-black-hole": _t(
        "13.4", "Neutron star with a black hole", "black holes", 5e-4, "C", ANY, _TREE, 261, "S+F", (10,),
        low=7.5e-5, high=7e-4, note="Gravitational-wave rate 9 to 86 per Gpc^3 per year."),
    "black-hole-black-hole": _t(
        "13.4", "Black hole with a black hole", "black holes", 2e-3, "C", ("young", "intermediate", "old", "halo",
        "cluster"), _TREE, 261, "S+F", (10, 11), low=1.7e-4, high=2.4e-2,
        note="More at low metallicity (halo, clusters, the Clouds)."),
    "dormant-black-hole-star": _t(
        "13.4", "Dormant black hole with a star", "black holes", 0.010, "S,C", ("intermediate", "old", "halo"),
        _TREE, 261, "S+F", (10,), low=1e-3, high=2.2e-2, note="Gaia BH1, BH2 and BH3."),
    "dormant-neutron-star-star": _t(
        "13.4", "Dormant neutron star with a star", "young neutron stars", None, "S", ("intermediate", "old", "halo"),
        _TREE, 261, "S+F", (10,),
        derived="part of the young-neutron-star companion share (pulsar-massive-companion)",
        note="Gaia NS1 and 21 astrometric candidates near 1.4 solar masses."),
    "compact-triple": _t(
        "13.4", "Compact triples and quadruples", "compact binaries", 0.30, "I", ANY, _TREE, 261, "S+F", (10,),
        by_class={"quadruple with two compact pairs": 1e-6}, theoretical=True,
        note="At most three compact members; the quadruple with two compact pairs is theoretical (under 1e-6)."),
    "cluster-compact-binary": _t(
        "13.4", "Black hole or neutron star in a cluster", "black holes in clusters", 0.10, "S", ("cluster",),
        "draw_cluster_compact_binary", 274, "F", (10, 11),
        by_class={"millisecond-pulsar-companion": 0.80},
        note="Cluster rates replace the field fractions; dynamical black-hole pairs come from the cluster density."),
    "intermediate-mass-black-hole-star": _t(
        "13.4", "Intermediate-mass black hole with a star", "intermediate-mass black holes", 0.01, "R",
        ("cluster", "halo"), _TREE, 261, "S", (10, 11), note="A captured star; no other companion."),
    "magnetar-binary": _t(
        "13.4", "Magnetar in a binary", "magnetars", 0.0, "R", ("young",), _TREE, 261, "F", (10, 11),
        theoretical=True, note="None confirmed; 0 by default."),
    "supermassive-black-hole-binary": _t(
        "13.4", "Supermassive black-hole binary", "galaxies", 0.0, "S", ("bulge",), "nucleus_text_flag", 155,
        "F", (10,), theoretical=True,
        note="Extragalactic; not a star system: the nucleus text flag only."),

    # --- 13.5 Substellar companions, planetary-mass companions and disks ---
    "brown-dwarf-close": _t(
        "13.5", "Brown dwarf close companion", "Sun-like stars", 0.006, "S", ANY, "draw_substellar_companion",
        275, "S+F", (10,), low=0.0, high=0.006, note="The brown-dwarf desert: under 0.6% inside 3 to 10 AU."),
    "brown-dwarf-wide": _t(
        "13.5", "Brown dwarf wide companion", "Sun-like stars", 0.032, "S", ANY, "draw_substellar_companion",
        275, "S+F", (10,), low=0.032, high=0.068, note="At 28 to 1,590 AU; 6.8% in one thesis."),
    "brown-dwarf-pair": _t(
        "13.5", "Brown dwarf pair or brown dwarf plus low-mass star", "brown dwarfs", 0.15, "R", ANY,
        "draw_substellar_companion", 275, "F", (10,)),
    "planetary-mass-companion": _t(
        "13.5", "Planetary-mass companion (a few Jupiter masses, wide)", "young stars", 0.03, "R", ("young",),
        "draw_substellar_companion", 275, "F", (10, 11), theoretical=True,
        note="Directly imaged companions at tens to hundreds of AU, a few per cent; a `brown-dwarf` member "
             "below 13 Jupiter masses, flagged `planetary-mass`."),
    "circumbinary-planet": _t(
        "13.5", "Circumbinary planet", "close eclipsing binaries", 0.10, "S", ANY, "draw_circumbinary",
        275, "F", (10,), low=0.10, high=0.47,
        note="At least 10% host a planet above 6 Earth radii at periods under 300 d if coplanar (47% if "
             "isotropic); Holman-Wiegert P-type zone from the node."),
    "circumbinary-disk": _t(
        "13.5", "Circumbinary disk", "young and evolved binaries", None, "S", ("young", "old"),
        "draw_circumbinary", 275, "F", (10, 11),
        derived="the disk fraction of multiple-with-disk (young) and the planetary-nebula binary share (evolved)",
        note="Disks but no planets; a descriptor on the node."),
    "s-type-planet": _t(
        "13.5", "Planet in a wide pair (S-type)", "wide pairs", None, "S", ANY,
        "planetgen.generation.system:StarSystem", None, "F", (10,), derived="the planet model at fill, per star"),
    "rogue-planet-pair": _t(
        "13.5", "Rogue planet pair (a rogue and its moon, or two rogues)", "rogue planets", 0.2, "S", ("young",),
        "planetgen.generation.phenomena.rogue:RoguePlanet", None, "F", (10,),
        note="The existing moon rule (`hosts.PLANET_SATELLITES`)."),

    # --- 13.6 Ejection and runaways ---
    "runaway-star": _t(
        "13.6", "Runaway star", "O stars", 0.20, "S", ("young", "cluster"), "draw_runaway_origin", 274,
        "S+F", (10, 11), low=0.10, high=0.30,
        by_class={"B": 0.075, "supernova channel": 2 / 3, "dynamical channel": 1 / 3},
        note="10% to 30% of O stars and 5% to 10% of B stars [R1558]; a bow-shock survey finds 23%. The "
             "supernova channel leaves a compact companion behind at the rate the survival rule gives."),
    "hypervelocity-star": _t(
        "13.6", "Hypervelocity star", "Galaxy count", None, "S", ("halo",),
        "planetgen.generation.run_sector:_add_hypervelocity_star", None, "S", (10,), unit="count", low=1e3, high=1e4,
        derived="tuning.PHENOMENON_DENSITY_PC3['hypervelocity-star'] (the existing scatter rate stays there)",
        note="Single: the Hills mechanism ejects one star of a disrupted binary."),
    "unbound-by-collapse": _t(
        "13.6", "Companion unbound by a collapse", "non-interacting massive binaries", 0.83, "C", ("young",),
        "draw_runaway_origin", 274, "F", (10,),
        note="Not stored; the survivor is an ordinary row, with a flag if it is a runaway."),
}
"""dict: Type key -> `SystemType`, one per row of the catalogue (sections 13.1
to 13.6). `row` is the catalogue row's first cell, exactly as the note spells
it; the completeness test reads the note and fails when a row has no entry."""


TUNING_KNOBS = {
    "tight_tail_weight": 6.0,
    "post_ce_fraction": 0.10,
    "accreting_fraction": {"low-mass-xray-binary": 0.01, "high-mass-xray-binary": 0.40},
    "merger_product_fraction": {"massive": 0.08, "field": 0.005, "white-dwarf": 0.10},
}
"""dict: The knobs of section 13.7. `tight_tail_weight` multiplies the density of
pairs with a period under a day for primaries under the massive limit, until
the contact share of FGK systems is inside 0.2% to 0.77% and the Algol share
of a constant-star-formation disk approaches 0.1% (GEN.270 set it to 6.0, the
plain draw being 1.0); `post_ce_fraction` is the post-common-envelope
share of white-dwarf companions (the Maoz 9.5% for double white dwarfs and a
matching value for white dwarf plus star); `accreting_fraction` is the share of
compact companions that accrete, by class (1% to 40%); `merger_product_fraction`
is the merger-product share of massive (8%), field (0.5%) and white-dwarf (10%)
stars. All are aligned against the audit counts of the note's section 15."""


def fraction_of_systems(key):
    """The type's default as a fraction of all systems, for a `count` entry
    (its count over `GALAXY_SYSTEMS`); a `fraction` entry's default is
    already a fraction of its own reference class and is returned as it is."""
    entry = SYSTEM_TYPE_RATES[key]
    if entry.default is None:
        raise ValueError(f"{key} has no default of its own: {entry.derived}")
    return entry.default / GALAXY_SYSTEMS if entry.unit == "count" else entry.default
