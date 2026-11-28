# tests/test_family_plausibility.py

"""
TEST.204: a family's drawn galaxies fall inside the ranges the galaxy types design states (`docs/design/galaxy/galaxy-types.md`, GEN.378), over many
seeds and by mass bin, and show the stated scaling slopes. Only the spiral family is built (`shape_draw.draw_galaxy` refuses another), so only its rows
are tested; a family that is built adds its rows to `FAMILIES` and its reference galaxies to `REFERENCES`.

A failure names the parameter, the family and the mass bin and counts how many seeds fell outside, as the census audits do (TEST.196). The ranges are the
design's: a law's median with its scatter, widened to four sigma, and the dials' own published ranges. A draw outside one is a finding for the item that
owns the law; the range is never widened here.

The reference galaxies (the Milky Way, Andromeda) are points the draw must be able to produce, never presets: each is drawn at its own mass and the draws
that land inside its measured range are counted.
"""

import math
import statistics
import types

import numpy as np
import pytest

from planetgen.galaxy import dials, mass_chain, shape_draw, structure as galaxyStructure
from tests.wide_masses import wide_mass

SEEDS = [number.to_bytes(16, "big") for number in range(1, 2501)]
DEFAULT_SEEDS = SEEDS[:120]
"""list: The seeds of the default-galaxy tests: an all-Auto draw redraws until the plan shows no grid warning (GEN.625), which is slow."""
SIGMAS = 4.0
"""float: The width of a drawn parameter's range, in standard deviations of its law's scatter: four sigma holds all but 6 in 100,000 draws of a Gaussian."""

SPIRAL = shape_draw.SPIRAL
MASS_BINS = ((3.0e9, 1.0e10), (1.0e10, 3.0e10), (3.0e10, 1.0e11), (1.0e11, 2.0e11))
"""tuple: The stellar mass bins, solar masses, of the spiral family's range (3e9 to 2e11, GEN.378 section 8)."""
SPIRAL_MASS_LIMITS = (3.0e9, 2.0e11)


def _bin_name(low, high):
    return f"{low:.0e} to {high:.0e} solar masses"


@pytest.fixture(scope="module")
def draws():
    """The seed's own draws over the family's mass range, the mass pinned (the seeds' masses follow the mass function between its limits)."""
    return [shape_draw.draw_galaxy(seed, stellar_mass_solar=wide_mass(seed)) for seed in SEEDS]


@pytest.fixture(scope="module")
def size_draws():
    """The seed's own size draws: a pinned flatness leaves the all-Auto window (GEN.625) out of the size's way."""
    return [shape_draw.draw_galaxy(seed, stellar_mass_solar=wide_mass(seed), flatness=shape_draw.KREGEL_Q) for seed in SEEDS]


# --- The ranges of the spiral family, one row per parameter -------------------------------------------------------------------------------------------

def _size_offset_dex(drawn):
    """How far the disk's scale length is from the size-mass relation's median at its mass, in dex."""
    median_pc = shape_draw.size_median_kpc(drawn.chain.stellar_mass_solar) / shape_draw.EXPONENTIAL_HALF_LIGHT_OVER_SCALE * 1000.0
    return math.log10(drawn.disk_scale_length_pc / median_pc)


def _size_limit_dex(drawn):
    return SIGMAS * shape_draw.size_sigma_ln(drawn.chain.stellar_mass_solar) / math.log(10.0)


def _bulge_offset_dex(drawn):
    """How far the bulge-to-disk ratio is from its stage's median, in dex (the stage's two medians and half widths, `shape_draw.LOG_BD_*`)."""
    bulge_over_disk = drawn.bulge_to_total / (1.0 - drawn.bulge_to_total)
    median, sigma = (shape_draw.LOG_BD_LATE, shape_draw.LOG_BD_LATE_SIGMA) if drawn.stage == "late" else (shape_draw.LOG_BD_EARLY, shape_draw.LOG_BD_EARLY_SIGMA)
    return (math.log10(bulge_over_disk) - median) / sigma


FLATNESS = dials.lever("flatness")
"""Row: The Structure dial's flatness (the thin disk's scale length over its scale height), whose range is the published one."""


def _bar_over_length(drawn):
    return drawn.bar_length_pc / drawn.disk_scale_length_pc


def _holds(value, low, high):
    return low <= value <= high


def _within(offset, limit):
    return -limit <= offset <= limit


# Each row: (parameter, what is read, the test of one drawn galaxy, the fixture that holds the draws). The test gives True inside the range.
ROWS = (
    ("stellar mass", lambda d: d.chain.stellar_mass_solar,
     lambda d: _holds(d.chain.stellar_mass_solar, *SPIRAL_MASS_LIMITS), "draws"),
    ("disk scale length, offset from the size-mass relation", _size_offset_dex,
     lambda d: _within(_size_offset_dex(d), _size_limit_dex(d)), "size_draws"),
    ("bulge-to-disk ratio, offset from its stage's median in sigmas", _bulge_offset_dex,
     lambda d: _within(_bulge_offset_dex(d), SIGMAS), "draws"),
    ("bar length in disk scale lengths", _bar_over_length,
     lambda d: _holds(_bar_over_length(d), *shape_draw.BAR_LENGTH_OVER_H_RANGE), "draws"),
    ("arm count", lambda d: d.arm_count,
     lambda d: d.arm_count in {arms for arms, _share in shape_draw.ARM_COUNT_SHARES}, "draws"),
    ("pitch angle", lambda d: d.pitch_angle_deg,
     lambda d: _holds(d.pitch_angle_deg, *shape_draw.PITCH_RANGE_DEG), "draws"),
    ("bulge-to-total", lambda d: d.bulge_to_total,
     lambda d: _holds(d.bulge_to_total, 0.0, dials.lever("bulge_share").maximum), "draws"),
    # The thin disk scale height has no row yet: the seeded draw lets 0.4 percent of seeds past `DERIVED_LIMITS`; GEN.1902 fixes the draw and adds the row.
    ("bulge scale radius", lambda d: d.bulge_scale_radius_pc,
     lambda d: dials.check_derived("bulge_scale_radius_pc", d.bulge_scale_radius_pc) is None, "draws"),
)


def _outside(rows, fixtures):
    """`{(parameter, mass bin): (outside, total)}` for the rows that have any seed outside."""
    found = {}
    for parameter, _read, inside, fixture in rows:
        for low, high in MASS_BINS:
            bin_draws = [d for d in fixtures[fixture] if low <= d.chain.stellar_mass_solar < high]
            outside = sum(1 for d in bin_draws if not inside(d))
            if outside:
                found[(parameter, _bin_name(low, high))] = (outside, len(bin_draws))
    return found


def test_every_spiral_parameter_falls_inside_its_range_in_every_mass_bin(draws, size_draws):
    found = _outside(ROWS, {"draws": draws, "size_draws": size_draws})
    assert not found, "the spiral family draws outside its range: " + "; ".join(
        f"{parameter} at {name}: {outside} of {total} seeds" for (parameter, name), (outside, total) in sorted(found.items()))


def test_every_mass_bin_has_seeds_so_that_no_row_passes_on_an_empty_bin(draws):
    for low, high in MASS_BINS:
        assert sum(1 for d in draws if low <= d.chain.stellar_mass_solar < high) >= 40, _bin_name(low, high)


def test_a_default_galaxy_is_a_spiral_of_the_family_limits():
    """With nothing pinned the family is the spiral and the mass the Auto window's, inside the family's mass limits."""
    masses = [shape_draw.draw_galaxy(seed).chain.stellar_mass_solar for seed in DEFAULT_SEEDS]
    assert all(SPIRAL_MASS_LIMITS[0] <= mass <= SPIRAL_MASS_LIMITS[1] for mass in masses)
    assert dials.lever("shape").choices == (SPIRAL,)


def test_a_default_galaxy_keeps_its_flatness_inside_the_structure_dials_published_range():
    """GEN.625: the all-Auto window redraws a flatness outside the Structure dial's range, so a default galaxy is one a user could have typed."""
    for seed in DEFAULT_SEEDS:
        drawn = shape_draw.draw_galaxy(seed)
        assert FLATNESS.minimum <= drawn.disk_scale_length_pc / drawn.disk_scale_height_pc <= FLATNESS.maximum


def test_the_flatness_of_a_galaxy_with_a_pinned_mass_follows_kregels_law_inside_the_dials_range(size_draws, draws):
    """The seed's own flatness (the size pinned, so the all-Auto window is out of the way) is Kregel's median with its scatter, and the share of the
    seeds outside the Structure dial's range is the share of a Gaussian of that width: a few in a thousand, not a different law."""
    flat = [shape_draw.draw_galaxy(seed, stellar_mass_solar=wide_mass(seed), compactness=0.0) for seed in SEEDS]
    ratios = [d.disk_scale_length_pc / d.disk_scale_height_pc for d in flat]
    logs = [math.log10(r) for r in ratios]
    assert statistics.mean(logs) == pytest.approx(math.log10(shape_draw.KREGEL_Q), abs=0.02)
    assert statistics.stdev(logs) == pytest.approx(shape_draw.KREGEL_SCATTER_DEX, rel=0.1)
    outside = sum(1 for r in ratios if not FLATNESS.minimum <= r <= FLATNESS.maximum) / len(ratios)
    assert outside < 0.01, f"{100 * outside:.2f} percent of the seeds with a pinned mass draw a flatness outside {FLATNESS.minimum:g} to {FLATNESS.maximum:g}"


# --- The scaling relations over the seeds -------------------------------------------------------------------------------------------------------

def test_the_disk_length_rises_with_the_mass_at_the_size_mass_relations_slope(size_draws):
    """Mass-size: the fitted slope of the log of the drawn scale length against the log of the mass is the relation's slope over the same masses."""
    log_mass = np.array([math.log10(d.chain.stellar_mass_solar) for d in size_draws])
    drawn = np.array([math.log10(d.disk_scale_length_pc) for d in size_draws])
    law = np.array([math.log10(shape_draw.size_median_kpc(d.chain.stellar_mass_solar)) for d in size_draws])
    assert np.polyfit(log_mass, drawn, 1)[0] == pytest.approx(np.polyfit(log_mass, law, 1)[0], abs=0.03)


def test_the_median_size_residual_is_nil_in_every_mass_bin(size_draws):
    for low, high in MASS_BINS:
        bin_draws = [d for d in size_draws if low <= d.chain.stellar_mass_solar < high]
        offsets = [_size_offset_dex(d) for d in bin_draws]
        assert abs(statistics.mean(offsets)) < 0.03, _bin_name(low, high)
        sampling = 4.0 / math.sqrt(2.0 * (len(offsets) - 1))     # the standard deviation of a Gaussian sample of this size, to four standard errors
        assert statistics.stdev(offsets) == pytest.approx(
            statistics.mean(_size_limit_dex(d) for d in bin_draws) / SIGMAS, rel=sampling), _bin_name(low, high)


def test_the_halo_follows_the_stellar_to_halo_relation_in_every_mass_bin(draws):
    """Stellar-to-halo: the median of the log of the halo over the stars is the relation's median (`mass_chain.median_halo_mass`) at the bin's masses,
    to the relation's own scatter over the root of the count."""
    for low, high in MASS_BINS:
        bin_draws = [d for d in draws if low <= d.chain.stellar_mass_solar < high and d.chain.halo_mass_solar is not None]
        offsets = [math.log10(d.chain.halo_mass_solar / mass_chain.median_halo_mass(d.chain.stellar_mass_solar)) for d in bin_draws]
        assert len(bin_draws) >= 40, _bin_name(low, high)
        assert abs(statistics.mean(offsets)) < 4.0 * statistics.stdev(offsets) / math.sqrt(len(offsets)) + 0.03, _bin_name(low, high)


def test_the_black_hole_of_a_classical_bulge_follows_the_bulge_mass_relation_across_the_seeds(draws):
    """M_BH against the bulge mass: the fitted slope over the galaxies with a classical bulge is the relation's 1.16 within the fit's error and the
    scatter is the relation's 0.29 dex (Kormendy and Ho 2013, as `structure.black_hole_mass_solar` draws it)."""
    classical = [d for d in draws if d.structure.classical_bulge and d.bulge_to_total * d.chain.stellar_mass_solar > 1.0e8]
    assert len(classical) > 500
    log_bulge = np.array([math.log10(d.bulge_to_total * d.chain.stellar_mass_solar) for d in classical])
    log_hole = np.array([math.log10(galaxyStructure.black_hole_mass_solar(d.structure, d.bulge_to_total)) for d in classical])
    slope = np.polyfit(log_bulge, log_hole, 1)[0]
    assert slope == pytest.approx(1.16, abs=0.05)
    residual = log_hole - (math.log10(0.49e9) + 1.16 * (log_bulge - 11.0))
    assert abs(float(residual.mean())) < 0.03
    assert float(residual.std()) == pytest.approx(0.29, abs=0.03)


# --- The reference galaxies ---------------------------------------------------------------------------------------------------------------------

# `reports/galaxy-size-sources-2026-10-21.md`, section 1 (read from the papers: Bland-Hawthorn and Gerhard 2016 [R0422], and for Andromeda Courteau 2011,
# Dorman 2013, Dalcanton 2023, Sick 2014, Tamm 2012, Kafle 2018, Patel 2017 [R0909, R0979, R1145, R4084]): each quantity as the range its sources span.
REFERENCES = {
    "Milky Way": {"stellar mass": (4.0e10, 6.0e10), "disk length pc": (2100.0, 3100.0), "disk height pc": (250.0, 350.0), "halo mass": (0.8e12, 1.4e12)},
    "Andromeda": {"stellar mass": (8.6e10, 12.6e10), "disk length pc": (4800.0, 5900.0), "disk height pc": (690.0, 850.0), "halo mass": (0.6e12, 2.8e12)},
}
REACH_AT_LEAST = 0.01
"""float: The share of the seeds that must land inside a reference quantity's range: it must be reachable, not typical."""
JOINT_REACH_AT_LEAST = 0.001


@pytest.fixture(scope="module")
def reference_draws():
    """The draws at each reference galaxy's own stellar mass (the middle of its range, in the log)."""
    found = {}
    for name, ranges in REFERENCES.items():
        low, high = ranges["stellar mass"]
        found[name] = [shape_draw.draw_galaxy(seed, stellar_mass_solar=math.sqrt(low * high)) for seed in SEEDS]
    return found


def _quantities(drawn):
    return {"disk length pc": drawn.disk_scale_length_pc, "disk height pc": drawn.disk_scale_height_pc, "halo mass": drawn.chain.halo_mass_solar}


@pytest.mark.parametrize("name", sorted(REFERENCES))
def test_each_quantity_of_a_reference_galaxy_is_reachable_by_the_draw_at_its_mass(name, reference_draws):
    found = reference_draws[name]
    for quantity, (low, high) in REFERENCES[name].items():
        if quantity == "stellar mass":
            continue
        share = sum(1 for d in found if _holds(_quantities(d)[quantity], low, high)) / len(found)
        assert share >= REACH_AT_LEAST, f"{name}: {quantity} {low:g} to {high:g} is reached by {100 * share:.2f} percent of the seeds at its mass"


@pytest.mark.parametrize("name", sorted(REFERENCES))
def test_a_reference_galaxys_length_and_height_are_reachable_together(name, reference_draws):
    ranges = REFERENCES[name]
    found = reference_draws[name]
    share = sum(1 for d in found if _holds(d.disk_scale_length_pc, *ranges["disk length pc"]) and _holds(d.disk_scale_height_pc, *ranges["disk height pc"])) / len(found)
    assert share >= JOINT_REACH_AT_LEAST, f"{name}: the disk length and height together are reached by {100 * share:.2f} percent of the seeds at its mass"


def test_the_reference_galaxies_are_inside_the_spiral_familys_mass_limits():
    for name, ranges in REFERENCES.items():
        assert SPIRAL_MASS_LIMITS[0] <= ranges["stellar mass"][0] and ranges["stellar mass"][1] <= SPIRAL_MASS_LIMITS[1], name


# --- The dials' ranges as plausibility bounds (GEN.601, GEN.609) --------------------------------------------------------------------------------

def _args(**given):
    return types.SimpleNamespace(**{**{row.dest: None for row in dials.LEVERS}, **given})


SHAPE_LEVERS = ("mass", "compactness", "density", "disk_scale_length", "flatness", "bulge_share", "bulge_size", "halo_ratio", "stellar_halo_share",
                "arm_count", "pitch_angle_deg")
"""tuple: The dials whose pin changes the shape the seed draws; each is swept from its floor to its ceiling."""
SWEEP_SEEDS = SEEDS[:3]


def _sweep_values(row):
    low, high = (row.min_abs if row.min_abs and row.minimum < 0.0 else row.minimum), row.maximum
    if row.kind == "int":
        return sorted({int(low), int(round((low + high) / 2.0)), int(high)})
    if low > 0.0:
        return [low, math.sqrt(low * high), high]
    return [low, (low + high) / 2.0, high]


def _plans():
    """`(lever, value, seed, outcome)` of every pin of the sweep: the plan's shape and draw, or the refusal (`ScaleError` or `ValueError`) it gave."""
    from planetgen.galaxy.shape_draw import ScaleError
    found = []
    for name in SHAPE_LEVERS:
        row = dials.lever(name)
        for value in _sweep_values(row):
            for seed in SWEEP_SEEDS:
                try:
                    outcome = shape_draw.plan_shape(_args(**{row.dest: value}), seed)
                except (ScaleError, ValueError) as exc:
                    outcome = exc
                found.append((name, value, seed, outcome))
    return found


@pytest.fixture(scope="module")
def sweep():
    return _plans()


def test_every_dial_pinned_from_its_floor_to_its_ceiling_gives_a_finite_galaxy_or_a_refusal_that_names_the_range(sweep):
    refused = 0
    for name, value, _seed, outcome in sweep:
        if isinstance(outcome, Exception):
            refused += 1
            assert str(outcome), f"{name} {value:g}: a refusal with no sentence"
            continue
        shape, drawn = outcome
        for field in (drawn.disk_scale_length_pc, drawn.disk_scale_height_pc, drawn.bulge_scale_radius_pc, drawn.bulge_to_total, drawn.pitch_angle_deg,
                      drawn.chain.stellar_mass_solar):
            assert math.isfinite(field), f"{name} {value:g}"
        assert shape is not None
    assert refused < len(sweep) / 2, "most of the sweep was refused: the dials' ranges are not ranges of galaxies"


def test_a_galaxy_the_sweep_makes_keeps_every_value_inside_the_published_bands(sweep):
    """The family stays inside its bands whatever dial is pinned: the disk's scale length and the mass inside the Scale dial's rows, the height and the
    bulge radius inside the derived limits, the bulge share inside the dial, the arms and the pitch inside the pattern dial."""
    mass, length, share = dials.lever("mass"), dials.lever("disk_scale_length"), dials.lever("bulge_share")
    arms, pitch = dials.lever("arm_count"), dials.lever("pitch_angle_deg")
    for name, value, _seed, outcome in sweep:
        if isinstance(outcome, Exception):
            continue
        _shape, drawn = outcome
        where = f"{name} pinned at {value:g}"
        assert mass.minimum <= drawn.chain.stellar_mass_solar <= mass.maximum or name in ("density",), where
        assert length.minimum <= drawn.disk_scale_length_pc <= length.maximum or name in ("density",), where
        assert dials.check_derived("disk_scale_height_pc", drawn.disk_scale_height_pc) is None, where
        assert dials.check_derived("bulge_scale_radius_pc", drawn.bulge_scale_radius_pc) is None, where
        assert share.minimum <= drawn.bulge_to_total <= share.maximum, where
        assert arms.minimum <= drawn.arm_count <= arms.maximum, where
        assert pitch.minimum <= drawn.pitch_angle_deg <= pitch.maximum, where


def test_the_halo_warning_fires_only_when_the_disk_is_a_decade_from_the_usual_share_of_the_virial_radius(sweep):
    """GEN.606: `halo_report` warns exactly when the disk's half-light radius over the halo's virial radius is more than the limit from the usual
    share, computed here from the plan's own numbers."""
    fired = quiet = 0
    for name, value, _seed, outcome in sweep:
        if isinstance(outcome, Exception):
            continue
        _shape, drawn = outcome
        _line, warning = shape_draw.halo_report(drawn)
        if drawn.chain.halo_mass_solar is None:
            assert warning is None, f"{name} {value:g}: a warning with no halo"
            continue
        ratio = drawn.disk_scale_length_pc * shape_draw.EXPONENTIAL_HALF_LIGHT_OVER_SCALE / 1000.0 / drawn.chain.virial_radius_kpc
        offset = abs(math.log10(ratio / shape_draw.DISK_HALF_LIGHT_OVER_R200))
        assert (warning is not None) == (offset > shape_draw.DISK_HALF_LIGHT_OVER_R200_DEX), f"{name} {value:g}: offset {offset:.2f} dex"
        fired += warning is not None
        quiet += warning is None
    assert quiet > 0


def test_the_grid_warning_fires_only_below_its_thresholds(sweep):
    """GEN.625: `dials.shape_warnings` names exactly the three shapes below its limits (a scale height, a bulge radius, a crest)."""
    for height, bulge, crest in ((20.0, 20.0, 50.0), (dials.WARN_HEIGHT_PC / 2, 20.0, 50.0), (20.0, dials.WARN_BULGE_PC / 2, 50.0), (20.0, 20.0, dials.WARN_CREST / 2),
                                 (dials.WARN_HEIGHT_PC, dials.WARN_BULGE_PC, dials.WARN_CREST)):
        notes = dials.shape_warnings(crest, 0.0, height, bulge)
        expected = (height < dials.WARN_HEIGHT_PC) + (bulge < dials.WARN_BULGE_PC) + (crest < dials.WARN_CREST)
        assert len(notes) == expected, (height, bulge, crest)
    for name, value, seed, outcome in sweep:
        if isinstance(outcome, Exception) or name not in ("flatness", "bulge_size"):
            continue
        _shape, drawn = outcome
        if drawn.disk_scale_height_pc >= dials.WARN_HEIGHT_PC and drawn.bulge_scale_radius_pc >= dials.WARN_BULGE_PC:
            assert not [n for n in dials.shape_warnings(1.0e9, 0.0, drawn.disk_scale_height_pc, drawn.bulge_scale_radius_pc) if "scale height" in n or "bulge radius" in n]
