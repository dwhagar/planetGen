"""GEN.495: the comet volatile population (`planetgen/physics/comet_population.py`, its data file
`physics/data/comet_population.json`, parts G1 sections 2.3, 3.1 and 3.2). The Solar System comets below are test targets
only: the draw must be able to reach them, and nothing in the file is centred on one."""

import copy
import json
import math
import concurrent.futures
import multiprocessing
import statistics

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.physics import comet_population as cp
from planetgen.physics.elements import DataError
from planetgen.util import draw

COMETS = {
    # percent of H2O; 67P bulk (Rubin 2019 Table 2), 103P (Dello Russo 2016, A'Hearn 2011), Hale-Bopp (Rob24 database)
    "67P": {"CO": 3.1, "CO2": 4.7, "CH4": 0.34, "C2H6": 0.29, "NH3": 0.67, "HCN": 0.14, "CH3OH": 0.21, "H2CO": 0.32,
            "H2S": 1.10},
    "103P": {"CO": 0.3, "CO2": 20.0},
    "Hale-Bopp": {"CO": 26.0, "CH4": 1.22, "HCN": 0.36, "H2S": 1.5},
}


def _document():
    with open(cp.DATA_FILE, encoding="utf-8") as handle:
        return json.load(handle)


def broken(change):
    data = copy.deepcopy(_document())
    change(data)
    return cp._parse(data)


def _sample(period_class, count, seed=1, **kwargs):
    with draw.bound(seed):
        return [cp.draw_volatiles(period_class, **kwargs) for _ in range(count)]


# --- the file states where each number comes from ---------------------------------------------------------------

def test_every_number_has_a_source_a_label_and_a_state():
    everything = list(cp.CARBON_OXIDES.values()) + list(cp.SPECIES.values())
    assert [s.name for s in everything][:2] == ["CO", "CO2"] and len(everything) == 9
    for species in everything:
        assert species.source in cp.SOURCES and species.label and species.state in (pt.VERIFIED, pt.PARTLY_VERIFIED)
        assert species.lo_percent < species.geo_mean_percent < species.hi_percent and species.sigma_dex > 0
        if species.state == pt.PARTLY_VERIFIED:
            assert species.limit, species.name


def test_no_solar_system_comet_is_in_the_data_file():
    text = json.dumps(_document()).lower()
    for name in ("67p", "churyumov", "hartley", "103p", "hale-bopp", "1p/", "borrelly"):
        assert name not in text, name
    assert "targets" not in _document()


def test_the_carbon_chain_flag_is_a_registered_table_with_every_class_gated():
    table = pt.REGISTRY.get("comet-carbon-chain")
    assert table is cp.CARBON_CHAIN_TABLE and set(table.classes()) == {"comet"}
    assert table.state_counts()[pt.OPEN_STATE] == 0
    for period, probability in cp.DEPLETED["probability"].items():
        shares = {row.name: p for row, p in table.probabilities(pt.Subject("comet", {"period_class": period}))}
        assert math.isclose(sum(shares.values()), 1.0, abs_tol=1e-12)
        assert math.isclose(shares[f"carbon-chain-depleted-{period}"], probability)
        assert len(shares) == 2


# --- the draw ---------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("period_class", cp.PERIOD_CLASSES)
def test_every_drawn_species_is_inside_its_observed_bounds(period_class):
    for v in _sample(period_class, 1500):
        assert list(v.percent)[:2] == ["CO", "CO2"] and set(v.percent) == {"CO", "CO2", *cp.SPECIES}
        for name, value in v.percent.items():
            s = cp.CARBON_OXIDES.get(name) or cp.SPECIES[name]
            assert s.lo_percent <= value <= s.hi_percent, (name, value)


def test_the_other_species_follow_their_truncated_log_normals_whatever_the_class():
    for name, species in cp.SPECIES.items():
        for period_class in cp.PERIOD_CLASSES:
            values = sorted(v.percent[name] for v in _sample(period_class, 3000, seed=7))
            ks = max(max(abs((i + 1) / len(values) - species.cdf(x)), abs(i / len(values) - species.cdf(x)))
                     for i, x in enumerate(values))
            assert ks < 0.035, (name, period_class, ks)


def test_co_and_co2_trade_off_and_the_class_moves_them():
    medians = {}
    for period_class in cp.PERIOD_CLASSES:
        sample = _sample(period_class, 4000, seed=3)
        co = [v.percent["CO"] for v in sample]
        co2 = [v.percent["CO2"] for v in sample]
        total = [a + b for a, b in zip(co, co2)]
        assert statistics.median(total) == pytest.approx(18.0, rel=0.12)
        medians[period_class] = (statistics.median(co), statistics.median(co2))
        low_co = [b for a, b in zip(co, co2) if a < statistics.median(co)]
        high_co = [b for a, b in zip(co, co2) if a >= statistics.median(co)]
        assert statistics.mean(low_co) > statistics.mean(high_co)
    assert medians["long_period"][0] > 3 * medians["jupiter_family"][0]
    assert medians["long_period"][1] < medians["jupiter_family"][1]
    assert medians["halley_type"] == pytest.approx(medians["jupiter_family"], rel=0.15)


def test_the_pooled_co_and_co2_match_the_population_table_within_its_scatter():
    sample = [v for c in cp.PERIOD_CLASSES for v in _sample(c, 4000, seed=11)]
    for name in ("CO", "CO2"):
        s = cp.CARBON_OXIDES[name]
        logs = [math.log10(v.percent[name]) for v in sample]
        assert 10 ** statistics.mean(logs) == pytest.approx(s.geo_mean_percent, rel=0.25)
        assert statistics.pstdev(logs) == pytest.approx(s.sigma_dex, abs=0.15)


def test_the_carbon_chain_flag_follows_the_class_fractions():
    for period_class, probability in cp.DEPLETED["probability"].items():
        sample = _sample(period_class, 6000, seed=5)
        share = sum(v.carbon_chain_depleted for v in sample) / len(sample)
        assert share == pytest.approx(probability, abs=0.025)
        assert {round(v.depleted_probability, 9) for v in sample} <= {round(probability, 9), round(1 - probability, 9)}


def test_an_interstellar_comet_draws_its_carbon_oxides_from_the_wider_ranges():
    sample = _sample("long_period", 3000, seed=13, interstellar=True)
    ranges = cp.INTERSTELLAR
    co = [v.percent["CO"] for v in sample]
    co2 = [v.percent["CO2"] for v in sample]
    assert min(co) >= 100 * ranges["co_over_h2o"]["lo"] and max(co) <= 100 * ranges["co_over_h2o"]["hi"]
    assert min(co2) >= 100 * ranges["co2_over_h2o"]["lo"] and max(co2) <= 100 * ranges["co2_over_h2o"]["hi"]
    assert max(co2) > cp.CARBON_OXIDES["CO2"].hi_percent and max(co) > cp.CARBON_OXIDES["CO"].hi_percent
    for name, species in cp.SPECIES.items():
        assert all(species.lo_percent <= v.percent[name] <= species.hi_percent for v in sample), name


def test_an_unknown_class_is_refused():
    with draw.bound(1):
        with pytest.raises(ValueError, match="period class"):
            cp.draw_volatiles("oort")


# --- the Solar System comets are reachable (test targets only) --------------------------------------------------

def test_each_measured_comet_lies_inside_the_hard_bounds_and_where_the_draw_has_density():
    for comet, values in COMETS.items():
        for name, value in values.items():
            s = cp.CARBON_OXIDES.get(name) or cp.SPECIES[name]
            assert s.lo_percent <= value <= s.hi_percent, (comet, name, value)
            if name not in cp.CARBON_OXIDES:
                assert s.cdf(value * 1.01) > 0.0 and s.cdf(value * 0.99) < 1.0, (comet, name, value)


def test_the_carbon_oxide_corners_103p_and_hale_bopp_are_reached_by_the_draw():
    sample = _sample("jupiter_family", 20000, seed=17) + _sample("long_period", 20000, seed=19)
    co_poor = [v for v in sample if v.percent["CO"] < 0.6 and v.percent["CO2"] > 14]
    co_rich = [v for v in sample if v.percent["CO"] > 14 and v.percent["CO2"] < 12]
    assert co_poor and co_rich



def test_67p_bulk_carbon_oxide_ratio_is_inside_the_central_99_percent_of_the_short_period_draw():
    # GEN.629: S = (CO + CO2)/H2O scatters as single comets do (0.2 dex), not as the median's uncertainty (0.1 dex)
    sample = _sample("jupiter_family", 10000, seed=29)
    totals = sorted(v.percent["CO"] + v.percent["CO2"] for v in sample)
    target = COMETS["67P"]["CO"] + COMETS["67P"]["CO2"]
    assert totals[50] < target < totals[-51]
    assert statistics.median(totals) == pytest.approx(cp.TOTAL_CARBON_OXIDE["median_percent"], rel=0.05)
    together = [v for v in sample if abs(v.percent["CO"] - 3.1) < 1.0 and abs(v.percent["CO2"] - 4.7) < 1.5]
    assert len(together) >= 10

# --- determinism ------------------------------------------------------------------------------------------------

def _worker(seed):
    with draw.bound(seed):
        out = []
        for period_class in cp.PERIOD_CLASSES:
            v = cp.draw_volatiles(period_class)
            out.append((sorted(v.percent.items()), v.carbon_chain_depleted))
        v = cp.draw_volatiles("long_period", interstellar=True)
        out.append((sorted(v.percent.items()), v.carbon_chain_depleted))
        return out


def test_a_draw_is_deterministic_from_its_seed_at_one_and_four_workers():
    seeds = list(range(8))
    serial = [_worker(s) for s in seeds]
    assert serial == [_worker(s) for s in seeds] and serial[0] != serial[1]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4, mp_context=multiprocessing.get_context("spawn")) as pool:
        assert list(pool.map(_worker, seeds)) == serial


def test_the_draw_reads_only_the_one_random_wrapper(monkeypatch):
    import random
    monkeypatch.setattr(random, "random", lambda: (_ for _ in ()).throw(AssertionError("module random used")))
    assert _worker(3) == _worker(3)


# --- the loader refuses -----------------------------------------------------------------------------------------

def test_a_species_without_a_source_or_with_a_state_it_does_not_know_is_refused():
    with pytest.raises(DataError, match="CH4.*source"):
        broken(lambda d: d["species"][0].update(source=""))
    with pytest.raises(DataError, match="CH4.*not in the file's sources"):
        broken(lambda d: d["species"][0].update(source="nobody99"))
    with pytest.raises(DataError, match="state"):
        broken(lambda d: d["species"][0].update(state="open"))
    with pytest.raises(DataError, match="H2S.*partly verified and names no limit"):
        broken(lambda d: d["species"][-1].update(limit=""))
    with pytest.raises(DataError, match="carbon_chain_depleted.*partly verified and names no limit"):
        broken(lambda d: d["carbon_chain_depleted"].update(limit=""))


def test_a_geometric_mean_outside_its_bounds_or_a_bad_number_is_refused():
    with pytest.raises(DataError, match="geometric mean.*not inside"):
        broken(lambda d: d["species"][0].update(geo_mean_percent=5.0))
    with pytest.raises(DataError, match="not a positive number"):
        broken(lambda d: d["species"][0].update(sigma_dex=0))
    with pytest.raises(DataError, match="not a fraction"):
        broken(lambda d: d["co_fraction"]["median"].update(long_period=1.0))
    with pytest.raises(DataError, match="probability for halley_type"):
        broken(lambda d: d["carbon_chain_depleted"]["probability"].pop("halley_type"))


def test_missing_parts_and_a_duplicated_or_misplaced_species_are_refused():
    with pytest.raises(DataError, match="states no interstellar"):
        broken(lambda d: d.pop("interstellar"))
    with pytest.raises(DataError, match="interstellar must hold"):
        broken(lambda d: d["interstellar"].pop("co2_over_h2o"))
    with pytest.raises(DataError, match="appears twice"):
        broken(lambda d: d["species"].append(copy.deepcopy(d["species"][0])))
    with pytest.raises(DataError, match="not drawn as independent|hyper-volatile block"):
        broken(lambda d: d["species"].append(dict(copy.deepcopy(d["species"][0]), name="CO")))
    with pytest.raises(DataError, match="carbon_oxides must hold"):
        broken(lambda d: d["carbon_oxides"].pop("CO2"))


# --- GEN.595: the median from the disc zone ---------------------------------------------------------------------

def test_a_zone_sets_the_median_carbon_oxides_and_the_class_no_longer_does():
    """A comet formed in a CO2-rich zone keeps CO2 over CO whatever its class, and one from a CO-rich zone the reverse; the scatter stays the population's."""
    co2_zone = _sample("long_period", 3000, seed=4, ice_ratios=(0.33, 0.025))
    co_zone = _sample("long_period", 3000, seed=4, ice_ratios=(0.025, 0.33))
    assert statistics.median(v.percent["CO2"] for v in co2_zone) == pytest.approx(33.0, rel=0.15)
    assert statistics.median(v.percent["CO"] for v in co2_zone) == pytest.approx(2.5, rel=0.2)
    assert statistics.median(v.percent["CO"] for v in co_zone) == pytest.approx(33.0, rel=0.3)
    assert statistics.median(v.percent["CO2"] for v in co_zone) == pytest.approx(2.5, rel=0.25)
    by_class = [statistics.median(v.percent["CO"] for v in _sample(c, 2000, seed=5, ice_ratios=(0.33, 0.1))) for c in cp.PERIOD_CLASSES]
    assert max(by_class) / min(by_class) < 1.2
    spread = statistics.pstdev(math.log10(v.percent["CO2"] + v.percent["CO"]) for v in co2_zone)
    assert spread == pytest.approx(cp.TOTAL_CARBON_OXIDE["sigma_dex"], rel=0.25)


def test_the_bounds_move_with_the_zones_median_and_a_zone_beyond_the_observed_range_is_drawn():
    """A Borisov-like zone (CO over water of 1.7) lies past the 30 percent the Solar System's comets reach; its draws stay as far from its median as the observed ones from theirs."""
    sample = _sample("long_period", 2000, seed=6, ice_ratios=(0.33, 1.67))
    co = [v.percent["CO"] for v in sample]
    assert statistics.median(co) == pytest.approx(167.0, rel=0.3)
    assert max(co) <= 167.0 * cp.CARBON_OXIDES["CO"].hi_percent / cp.CARBON_OXIDES["CO"].geo_mean_percent * 1.0001


def test_a_star_with_no_carbon_freezes_no_carbon_oxide_and_no_ratios_keeps_the_population():
    assert cp.draw_volatiles("long_period", ice_ratios=(0.0, 0.0)).percent["CO"] == 0.0
    with draw.bound(8):
        plain = cp.draw_volatiles("jupiter_family")
    with draw.bound(8):
        assert cp.draw_volatiles("jupiter_family", ice_ratios=None).percent == plain.percent
