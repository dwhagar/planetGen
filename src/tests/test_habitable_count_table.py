"""The habitable bodies per system that sizes the directive fill (GEN.1245): `physics/data/habitable_count.json`, written by
`scripts/measure_habitable_count.py`, read by `generation.directives`."""
import json

import pytest

from planetgen import tuning
from planetgen.generation import directives


def _file():
    with open(directives.COUNT_FILE, encoding="utf-8") as handle:
        return json.load(handle)


def test_the_table_is_the_measured_histogram_and_names_what_it_was_measured_on():
    data = _file()
    assert data["systems"] >= 2000, "a few thousand systems, or the tail is noise"
    assert sum(data["counts"]) == data["systems"] and all(n >= 0 for n in data["counts"])
    assert isinstance(data["epoch"], int) and isinstance(data["seed"], int)
    assert directives.HABITABLE_COUNT_PMF == pytest.approx([n / data["systems"] for n in data["counts"]])
    assert sum(directives.HABITABLE_COUNT_PMF) == pytest.approx(1.0)
    assert directives.HABITABLE_COUNT_PMF[-1] > 0.0, "the list ends at the largest count met"


def test_the_table_has_the_spread_of_a_life_class_body_count():
    """Since epoch 121 a habitable body is a body of a life class: most systems hold at least one, a few hold many."""
    pmf = directives.HABITABLE_COUNT_PMF
    assert 0.1 < pmf[0] < 0.5, "systems with none"
    assert pmf[0] < 1.0 - pmf[0], "most systems hold one"
    assert len(pmf) >= 6, "a system can hold six or more"
    assert 0.8 < sum(i * share for i, share in enumerate(pmf)) < 4.0, "the mean count per system"


def test_a_bad_file_is_refused(tmp_path):
    path = tmp_path / "habitable_count.json"
    path.write_text(json.dumps({"counts": [3, 2], "systems": 9}))
    with pytest.raises(Exception, match="does not hold a histogram"):
        directives._load_count_pmf(str(path))


def test_the_compound_tail_agrees_with_drawing_systems_from_the_table():
    """The recursion of the impossible check and a straight sum of draws from the same table give one law."""
    pmf = directives.HABITABLE_COUNT_PMF
    systems, need = 3, 6
    exact = directives._fixed_tail(systems, need)
    dist = {0: 1.0}
    for _ in range(systems):
        step = {}
        for have, mass in dist.items():
            for k, share in enumerate(pmf):
                step[have + k] = step.get(have + k, 0.0) + mass * share
        dist = step
    total = sum(mass for have, mass in dist.items() if have >= need)
    assert exact == pytest.approx(total, rel=1e-9)
    # the compound-Poisson law at a mean of `systems` systems puts the same order on the same minimum
    assert directives._compound_poisson_tail(float(systems), need) > 0.0


def test_a_request_the_forced_fill_can_reach_is_never_called_impossible():
    """The fill builds as many systems as the minimum asks, so a sector holding that many systems is never refused."""
    for want in (1, 6, 12, 40):
        assert directives.impossible(directives.parse([f"habitable>={want}"]), systems=want) is None


def test_six_in_three_systems_is_common_not_rare():
    """The reason for GEN.1245: the old table put six in three systems at about 1 in 10,000; the measured table puts it near one
    in four to one in ten, which is what the sectors draw (TEST.555 part 3)."""
    assert 0.02 < directives._fixed_tail(3, 6) < 0.5
