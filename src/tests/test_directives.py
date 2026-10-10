"""GEN.96: generation directives -- parsing, the redraw loop, the Generate page."""
import pytest

from planetgen.generation import directives
from planetgen.util import draw
from planetgen.web import generate_page


class _Star:
    def __init__(self, type_, yerkes="V"):
        self.type, self.yerkes_class = type_, yerkes


class _System:
    def __init__(self, type_, hab=0, yerkes="V"):
        self.star, self.hab_count = _Star(type_, yerkes), hab


class _Entry:
    def __init__(self, *args, **kwargs):
        self.star_system = _System(*args, **kwargs)


class _Sector:
    def __init__(self, entries):
        self.entries = entries


def test_parse_reads_each_kind_and_keeps_the_larger_repeat():
    d = directives.parse(["systems>=5", "habitable >= 2", "type:g>=3", "type:WD>=1", "systems>=7"])
    assert d.minimums == {"systems": 7, "habitable": 2, "type:G": 3, "type:wd": 1}


@pytest.mark.parametrize("text", ["systems=5", "type:Z>=1", "planets>=2", "systems>=-1", ""])
def test_parse_refuses_what_it_cannot_read(text):
    with pytest.raises(directives.DirectiveError):
        directives.parse([text])


def test_digest_ignores_order():
    assert directives.parse(["systems>=5", "habitable>=1"]).digest() == directives.parse(
        ["habitable>=1", "systems>=5"]).digest()


def test_measure_counts_systems_habitable_and_types():
    sector = _Sector([_Entry("G2", hab=1), _Entry("K5", hab=2), _Entry("G0"), _Entry("D", yerkes="D")])
    got = directives.measure(sector)
    assert got["systems"] == 4 and got["habitable"] == 3
    assert got["type:G"] == 2 and got["type:K"] == 1 and got["type:wd"] == 1


def test_no_directive_is_one_draw():
    calls = []
    sector, outcome = directives.generate(directives.Directive(), lambda: calls.append(1) or _Sector([]))
    assert len(calls) == 1 and outcome.status == directives.MET_NATURALLY


def test_redraws_until_met_and_is_repeatable():
    seed = bytes(range(16))

    def build():
        return _Sector([_Entry("G2") for _ in range(draw.randint(0, 6))])

    d = directives.parse(["systems>=5"])
    results = []
    for _ in range(2):
        with draw.bound(1234):
            results.append(directives.generate(d, build, galaxy_seed=seed, address=(3, 0, 1)))
    for sector, outcome in results:
        assert len(sector.entries) >= 5 and outcome.status in (directives.MET_NATURALLY, directives.MET_AFTER_ATTEMPTS)
    assert len(results[0][0].entries) == len(results[1][0].entries)
    assert results[0][1].attempts == results[1][1].attempts


def test_unmet_keeps_the_closest_draw():
    sizes = iter([1, 3, 2, 3])
    sector, outcome = directives.generate(
        directives.parse(["systems>=9"]), lambda: _Sector([_Entry("G2")] * next(sizes)), max_attempts=4)
    assert outcome.status == directives.UNMET and outcome.attempts == 4
    assert len(sector.entries) == 3 and outcome.missed == {"systems": (9, 3)}
    assert "systems 3/9" in outcome.describe()


def test_generate_page_builds_directive_argv():
    assert generate_page.directive_argv({"directive_systems": "5", "directive_type_g": "2", "directive_habitable": ""}) == [
        "--directive", "systems>=5", "--directive", "type:G>=2"]
    with pytest.raises(generate_page.FormError):
        generate_page.directive_argv({"directive_systems": "-3"})
