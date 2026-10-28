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
        self.star, self.habitable_class_count = _Star(type_, yerkes), hab


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


def test_a_habitable_minimum_missed_for_k_draws_takes_the_forced_fill_and_says_so():
    """GEN.180: a directive that drawing cannot meet is met by the forced fill, marked met_forced."""
    counts = []

    def build_forced(count):
        counts.append(count)
        return _Sector([_Entry("G2", hab=1) for _ in range(count)])

    sector, outcome = directives.generate(
        directives.parse(["habitable>=4"]), lambda: _Sector([_Entry("G2", hab=0)] * 2), max_attempts=5,
        build_forced=build_forced)
    assert outcome.status == directives.MET_FORCED and counts == [4] and outcome.missed == {}
    assert len(sector.entries) == 4 and "forced fill" in outcome.describe() and "met_forced" in outcome.describe()
    assert outcome.attempts == 6, "the K draws, then the forced one"


def test_the_forced_fill_is_not_tried_when_the_draws_met_it_or_for_other_minimums():
    forced = []
    build = lambda: _Sector([_Entry("G2", hab=1)] * 3)       # noqa: E731
    sector, outcome = directives.generate(directives.parse(["habitable>=2"]), build, build_forced=forced.append)
    assert outcome.status == directives.MET_NATURALLY and not forced
    sector, outcome = directives.generate(directives.parse(["systems>=9"]), build, max_attempts=3, build_forced=forced.append)
    assert outcome.status == directives.UNMET and not forced, "only habitable worlds can be forced"


def test_a_forced_fill_that_still_misses_another_minimum_is_unmet_with_the_nearest_draw():
    sector, outcome = directives.generate(
        directives.parse(["habitable>=2", "systems>=9"]), lambda: _Sector([_Entry("G2")] * 3), max_attempts=2,
        build_forced=lambda count: _Sector([_Entry("G2", hab=1)] * count))
    assert outcome.status == directives.UNMET and "habitable" not in outcome.missed
    assert outcome.missed == {"systems": (9, 2)} and len(sector.entries) == 2, "the forced sector is nearer"


def test_the_forced_fill_is_repeatable_for_a_seed():
    seed = bytes(range(16))

    def build_forced(count):
        return _Sector([_Entry("G2", hab=1) for _ in range(count + draw.randint(0, 3))])

    sizes = []
    for _ in range(2):
        with draw.bound(99):
            sector, outcome = directives.generate(
                directives.parse(["habitable>=3"]), lambda: _Sector([_Entry("G2")]), galaxy_seed=seed,
                address=(1, 2, 3), max_attempts=3, build_forced=build_forced)
        sizes.append(len(sector.entries))
        assert outcome.status == directives.MET_FORCED
    assert sizes[0] == sizes[1]


def test_the_compound_poisson_tail_matches_the_tables_of_the_design_note():
    """sampling-backfill-and-resume.md 8.1: P(at least N habitable worlds) at 6.3 and 2.0 expected systems."""
    assert directives._compound_poisson_tail(6.3, 3) == pytest.approx(0.327, abs=0.003)
    assert directives._compound_poisson_tail(6.3, 5) == pytest.approx(0.104, abs=0.003)
    assert directives._compound_poisson_tail(2.0, 2) == pytest.approx(0.155, abs=0.003)
    assert directives._compound_poisson_tail(0.35, 1) == pytest.approx(0.082, abs=0.003)
    assert directives._poisson_tail(6.3, 0) == 1.0 and directives._poisson_tail(0.0, 1) == 0.0
    assert directives._poisson_tail(6.3, 500) == 0.0, "an underflowing tail ends, it does not loop"


@pytest.mark.parametrize("minimums, systems, mean, quoted", [
    (["systems>=5"], 3, None, "systems >= 5 cannot be met in a sector of 3 systems"),
    (["habitable>=6"], 1, None, "habitable >= 6 cannot be met"),
    (["type:G>=40"], 10, None, "type:G >= 40 cannot be met"),
    (["systems>=500"], None, 6.3, "about 6.3 systems"),
])
def test_an_impossible_request_is_refused_with_the_reason(minimums, systems, mean, quoted):
    reason = directives.impossible(directives.parse(minimums), systems=systems, mean=mean)
    assert reason and quoted in reason and "raise the number of systems" in reason


@pytest.mark.parametrize("minimums, systems, mean", [
    (["systems>=5"], 5, None), (["habitable>=2"], 1, None), (["habitable>=10"], None, 6.3),
    (["systems>=5"], None, 0.35), (["type:wd>=3", "habitable>=3"], 4, None), ([], 1, None)])
def test_a_request_a_sector_can_reach_is_not_refused(minimums, systems, mean):
    assert directives.impossible(directives.parse(minimums), systems=systems, mean=mean) is None
    assert directives.impossible(directives.parse(["systems>=5000"])) is None, "no count known: nothing is judged"


def test_the_command_line_refuses_an_impossible_directive_before_drawing(capsys, monkeypatch):
    from planetgen.cli import generate as generate_cli
    for argv, quoted in ((["sector", "--num-systems", "3", "--directive", "systems>=9", "--estimate-only"],
                          "systems >= 9 cannot be met in a sector of 3 systems"),
                         (["sectors", "--shell", "1", "--directive", "systems>=99999", "--estimate-only"],
                          "cannot be met")):
        monkeypatch.setattr("sys.argv", ["planetgen", *argv])
        with pytest.raises(SystemExit) as stopped:
            generate_cli.main()
        assert stopped.value.code == 2
        assert quoted in capsys.readouterr().err
