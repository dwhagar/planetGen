"""ADM.29: span fills -- ranges, counts by prefix sums, the CLI and the Generate page."""
import argparse

import pytest

from planetgen.galaxy import spherical
from planetgen.galaxy.skeleton import GalaxyBounds
from planetgen.galaxy.span import Span, SpanError, parse_range
from planetgen.galaxy.spherical import wedge_count
from planetgen.web import generate_page


def _bounds():
    # the equator band reaches shell 6, bands +-1 reach shell 4, bands +-2 reach shell 3, +-3 reach shell 3
    return GalaxyBounds([(0, 0), (1, 1), (2, 2), (3, 3), (4, 3), (5, 2), (6, 1)], 4.0)


def test_parse_range():
    assert parse_range("3:7", "x") == (3, 7)
    assert parse_range("4", "x") == (4, 4)
    assert parse_range("-2:2", "x") == (-2, 2)
    for bad in ("", "1:2:3", "a:b", "3:"):
        with pytest.raises(SpanError):
            parse_range(bad, "x")


@pytest.mark.parametrize("kwargs", [
    {"shells": (5, 2)}, {"bands": (3, 1)}, {"wedges": (0, 3)}, {"wedges": (0, 3), "shells": (1, 2)},
    {"shells": (2, 2), "bands": (0, 0), "wedges": (0, 10 ** 6)}, {"shells": (-1, 2)},
])
def test_span_refuses_nonsense(kwargs):
    with pytest.raises(SpanError):
        Span(**kwargs)


def test_count_matches_the_addresses_it_lists():
    bounds = _bounds()
    arc_wrap = (wedge_count(3, 0) - 2, 1)
    for span in (Span(), Span(shells=(1, 4)), Span(bands=(0, 1)), Span(shells=(2, 5), bands=(-1, 1)),
                 Span(bands=(4, 4)), Span(shells=(3, 3), bands=(0, 0), wedges=(2, 6)),
                 Span(shells=(3, 3), bands=(0, 0), wedges=arc_wrap)):
        listed = list(span.addresses(bounds))
        assert span.count(bounds) == len(listed) == len(set(listed))
        assert all(bounds.contains(shell, band) for shell, band, _wedge in listed)
        assert all(0 <= wedge < wedge_count(shell, band) for shell, band, wedge in listed)


def test_whole_galaxy_is_the_outline_cell_count():
    bounds = _bounds()
    assert Span().count(bounds) == bounds.cell_count()
    assert bounds.cell_count() == sum(spherical.cells_in_band_range(shell, high) for shell, high in
                                      [(0, 0), (1, 1), (2, 2), (3, 3), (4, 3), (5, 2), (6, 1)])


def test_wedge_arc_wraps_through_zero():
    total = wedge_count(3, 0)
    wedges = [wedge for _s, _b, wedge in Span(shells=(3, 3), bands=(0, 0), wedges=(total - 2, 1)).addresses(_bounds())]
    assert wedges == [total - 2, total - 1, 0, 1]


def test_shells_beyond_a_bands_edge_are_left_out():
    listed = list(Span(shells=(5, 6), bands=(3, 3)).addresses(_bounds()))
    assert listed == []  # band 3 ends at shell 4


def test_a_span_lists_bands_north_to_south_a_shell_at_a_time():
    listed = list(Span(shells=(1, 2), bands=(-1, 1)).addresses(_bounds()))
    assert [shell for shell, _b, _w in listed] == sorted(shell for shell, _b, _w in listed)
    shell_one = [band for shell, band, _w in listed if shell == 1]
    assert shell_one == sorted(shell_one, reverse=True)


def _parse(argv):
    from planetgen.cli.generate import add_galaxy_arguments, add_shared_generation_options, validate_galaxy_args, validate_shared_generation_args
    parser = argparse.ArgumentParser(prefix_chars="-+")
    add_shared_generation_options(parser)
    add_galaxy_arguments(parser)
    args = parser.parse_args(argv)
    validate_shared_generation_args(args, parser)
    args.command = "galaxy"
    validate_galaxy_args(args, parser)
    return args


def test_cli_span_options_build_a_span():
    args = _parse(["--shells", "3:5", "--bands=-1:1"])
    assert args.span.shells == (3, 5) and args.span.bands == (-1, 1) and args.span.wedges is None


@pytest.mark.parametrize("argv", [["--wedges", "1:2"], ["--shells", "1:3", "--wedges", "2"], ["--shells", "9:2"],
                                  ["--bands", "1", "--shell", "3"]])
def test_cli_refuses_a_bad_span(argv):
    with pytest.raises(SystemExit):
        _parse(argv)


def test_page_builds_span_argv():
    argv, description = generate_page.span_argv({"span_rings": "3:5", "span_layers": "-1:1", "span_limit": "40"})
    assert argv == ["--rings=3:5", "--layers=-1:1", "--limit", "40"]
    assert "rings 3 to 5" in description and "layers -1 to 1" in description
    assert generate_page.span_argv({"span_layers": "0", "whole_span": "1"})[0] == ["--layers=0", "--yes"]
    for bad in ({}, {"span_slots": "1:2"}, {"span_rings": "x"}):
        with pytest.raises(generate_page.FormError):
            generate_page.span_argv(bad)


def test_ball_extent_holds_the_face_neighbours_at_one_sector():
    import math
    from planetgen.generation import run_galaxy
    edge = 4.0
    radius = run_galaxy.ball_extent(1, edge)
    assert radius == pytest.approx(1.75 * edge)
    for shell, band, wedge in ((5, 0, 0), (5, 2, 3), (20, 0, 3), (60, 2, 0)):
        center = spherical.sector_position_pc(shell, band, wedge, edge)
        face = spherical.neighbor_addresses(shell, band, wedge)
        reached = {c[:3] for c in spherical.enumerate_sectors_within_radius(center, radius, edge)}
        assert set(face) <= reached
        assert all(math.dist(center, spherical.sector_position_pc(*a, edge)) <= radius for a in face)


def test_page_center_argv_takes_a_radius_in_sectors():
    argv, description = generate_page.center_argv({"center_by": "sector", "center_sector": "7",
                                                   "center_cylinder_sectors": "3", "center_cylinder_layers": "1"})
    assert argv == ["--center-sector", "7", "--cylinder-sectors", "3", "--cylinder-layers", "1"]
    assert "3 sectors across" in description
