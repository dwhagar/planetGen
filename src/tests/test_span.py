"""ADM.29: span fills -- ranges, counts by prefix sums, the CLI and the Generate page."""
import argparse

import pytest

from planetgen.galaxy.geometry import ring_sector_count
from planetgen.galaxy.skeleton import GalaxyBounds
from planetgen.galaxy.span import Span, SpanError, parse_range
from planetgen.web import generate_page


def _bounds():
    # layer 0 reaches ring 6, layers +-1 reach ring 3, layers +-2 reach ring 1
    return GalaxyBounds({0: 6, 1: 3, -1: 3, 2: 1, -2: 1}, 4.0)


def test_parse_range():
    assert parse_range("3:7", "x") == (3, 7)
    assert parse_range("4", "x") == (4, 4)
    assert parse_range("-2:2", "x") == (-2, 2)
    for bad in ("", "1:2:3", "a:b", "3:"):
        with pytest.raises(SpanError):
            parse_range(bad, "x")


@pytest.mark.parametrize("kwargs", [
    {"rings": (5, 2)}, {"layers": (3, 1)}, {"slots": (0, 3)}, {"slots": (0, 3), "rings": (1, 2)},
    {"rings": (2, 2), "slots": (0, 10 ** 6)}, {"rings": (-1, 2)},
])
def test_span_refuses_nonsense(kwargs):
    with pytest.raises(SpanError):
        Span(**kwargs)


def test_count_matches_the_addresses_it_lists():
    bounds = _bounds()
    for span in (Span(), Span(rings=(1, 4)), Span(layers=(0, 1)), Span(rings=(2, 5), layers=(-1, 1)),
                 Span(layers=(2, 2)), Span(rings=(3, 3), slots=(2, 6)), Span(rings=(3, 3), slots=(ring_sector_count(3) - 2, 1))):
        listed = list(span.addresses(bounds))
        assert span.count(bounds) == len(listed) == len(set(listed))
        assert all(bounds.contains(ring, layer) for ring, layer, _slot in listed)


def test_whole_galaxy_is_the_outline_cell_count():
    bounds = _bounds()
    assert Span().count(bounds) == bounds.cell_count()


def test_slot_arc_wraps_through_zero():
    total = ring_sector_count(3)
    slots = [slot for _r, layer, slot in Span(rings=(3, 3), layers=(0, 0), slots=(total - 2, 1)).addresses(_bounds())]
    assert slots == [total - 2, total - 1, 0, 1]


def test_rings_beyond_a_layers_edge_are_left_out():
    listed = list(Span(rings=(2, 5), layers=(2, 2)).addresses(_bounds()))
    assert listed == []  # layer 2 ends at ring 1


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
    args = _parse(["--rings", "3:5", "--layers=-1:1"])
    assert args.span.rings == (3, 5) and args.span.layers == (-1, 1) and args.span.slots is None


@pytest.mark.parametrize("argv", [["--slots", "1:2"], ["--rings", "1:3", "--slot", "2"], ["--rings", "9:2"],
                                  ["--layers", "1", "--column"]])
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


def test_cylinder_extent_holds_the_face_neighbours_at_one_sector():
    import math
    from planetgen.galaxy.geometry import enumerate_sectors_within_radius, sector_address_at
    from planetgen.generation import run_galaxy
    edge = 4.0
    horizontal, layers, sphere = run_galaxy.cylinder_extent(1, 1, edge)
    assert horizontal == pytest.approx(1.385 * edge) and layers == 1
    assert sphere == pytest.approx(math.hypot(horizontal, 1.5 * edge))
    center = (60.0, 25.0, 0.0)
    here = sector_address_at(center, edge)
    can = [c for c in enumerate_sectors_within_radius(center, sphere, edge)
           if abs(c[1] - here[1]) <= layers and math.hypot(c[3] - center[0], c[4] - center[1]) <= horizontal]
    plane = [c for c in can if c[1] == here[1]]
    assert 5 <= len(plane) <= 7 and len(can) == len(plane) * 3


def test_cli_cylinder_options():
    args = _parse(["--center-sector", "3", "--cylinder-sectors", "2"])
    assert args.cylinder_sectors == 2 and args.cylinder_layers == 2 and args.radius_pc > 0
    assert _parse(["--center-sector", "3", "--cylinder-sectors", "2", "--cylinder-layers", "0"]).cylinder_layers == 0
    for bad in (["--cylinder-sectors", "2"], ["--center-sector", "3", "--cylinder-sectors", "2", "--radius-pc", "5"],
                ["--center-sector", "3", "--cylinder-layers", "1"], ["--center-sector", "3", "--cylinder-sectors", "999"],
                ["--rings", "1", "--cylinder-sectors", "2"]):
        with pytest.raises(SystemExit):
            _parse(bad)


def test_page_center_argv_takes_a_radius_in_sectors():
    argv, description = generate_page.center_argv({"center_by": "sector", "center_sector": "7",
                                                   "center_cylinder_sectors": "3", "center_cylinder_layers": "1"})
    assert argv == ["--center-sector", "7", "--cylinder-sectors", "3", "--cylinder-layers", "1"]
    assert "3 sectors across" in description
