# tests/test_ux_compact_pairs.py

"""
UX.186: the pages and the API name the compact binaries of a stored member tree (GEN.272): the kind of an accreting pair
(high-mass, Be/X-ray or low-mass X-ray binary), the microquasar flag, and the spider kind (black widow or redback) of an
ablating pulsar. They are read from what the tree stores (node interaction and gas-source bits, member flag bits, the
companion's mass). GEN.952: the members keep their radius, so an accreting pair also carries its X-ray luminosity, the figure
`MemberTree.x_ray_luminosities` works out when the tree is drawn; a spider pair, and a tree stored without radii, have none.
"""

import pytest

from planetgen.db import query, scatter_systems, store
from planetgen.db.compact_pairs import compact_pairs
from planetgen.generation import multistar
from planetgen.web.lib.fmt import format_power_w
from planetgen.web.scatter_members_view import members_view
from tests.test_api import client  # noqa: F401
from tests.test_galactic_motion import SCATTER_RADIUS_PC, _insert_scatter_point, _seed_scatter_grid
from tests.test_gen_system_motion import EPOCH
from tests.test_uncharted_phenomenon import db_app  # noqa: F401

POINT = (SCATTER_RADIUS_PC, 0.0, 0.0)


RADIUS_AU = {"neutron-star": 6.7e-11, "black-hole": 2.0e-10}
"""dict: A compact object's radius in AU for the tests (10 km; a few Schwarzschild radii); a star's is 0.01 AU."""


def _member(number, kind, mass, flags=0):
    return multistar.Member(number, 1, number, kind, None if kind in multistar.COMPACT_KINDS else "V", flags, mass, mass,
                            RADIUS_AU.get(kind, 0.01), number + 1)


def _node(interaction, flags=0):
    return multistar.Node(1, None, 0, 10.0, 0.5, 0.5 * 149597870.7, 0.1, 0.0, 0.0, 0.0, 0.0, 0.01, "close", interaction, flags)


def _tree(compact, companion, node):
    return multistar.MemberTree((compact, companion), (node,), 5e9)


CASES = {
    "high-mass-xray-binary": (_member(0, "neutron-star", 1.4), _member(1, "main-sequence", 25.0),
                              _node("accreting", multistar.NODE_WIND_FED), "High-mass X-ray binary"),
    "be-xray-binary": (_member(0, "neutron-star", 1.4), _member(1, "main-sequence", 12.0),
                       _node("accreting", multistar.NODE_DISC_FED), "Be/X-ray binary"),
    "low-mass-xray-binary": (_member(0, "black-hole", 8.0), _member(1, "main-sequence", 0.8),
                             _node("accreting"), "Low-mass X-ray binary"),
    "black-widow": (_member(0, "neutron-star", 1.4, multistar.MEMBER_ABLATING), _member(1, "main-sequence", 0.03),
                    _node("detached"), "Black widow pulsar"),
    "redback": (_member(0, "neutron-star", 1.4, multistar.MEMBER_ABLATING), _member(1, "main-sequence", 0.3),
                _node("detached"), "Redback pulsar"),
}


def _rows(tree):
    members = [{"member_no": m.member_no, "kind": m.kind, "mass_solar": m.mass_solar, "parent_node_no": m.parent_node_no,
                "flags": m.flags, "radius_au": m.radius_au} for m in tree.members]
    nodes = [{"node_no": n.node_no, "interaction": n.interaction, "flags": n.flags, "separation_au": n.separation_au,
              "eccentricity": n.eccentricity} for n in tree.nodes]
    return members, nodes


@pytest.mark.parametrize("kind", sorted(CASES))
def test_each_kind_of_compact_pair_is_named_from_what_the_tree_stores(kind):
    compact, companion, node, _label = CASES[kind]
    [pair] = compact_pairs(*_rows(_tree(compact, companion, node)))
    assert pair["kind"] == kind and (pair["compact"], pair["companion"]) == (0, 1) and pair["microquasar"] is False


def test_an_accreting_pair_carries_the_x_ray_luminosity_the_drawn_tree_works_out_and_a_spider_has_none():
    for kind, (compact, companion, node, _label) in CASES.items():
        tree = _tree(compact, companion, node)
        [pair] = compact_pairs(*_rows(tree))
        if kind in ("black-widow", "redback"):
            assert pair["x_ray_luminosity_w"] is None
        else:
            assert pair["x_ray_luminosity_w"] == tree.x_ray_luminosities()[node.node_no] > 0.0


def test_a_member_stored_without_a_radius_leaves_the_luminosity_empty():
    compact, companion, node, _label = CASES["high-mass-xray-binary"]
    members, nodes = _rows(_tree(compact, companion, node))
    members[1]["radius_au"] = None
    assert compact_pairs(members, nodes)[0]["x_ray_luminosity_w"] is None


def test_a_pair_that_is_not_compact_or_not_accreting_is_left_out():
    star = _member(0, "main-sequence", 1.0)
    assert compact_pairs(*_rows(_tree(star, _member(1, "main-sequence", 0.8), _node("accreting")))) == []
    assert compact_pairs(*_rows(_tree(_member(0, "neutron-star", 1.4), _member(1, "main-sequence", 1.0),
                                      _node("detached")))) == []


def test_the_microquasar_flag_is_the_compact_objects_not_the_companions():
    node = _node("accreting", multistar.NODE_WIND_FED)
    jets = _tree(_member(0, "black-hole", 9.0, multistar.MEMBER_MICROQUASAR), _member(1, "main-sequence", 30.0), node)
    assert compact_pairs(*_rows(jets))[0]["microquasar"] is True
    flagged_companion = _tree(_member(0, "black-hole", 9.0), _member(1, "main-sequence", 30.0, multistar.MEMBER_MICROQUASAR), node)
    assert compact_pairs(*_rows(flagged_companion))[0]["microquasar"] is False


def _stored(mysql_config, tree):
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "bright_stars", POINT)
            scatter_systems.store_member_tree(conn, "bright_stars", row_id, tree, EPOCH)
        return row_id, query.scatter_member_tree(conn, "bright_stars", row_id)
    finally:
        conn.close()


@pytest.mark.parametrize("kind", sorted(CASES))
def test_the_uncharted_star_page_and_the_api_name_each_kind(kind, db_app, client, mysql_config):  # noqa: F811
    compact, companion, node, label = CASES[kind]
    row_id, found = _stored(mysql_config, _tree(compact, companion, node))
    assert [pair["kind"] for pair in found["compact_pairs"]] == [kind]
    html = db_app.test_client().get(f"/systems/uncharted/{row_id}").get_data(as_text=True)
    assert 'id="pairs-heading"' in html and label in html
    assert [pair["kind"] for pair in client.get(f"/api/uncharted-systems/{row_id}").get_json()["compact_pairs"]] == [kind]


def test_the_stored_radius_gives_the_page_and_the_api_the_x_ray_luminosity(db_app, client, mysql_config):  # noqa: F811
    compact, companion, node, _label = CASES["high-mass-xray-binary"]
    tree = _tree(compact, companion, node)
    row_id, found = _stored(mysql_config, tree)
    assert [m["radius_au"] for m in found["members"]] == [compact.radius_au, companion.radius_au]
    watts = tree.x_ray_luminosities()[node.node_no]
    assert found["compact_pairs"][0]["x_ray_luminosity_w"] == pytest.approx(watts)
    assert members_view(found)["pairs"][0][5] == format_power_w(found["compact_pairs"][0]["x_ray_luminosity_w"])
    html = db_app.test_client().get(f"/systems/uncharted/{row_id}").get_data(as_text=True)
    assert "X-ray luminosity" in html
    api = client.get(f"/api/uncharted-systems/{row_id}").get_json()
    assert api["compact_pairs"][0]["x_ray_luminosity_w"] == pytest.approx(watts)


def test_the_page_says_when_an_x_ray_binary_is_a_microquasar(db_app, mysql_config):  # noqa: F811
    node = _node("accreting", multistar.NODE_WIND_FED)
    jets = _tree(_member(0, "black-hole", 9.0, multistar.MEMBER_MICROQUASAR), _member(1, "main-sequence", 30.0), node)
    quiet = _tree(_member(0, "black-hole", 9.0), _member(1, "main-sequence", 30.0), node)
    assert members_view(_stored(mysql_config, jets)[1])["pairs"][0][4] == "Yes"
    assert members_view(_stored(mysql_config, quiet)[1])["pairs"][0][4] == "No"


def test_a_system_with_no_compact_pair_has_no_compact_binaries_panel(db_app, mysql_config):  # noqa: F811
    star = _tree(_member(0, "main-sequence", 1.0), _member(1, "main-sequence", 0.8), _node("detached"))
    row_id, found = _stored(mysql_config, star)
    assert found["compact_pairs"] == [] and members_view(found)["pairs"] == []
    assert 'id="pairs-heading"' not in db_app.test_client().get(f"/systems/uncharted/{row_id}").get_data(as_text=True)
