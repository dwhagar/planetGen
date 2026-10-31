# planetgen/physics/hierarchy_stability.py

"""
Stability of a tree of stars level by level (GEN.1016, GEN.1006 part P10; docs/design/systems/quadruple-stability-3p1.md)
========================================================================================================================

`physics/stability.py` tests one triple. A tree of four or more stars is several triples nested in each other, and
two shapes of four stars need more than that:

- **A triple level** is a node with a sub-node for one child: the sub-node's two sides are the inner pair, the other child
  is the third body. It must clear `level_margin` times the axis ratio the triple criterion asks for (`stability.criterion`
  picks Equation 4, Mardling and Aarseth or Eggleton and Kiseleva by the masses). The margin is 1.0, none, until a 2+2
  network is ported: the nested test is right for 88 percent of 2+2 systems (Vynatheya, Mardling and Hamers 2023, Table 3).
- **A 3+1 quadruple** is an inner pair, a third star about it and a fourth star about all three. Two nested triples
  are right for 66 percent of such systems and call too many of them stable, because the outer star pumps the eccentricity of
  the middle orbit. Both of its triple levels must pass at their own margin, and then the network of the same paper
  (`quadruple_stability`) must call it stable with both axis ratios loosened by `axis_margin`, which keeps a margin on the
  systems the network calls stable (its precision for "stable" was 0.55 on our own N-body check without it, 0.89 with 1.3). Inputs
  outside the range the authors documented (`quadruple_stability.outside_range`, or a hierarchy the network was not trained on) fall
  back to both triples passing at `fallback_margin` and the verdict is flagged `confident=False`.
- **A 2+2 quadruple** is two pairs about each other. Its two nested triples (one pair with the other pair as the third body) must
  pass at `level_margin`, then the 2+2 network of the same paper (`quadruple_stability`, GEN.1031) must call it stable with both
  `a_in / a_out` loosened by its `axis_margin` (1.5, tuned to the paper's precision on our own 2+2 REBOUND samples, GEN.1155; L2 in the
  registry); outside its documented range, or with overlapping apsides, the fallback is the 3+1's: both triples at `fallback_margin`,
  `confident=False`.

All of it is a pure function of the tree's stored orbits, so the generator can ask while it draws (`quadruple_verdict`
takes the numbers) and a test can ask of a stored system (`verdicts` takes the `OrbitNode`s and `TreeMember`s). The verifier
(`physics/clearance.py`) does not call it: it checks the generator with other criteria on purpose.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, stability
from planetgen.physics import quadruple_stability as quad


@dataclass(frozen=True)
class Verdict:
    """One level's answer. `shape` is `"triple"`, `"3+1"` or `"2+2"`; `node_no` the outermost node of the level; `criterion` what
    decided it (`"eq4"`, `"ma01"`, `"ek95"`, `"network"`, `"nested"`); `confident` false when the answer is an extrapolation or a
    stand-in for a criterion that is not built; `detail` says why in a few words."""
    node_no: int
    shape: str
    stable: bool
    criterion: str
    confident: bool
    detail: str = ""


def _mutual(inner, outer):
    return stability.mutual_inclination(math.radians(inner.inclination_deg), math.radians(inner.ascending_node_deg),
                                        math.radians(outer.inclination_deg), math.radians(outer.ascending_node_deg))


def triple_verdict(m1, m2, m3, a_in, a_out, e_in, e_out, mutual, margin=None):
    """Whether the triple clears `margin` (default `level_margin`) times its criterion's axis ratio; masses in solar masses, axes in AU,
    `mutual` in radians."""
    margin = quad.number("level_margin") if margin is None else margin
    method = stability.criterion(m1, m2, m3, a_in / a_out)
    required = stability.required_axis_ratio(m1, m2, m3, e_in, e_out, mutual, method, alpha=a_in / a_out)
    ratio = a_out / a_in
    return ratio > required * margin, method, ratio, required


def quadruple_verdict(masses, axes_au, eccentricities, inclinations, axis_margin=None, fallback_margin=None, level_margin=None):
    """
    The verdict of a 3+1 quadruple from numbers: `masses` `(m1, m2, m3, m4)` in solar masses (`m1 >= m2` the inner pair, `m3` the middle
    star, `m4` the outer one), `axes_au` `(a_in, a_mid, a_out)`, `eccentricities` `(e_in, e_mid, e_out)` and `inclinations`
    `(i_in_mid, i_in_out, i_mid_out)` in radians. Returns `(stable, criterion, confident, detail)`.
    """
    axis_margin = quad.number("axis_margin") if axis_margin is None else axis_margin
    fallback_margin = quad.number("fallback_margin") if fallback_margin is None else fallback_margin
    level_margin = quad.number("level_margin") if level_margin is None else level_margin
    m1, m2, m3, m4 = masses
    a_in, a_mid, a_out = axes_au
    e_in, e_mid, e_out = eccentricities
    i_in_mid, i_in_out, i_mid_out = inclinations
    if not a_in < a_mid < a_out:
        return False, "nested", True, "the three axes do not increase outward"
    inner_ok, _, _, _ = triple_verdict(m1, m2, m3, a_in, a_mid, e_in, e_mid, i_in_mid, level_margin)
    outer_ok, _, _, _ = triple_verdict(m1 + m2, m3, m4, a_mid, a_out, e_mid, e_out, i_mid_out, level_margin)
    if not (inner_ok and outer_ok):
        return False, "nested", True, "an inner or an outer triple fails its own criterion"
    ratios = quad.ratios(m1, m2, m3, m4, a_in, a_mid, a_out)
    beyond = quad.outside_range(*ratios, e_in, e_mid, e_out)
    if beyond or not quad.hierarchy_holds(ratios[3], ratios[4], e_in, e_mid, e_out):
        loose_in = triple_verdict(m1, m2, m3, a_in, a_mid, e_in, e_mid, i_in_mid, fallback_margin)[0]
        loose_out = triple_verdict(m1 + m2, m3, m4, a_mid, a_out, e_mid, e_out, i_mid_out, fallback_margin)[0]
        why = "outside the network's range (" + ", ".join(beyond) + ")" if beyond else "the apsides overlap, which the network was not trained on"
        return loose_in and loose_out, "nested", False, why
    q_in, q_mid, q_out, alpha_in_mid, alpha_mid_out = ratios
    stable = bool(quad.is_stable(q_in, q_mid, q_out, alpha_in_mid * axis_margin, alpha_mid_out * axis_margin,
                                 e_in, e_mid, e_out, i_in_mid, i_in_out, i_mid_out))
    return stable, "network", True, f"axis ratios loosened {axis_margin:g} times"


def quadruple_verdict_2p2(masses, axes_au, eccentricities, inclinations, axis_margin=None, fallback_margin=None, level_margin=None):
    """
    The verdict of a 2+2 quadruple from numbers: `masses` `(m1, m2, m3, m4)` in solar masses (`m1 >= m2` and `m3 >= m4` the two binaries, `m1 + m2 >=
    m3 + m4`), `axes_au` `(a_in1, a_in2, a_out)`, `eccentricities` `(e_in1, e_in2, e_out)` and `inclinations` `(i_in1_in2, i_in1_out, i_in2_out)`
    in radians. Returns `(stable, criterion, confident, detail)` as `quadruple_verdict` does.
    """
    axis_margin = quad.number_2p2("axis_margin") if axis_margin is None else axis_margin
    fallback_margin = quad.number("fallback_margin") if fallback_margin is None else fallback_margin
    level_margin = quad.number("level_margin") if level_margin is None else level_margin
    m1, m2, m3, m4 = masses
    a_in1, a_in2, a_out = axes_au
    e_in1, e_in2, e_out = eccentricities
    i_in1_in2, i_in1_out, i_in2_out = inclinations
    if not (a_in1 < a_out and a_in2 < a_out):
        return False, "nested", True, "an inner axis is not inside the outer one"
    first = triple_verdict(m1, m2, m3 + m4, a_in1, a_out, e_in1, e_out, i_in1_out, level_margin)[0]
    second = triple_verdict(m3, m4, m1 + m2, a_in2, a_out, e_in2, e_out, i_in2_out, level_margin)[0]
    if not (first and second):
        return False, "nested", True, "a binary with the other pair as its third body fails its own criterion"
    ratios = quad.ratios_2p2(m1, m2, m3, m4, a_in1, a_in2, a_out)
    beyond = quad.outside_range_2p2(*ratios, e_in1, e_in2, e_out)
    if beyond or not quad.hierarchy_holds_2p2(ratios[3], ratios[4], e_in1, e_in2, e_out):
        loose_first = triple_verdict(m1, m2, m3 + m4, a_in1, a_out, e_in1, e_out, i_in1_out, fallback_margin)[0]
        loose_second = triple_verdict(m3, m4, m1 + m2, a_in2, a_out, e_in2, e_out, i_in2_out, fallback_margin)[0]
        why = "outside the network's range (" + ", ".join(beyond) + ")" if beyond else "the apsides overlap, which the network was not trained on"
        return loose_first and loose_second, "nested", False, why
    q_in1, q_in2, q_out, alpha_in1_out, alpha_in2_out = ratios
    stable = bool(quad.is_stable_2p2(q_in1, q_in2, q_out, alpha_in1_out * axis_margin, alpha_in2_out * axis_margin,
                                     e_in1, e_in2, e_out, i_in1_in2, i_in1_out, i_in2_out))
    return stable, "network", True, f"axis ratios loosened {axis_margin:g} times"


class _Tree:
    """The children of every node: `kids[node_no] = {slot: ("member", member) | ("node", node)}`."""

    def __init__(self, nodes, members):
        self.nodes = {node.node_no: node for node in nodes}
        self.members = {member.member_no: member for member in members}
        self.kids = {}
        for member in self.members.values():
            self.kids.setdefault(member.parent_node_no, {})[member.slot] = ("member", member)
        for node in self.nodes.values():
            if node.parent_node_no is not None:
                self.kids.setdefault(node.parent_node_no, {})[node.slot] = ("node", node)

    def children(self, node_no):
        found = self.kids.get(node_no, {})
        return found.get(0), found.get(1)

    def mass(self, child):
        return child[1].mass_solar

    def pair(self, node_no):
        """The two children of a node that holds two members, as masses `(heavy, light)`, or None."""
        low, high = self.children(node_no)
        if low is None or high is None or low[0] != "member" or high[0] != "member":
            return None
        return tuple(sorted((low[1].mass_solar, high[1].mass_solar), reverse=True))


def _axis_au(node):
    return node.separation_km / constants.AU_TO_KM


def _levels(tree):
    """Every `(outer node, inner node, third child)` where a child of the outer node is a node with two children."""
    out = []
    for outer in tree.nodes.values():
        low, high = tree.children(outer.node_no)
        if low is None or high is None:
            continue
        for inner_side, third in ((low, high), (high, low)):
            if inner_side[0] == "node":
                out.append((outer, inner_side[1], third))
    return out


def _triple(tree, outer, inner, third, margin):
    first, second = tree.children(inner.node_no)
    if first is None or second is None:
        return None
    m1, m2 = sorted((tree.mass(first), tree.mass(second)), reverse=True)
    stable, method, ratio, required = triple_verdict(m1, m2, tree.mass(third), _axis_au(inner), _axis_au(outer), inner.eccentricity,
                                                     outer.eccentricity, _mutual(inner, outer), margin)
    detail = f"a_out/a_in {ratio:.4g} against {required:.4g} by {method}" + ("" if margin == 1.0 else f" times {margin:g}")
    return Verdict(outer.node_no, "triple", stable, method, True, detail)


def _chain_3p1(tree, outer):
    """`(pair node, middle node, outer node)` where `outer` holds a star and a node, that node a star and a pair of stars; else None."""
    for outer_side in tree.children(outer.node_no):
        if outer_side is None or outer_side[0] != "node":
            continue
        middle = outer_side[1]
        for middle_side in tree.children(middle.node_no):
            if middle_side is None or middle_side[0] != "node":
                continue
            inner = middle_side[1]
            if tree.pair(inner.node_no) is None:
                continue
            others = [side for side in tree.children(outer.node_no) if side is not outer_side] + [
                side for side in tree.children(middle.node_no) if side is not middle_side]
            if len(others) == 2 and all(side is not None and side[0] == "member" for side in others):
                return inner, middle, outer
    return None


def verdicts(nodes, members, level_margin=None, axis_margin=None, fallback_margin=None):
    """
    Every level's verdict of a stored tree (`OrbitNode`s and `TreeMember`s of `physics/system_tree`): one `"triple"` for each node with a
    sub-node child, a `"3+1"` for each inner pair, middle star and outer star, and a `"2+2"` for a node of two pairs. A `group` node is
    not hierarchical and has no verdict. Empty for one or two stars.
    """
    level_margin = quad.number("level_margin") if level_margin is None else level_margin
    tree = _Tree([node for node in nodes if node.kind != "group"], members)
    found = []
    for outer, inner, third in _levels(tree):
        verdict = _triple(tree, outer, inner, third, level_margin)
        if verdict is not None:
            found.append(verdict)
    for outer in tree.nodes.values():
        chain = _chain_3p1(tree, outer)
        if chain is not None:
            found.append(_verdict_3p1(tree, chain, level_margin, axis_margin, fallback_margin))
        low, high = tree.children(outer.node_no)
        if low is not None and high is not None and low[0] == "node" and high[0] == "node" \
                and tree.pair(low[1].node_no) and tree.pair(high[1].node_no):
            found.append(_verdict_2p2(tree, outer, low[1], high[1], level_margin, axis_margin, fallback_margin))
    return found


def _verdict_3p1(tree, chain, level_margin, axis_margin, fallback_margin):
    inner, middle, outer = chain
    m1, m2 = tree.pair(inner.node_no)
    star3 = next(side for side in tree.children(middle.node_no) if side[0] == "member")
    star4 = next(side for side in tree.children(outer.node_no) if side[0] == "member")
    stable, criterion, confident, detail = quadruple_verdict(
        (m1, m2, tree.mass(star3), tree.mass(star4)), (_axis_au(inner), _axis_au(middle), _axis_au(outer)),
        (inner.eccentricity, middle.eccentricity, outer.eccentricity),
        (_mutual(inner, middle), _mutual(inner, outer), _mutual(middle, outer)), axis_margin, fallback_margin, level_margin)
    return Verdict(outer.node_no, "3+1", stable, criterion, confident, detail)


def _verdict_2p2(tree, outer, left, right, level_margin, axis_margin, fallback_margin):
    heavy, light = sorted((left, right), key=lambda pair: -pair.mass_solar)
    m1, m2 = tree.pair(heavy.node_no)
    m3, m4 = tree.pair(light.node_no)
    stable, criterion, confident, detail = quadruple_verdict_2p2(
        (m1, m2, m3, m4), (_axis_au(heavy), _axis_au(light), _axis_au(outer)),
        (heavy.eccentricity, light.eccentricity, outer.eccentricity),
        (_mutual(heavy, light), _mutual(heavy, outer), _mutual(light, outer)), axis_margin, fallback_margin, level_margin)
    return Verdict(outer.node_no, "2+2", stable, criterion, confident, detail)


def tree_is_stable(nodes, members, **margins):
    """Whether every level of the tree passes (True for one or two stars)."""
    return all(verdict.stable for verdict in verdicts(nodes, members, **margins))
