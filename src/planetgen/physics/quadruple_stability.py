# planetgen/physics/quadruple_stability.py

"""
Stability of a 3+1 or a 2+2 quadruple star system (GEN.1029 and GEN.1031, docs/design/systems/quadruple-stability-3p1.md; the criterion GEN.1016 puts in the multi-star limits)
=================================================================================================================================================================================

A 3+1 quadruple is an inner binary (`m1 >= m2`), a third star `m3` on the intermediate orbit around the binary, and a fourth star `m4` on the
outer orbit around those three. Testing it as two nested triples (`physics/stability.py`) is right for 66 percent of such systems and calls
too many of them stable (precision 0.59), because the outer star excites the eccentricity of the intermediate orbit and that destabilises the
inner triple, which a triple test cannot see. Vynatheya, Mardling and Hamers 2023 (arXiv 2301.09930) trained a small neural network on 5e5
N-body integrations and got 93 percent. This module is that network as plain numbers: 11 inputs, four hidden layers of 50 logistic units, one
logistic output, the probability that the system is *unstable*. The 8,301 numbers (weights and biases) are in `data/quadruple_stability_3p1_weights.json`, made from the
authors' scikit-learn model file by `scripts/build_quadruple_stability.py` (read without running the pickle); the forward pass here is numpy
only and gives the scikit-learn probabilities to 1e-8.

The inputs, all dimensionless (the problem is scale free):

    q_in = m2/m1,  q_mid = m3/(m1+m2),  q_out = m4/(m1+m2+m3)
    alpha_in_mid = a_in/a_mid,  alpha_mid_out = a_mid/a_out
    e_in, e_mid, e_out
    i_in_mid, i_in_out, i_mid_out    (mutual inclinations of the orbit pairs, radians; the network gets them over pi)

A 2+2 quadruple is two binaries (`m1 >= m2` and `m3 >= m4`, with `m1 + m2 >= m3 + m4`) that orbit each other. The paper's second network is the same
shape (the same 11 inputs, four hidden layers of 50 logistic units, 8,301 numbers, in `data/quadruple_stability_2p2_weights.json`); its inputs are

    q_in1 = m2/m1,  q_in2 = m4/m3,  q_out = (m3+m4)/(m1+m2)
    alpha_in1_out = a_in1/a_out,  alpha_in2_out = a_in2/a_out
    e_in1, e_in2, e_out
    i_in1_in2, i_in1_out, i_in2_out

and the functions for it carry the suffix `_2p2` (`ratios_2p2`, `unstable_probability_2p2`, `is_stable_2p2`, `hierarchy_holds_2p2`, `outside_range_2p2`). The authors' scores
for 2+2 are 94 percent overall, precision 0.94 for stable and 0.95 for unstable, recall 0.94 and 0.95 (their table 3) against 83 percent for two nested
Mardling-Aarseth triples. The paper explains why the nested test does better for a 2+2: it treats each inner binary as a point mass, which is
consistent when both binaries are tight, so the network's gain over the nested test is smaller here.

What "stable" means here: the system stays bound for 100 outer orbits and a twin with its inner axis larger by 1e-6 au does not diverge from
it by more than 1 percent in the inner axis (`outer_orbits`, `divergence_limit`). It is a chaos test over a short span, not a guarantee for a
stellar lifetime, and secular evolution on longer times is not in it.

Limits (what the number is not): the network was trained on masses within a factor 10 of each other, eccentricities to 0.95, isotropic orbital
angles and systems near the Mardling-Aarseth boundary of the two nested triples; `outside_range` names the inputs outside what is
documented, and the verdict there is extrapolation. The arguments of periapsis and the mean anomalies are not inputs. The authors' scores:
93 percent overall, precision 0.91 for stable and 0.95 for unstable, recall 0.91 and 0.95 (their table 4); the check in the design note on an independent
sample is in its section 3.
"""

import json
import os

import numpy as np

from planetgen.physics import elements
from planetgen.util import detmath

INPUT_NAMES = ("q_in", "q_mid", "q_out", "alpha_in_mid", "alpha_mid_out", "e_in", "e_mid", "e_out", "i_in_mid", "i_in_out", "i_mid_out")
"""tuple: The eleven inputs of the 3+1 network, in the order it takes them."""

INPUT_NAMES_2P2 = ("q_in1", "q_in2", "q_out", "alpha_in1_out", "alpha_in2_out", "e_in1", "e_in2", "e_out", "i_in1_in2", "i_in1_out", "i_in2_out")
"""tuple: The eleven inputs of the 2+2 network, in the order it takes them."""


class QuadrupleStabilityError(ValueError):
    """A rule file that is refused."""


def _load(kind, input_names):
    """`(rule document, coefficient matrices, bias vectors)` of the network of `kind` (`"3p1"` or `"2p2"`), refusing a file that does not hold what it claims."""
    data_file = os.path.join(elements.DATA_DIRECTORY, f"quadruple_stability_{kind}.json")
    weights_file = os.path.join(elements.DATA_DIRECTORY, f"quadruple_stability_{kind}_weights.json")
    with open(data_file, encoding="utf-8") as handle:
        document = json.load(handle)
    for name, entry in document["numbers"].items():
        if entry.get("tier") not in tuple(f"T{n}" for n in range(1, 10)) or entry.get("tweak") not in tuple(f"L{n}" for n in range(5)) \
                or not entry.get("basis") or "value" not in entry:
            raise QuadrupleStabilityError(f"{name}: every number needs a value, a tier T1..T9, a tweak L0..L4 and a basis, got {entry!r}")
    with open(weights_file, encoding="utf-8") as handle:
        network = json.load(handle)
    coefs = [np.array(c, dtype=float) for c in network["coefs"]]
    intercepts = [np.array(b, dtype=float) for b in network["intercepts"]]
    layers = network["layers"]
    if network["activation"] != "logistic" or layers[0] != len(input_names) or layers[-1] != 1 \
            or [c.shape for c in coefs] != [(layers[i], layers[i + 1]) for i in range(len(layers) - 1)] \
            or [b.shape for b in intercepts] != [(n,) for n in layers[1:]]:
        raise QuadrupleStabilityError(f"{weights_file}: the network does not match its layers {layers}")
    return document, coefs, intercepts


DOCUMENT, _COEFS, _INTERCEPTS = _load("3p1", INPUT_NAMES)
NUMBERS = DOCUMENT["numbers"]
VERSION = DOCUMENT["version"]
DOCUMENT_2P2, _COEFS_2P2, _INTERCEPTS_2P2 = _load("2p2", INPUT_NAMES_2P2)
NUMBERS_2P2 = DOCUMENT_2P2["numbers"]


def number(name):
    """The value of the 3+1 rule file's number `name`."""
    return NUMBERS[name]["value"]


def number_2p2(name):
    """The value of the 2+2 rule file's number `name`."""
    return NUMBERS_2P2[name]["value"]


def _forward(coefs, intercepts, inputs):
    """The logistic network on `inputs` (eleven numbers or arrays that broadcast, inclinations already over pi): the probability of 'unstable'."""
    activation = np.stack(np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in inputs]), axis=-1)
    for weights, bias in zip(coefs, intercepts):
        activation = 1.0 / (1.0 + np.exp(-(detmath.matmul(activation, weights) + bias)))
    return activation[..., 0]


def ratios(m1, m2, m3, m4, a_in, a_mid, a_out):
    """`(q_in, q_mid, q_out, alpha_in_mid, alpha_mid_out)` of a 3+1 system from its masses (`m1 >= m2`) and its three semi-major axes (any one length unit)."""
    if not m1 >= m2 > 0 or not m3 > 0 or not m4 > 0:
        raise ValueError("masses must be positive and the inner pair ordered m1 >= m2")
    if not 0 < a_in < a_mid < a_out:
        raise ValueError("the axes must increase from the inner binary to the outer star")
    return m2 / m1, m3 / (m1 + m2), m4 / (m1 + m2 + m3), a_in / a_mid, a_mid / a_out


def unstable_probability(q_in, q_mid, q_out, alpha_in_mid, alpha_mid_out, e_in, e_mid, e_out, i_in_mid=0.0, i_in_out=0.0, i_mid_out=0.0):
    """
    The network's probability that the system is unstable, 0 to 1. Every argument is a number or an array (they broadcast); angles in radians.
    """
    return _forward(_COEFS, _INTERCEPTS, (q_in, q_mid, q_out, alpha_in_mid, alpha_mid_out, e_in, e_mid, e_out,
                                          np.asarray(i_in_mid) / np.pi, np.asarray(i_in_out) / np.pi, np.asarray(i_mid_out) / np.pi))


def is_stable(q_in, q_mid, q_out, alpha_in_mid, alpha_mid_out, e_in, e_mid, e_out, i_in_mid=0.0, i_in_out=0.0, i_mid_out=0.0):
    """Whether the network calls the system stable (its probability of 'unstable' is below `stable_below`); an array for array inputs."""
    return unstable_probability(q_in, q_mid, q_out, alpha_in_mid, alpha_mid_out, e_in, e_mid, e_out, i_in_mid, i_in_out, i_mid_out) < number("stable_below")


def hierarchy_holds(alpha_in_mid, alpha_mid_out, e_in, e_mid, e_out):
    """Whether each orbit's apoapsis is inside the next orbit's periapsis, the hierarchy the training sample kept."""
    return bool(alpha_in_mid * (1.0 + e_in) < 1.0 - e_mid and alpha_mid_out * (1.0 + e_mid) < 1.0 - e_out)


def outside_range(q_in, q_mid, q_out, alpha_in_mid, alpha_mid_out, e_in, e_mid, e_out):
    """The names of the inputs outside what the authors documented (an empty tuple inside); the verdict outside is extrapolation."""
    outside = []
    if q_in < number("q_in_min") or q_in > 1.0:
        outside.append("q_in")
    if q_mid > number("q_mid_max"):
        outside.append("q_mid")
    if q_out > number("q_out_max"):
        outside.append("q_out")
    if not 0.0 < alpha_in_mid < 1.0:
        outside.append("alpha_in_mid")
    if not 0.0 < alpha_mid_out < 1.0:
        outside.append("alpha_mid_out")
    for name, value in (("e_in", e_in), ("e_mid", e_mid), ("e_out", e_out)):
        if not 0.0 <= value <= number("eccentricity_max"):
            outside.append(name)
    return tuple(outside)


def ratios_2p2(m1, m2, m3, m4, a_in1, a_in2, a_out):
    """`(q_in1, q_in2, q_out, alpha_in1_out, alpha_in2_out)` of a 2+2 system from its masses (`m1 >= m2`, `m3 >= m4`, `m1 + m2 >= m3 + m4`) and its three semi-major axes."""
    if not m1 >= m2 > 0 or not m3 >= m4 > 0:
        raise ValueError("masses must be positive and each binary ordered, m1 >= m2 and m3 >= m4")
    if not m1 + m2 >= m3 + m4:
        raise ValueError("the first binary must be the heavier, m1 + m2 >= m3 + m4")
    if not 0 < a_in1 < a_out or not 0 < a_in2 < a_out:
        raise ValueError("both inner axes must be smaller than the outer axis")
    return m2 / m1, m4 / m3, (m3 + m4) / (m1 + m2), a_in1 / a_out, a_in2 / a_out


def unstable_probability_2p2(q_in1, q_in2, q_out, alpha_in1_out, alpha_in2_out, e_in1, e_in2, e_out, i_in1_in2=0.0, i_in1_out=0.0, i_in2_out=0.0):
    """The 2+2 network's probability that the system is unstable, 0 to 1. Every argument is a number or an array (they broadcast); angles in radians."""
    return _forward(_COEFS_2P2, _INTERCEPTS_2P2, (q_in1, q_in2, q_out, alpha_in1_out, alpha_in2_out, e_in1, e_in2, e_out,
                                                 np.asarray(i_in1_in2) / np.pi, np.asarray(i_in1_out) / np.pi, np.asarray(i_in2_out) / np.pi))


def is_stable_2p2(q_in1, q_in2, q_out, alpha_in1_out, alpha_in2_out, e_in1, e_in2, e_out, i_in1_in2=0.0, i_in1_out=0.0, i_in2_out=0.0):
    """Whether the 2+2 network calls the system stable (its probability of 'unstable' is below `stable_below`); an array for array inputs."""
    return unstable_probability_2p2(q_in1, q_in2, q_out, alpha_in1_out, alpha_in2_out, e_in1, e_in2, e_out, i_in1_in2, i_in1_out, i_in2_out) < number_2p2("stable_below")


def hierarchy_holds_2p2(alpha_in1_out, alpha_in2_out, e_in1, e_in2, e_out):
    """Whether each inner binary's apoapsis is inside the outer orbit's periapsis, the hierarchy the training sample kept."""
    return bool(alpha_in1_out * (1.0 + e_in1) < 1.0 - e_out and alpha_in2_out * (1.0 + e_in2) < 1.0 - e_out)


def outside_range_2p2(q_in1, q_in2, q_out, alpha_in1_out, alpha_in2_out, e_in1, e_in2, e_out):
    """The names of the 2+2 inputs outside what the authors documented (an empty tuple inside); the verdict outside is extrapolation."""
    outside = []
    for name, value in (("q_in1", q_in1), ("q_in2", q_in2)):
        if value < number_2p2("q_in_min") or value > 1.0:
            outside.append(name)
    if not 0.0 < q_out <= number_2p2("q_out_max"):
        outside.append("q_out")
    for name, value in (("alpha_in1_out", alpha_in1_out), ("alpha_in2_out", alpha_in2_out)):
        if not 0.0 < value < 1.0:
            outside.append(name)
    for name, value in (("e_in1", e_in1), ("e_in2", e_in2), ("e_out", e_out)):
        if not 0.0 <= value <= number_2p2("eccentricity_max"):
            outside.append(name)
    return tuple(outside)
