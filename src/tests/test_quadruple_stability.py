# tests/test_quadruple_stability.py

"""GEN.1029 and GEN.1031: the 3+1 and 2+2 quadruple stability networks of Vynatheya, Mardling and Hamers 2023, against the scikit-learn models they were made from."""

import os
import pathlib
import pickle
import sys

import numpy as np
import pytest

from planetgen.physics import quadruple_stability as quad

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_quadruple_stability as builder  # noqa: E402

# Probabilities of 'unstable' from the scikit-learn model mlp_model_3p1_ghost_v1.2.2.pkl itself (predict_proba), read 10-30.
SKLEARN = (
    ((1.0, 0.5, 0.33, 0.2, 0.2, 0, 0, 0, 0, 0, 0), 0.0009255504192062601),
    ((0.5, 1.2, 0.8, 0.1, 0.15, 0.3, 0.2, 0.1, 0.5, 1.0, 0.7), 0.017691511417149492),
    ((0.3, 2.5, 1.5, 0.35, 0.3, 0.6, 0.5, 0.4, 1.5, 2.0, 2.5), 0.9999999706663643),
    ((0.9, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1, 0.2, 0.3, 0.1), 1.042854800049108e-06),
    ((0.2, 0.3, 0.2, 0.25, 0.4, 0.0, 0.7, 0.5, 0.0, 0.0, 0.0), 0.9999995836619994),
)


# The same for mlp_model_2p2_ghost_v1.2.2.pkl: an MLPClassifier given the state of the pickle (read without running its code) and asked for predict_proba, 10-30.
SKLEARN_2P2 = (
    ((1.0, 1.0, 1.0, 0.2, 0.2, 0, 0, 0, 0, 0, 0), 0.004994214692352674),
    ((0.5, 0.8, 0.7, 0.1, 0.15, 0.3, 0.2, 0.1, 0.5, 1.0, 0.7), 0.0006070406302067851),
    ((0.3, 0.4, 0.9, 0.35, 0.3, 0.6, 0.5, 0.4, 1.5, 2.0, 2.5), 0.9998603646350366),
    ((0.9, 0.9, 0.5, 0.05, 0.05, 0.1, 0.1, 0.1, 0.2, 0.3, 0.1), 6.146184382850036e-10),
    ((0.2, 0.5, 0.6, 0.25, 0.4, 0.0, 0.7, 0.5, 0.0, 0.0, 0.0), 0.9999812592135476),
)


@pytest.mark.parametrize("inputs, probability", SKLEARN)
def test_the_numpy_forward_pass_gives_the_scikit_learn_probability(inputs, probability):
    assert quad.unstable_probability(*inputs) == pytest.approx(probability, abs=1e-6)


def test_the_weights_are_the_authors_network_of_8301_numbers():
    sizes = [c.size for c in quad._COEFS] + [b.size for b in quad._INTERCEPTS]
    assert sum(sizes) == 8301
    assert [c.shape for c in quad._COEFS] == [(11, 50), (50, 50), (50, 50), (50, 50), (50, 1)]
    assert quad.INPUT_NAMES[:3] == ("q_in", "q_mid", "q_out") and len(quad.INPUT_NAMES) == 11


@pytest.mark.parametrize("inputs, probability", SKLEARN_2P2)
def test_the_2p2_forward_pass_gives_the_scikit_learn_probability(inputs, probability):
    assert quad.unstable_probability_2p2(*inputs) == pytest.approx(probability, abs=1e-6)


def test_the_2p2_weights_are_the_authors_network_of_8301_numbers():
    sizes = [c.size for c in quad._COEFS_2P2] + [b.size for b in quad._INTERCEPTS_2P2]
    assert sum(sizes) == 8301
    assert [c.shape for c in quad._COEFS_2P2] == [(11, 50), (50, 50), (50, 50), (50, 50), (50, 1)]
    assert quad.INPUT_NAMES_2P2[:3] == ("q_in1", "q_in2", "q_out") and len(quad.INPUT_NAMES_2P2) == 11
    assert not np.allclose(quad._COEFS[0], quad._COEFS_2P2[0])


def test_the_authors_example_system_is_stable_and_a_loose_one_is_not():
    """The command line example of the authors' README (q 1, 0.5, 0.33; alpha 0.2, 0.2; circular, coplanar) prints 'ML stable'."""
    assert quad.is_stable(1.0, 0.5, 0.33, 0.2, 0.2, 0.0, 0.0, 0.0)
    assert quad.is_stable(1.0, 0.5, 0.33, 0.01, 0.01, 0.0, 0.0, 0.0)
    assert not quad.is_stable(1.0, 0.5, 0.33, 0.4, 0.4, 0.0, 0.0, 0.0)


def test_arrays_broadcast_and_a_tighter_system_is_never_more_stable_along_the_axis_ratio():
    alpha = np.linspace(0.02, 0.45, 40)
    stable = quad.is_stable(0.8, 0.6, 0.4, alpha, alpha, 0.1, 0.1, 0.1)
    assert stable.shape == alpha.shape and stable[0] and not stable[-1]
    assert np.all(np.diff(stable.astype(int)) <= 0)


def test_the_2p2_authors_example_is_stable_and_a_loose_one_is_not():
    """The command line example of the authors' README for 2+2 (q 1, 1, 1; alpha 0.2, 0.2; circular, coplanar)."""
    assert quad.is_stable_2p2(1.0, 1.0, 1.0, 0.2, 0.2, 0.0, 0.0, 0.0)
    assert quad.is_stable_2p2(1.0, 1.0, 1.0, 0.01, 0.01, 0.0, 0.0, 0.0)
    assert not quad.is_stable_2p2(1.0, 1.0, 1.0, 0.4, 0.4, 0.0, 0.0, 0.0)


def test_2p2_arrays_broadcast_and_a_tighter_system_is_never_more_stable_along_the_axis_ratio():
    alpha = np.linspace(0.02, 0.45, 40)
    stable = quad.is_stable_2p2(0.8, 0.6, 0.7, alpha, alpha, 0.1, 0.1, 0.1)
    assert stable.shape == alpha.shape and stable[0] and not stable[-1]
    assert np.all(np.diff(stable.astype(int)) <= 0)


def test_probabilities_stay_between_zero_and_one():
    rng = np.random.default_rng(3)
    p = quad.unstable_probability(rng.uniform(0.1, 1, 500), rng.uniform(0.05, 3, 500), rng.uniform(0.05, 2, 500), rng.uniform(0.01, 0.5, 500),
                                  rng.uniform(0.01, 0.5, 500), rng.uniform(0, 0.95, 500), rng.uniform(0, 0.95, 500), rng.uniform(0, 0.95, 500),
                                  rng.uniform(0, np.pi, 500), rng.uniform(0, np.pi, 500), rng.uniform(0, np.pi, 500))
    assert np.all((p >= 0) & (p <= 1))


def test_the_ratios_of_masses_and_axes():
    assert quad.ratios(2.0, 1.0, 1.5, 2.25, 0.1, 0.5, 5.0) == pytest.approx((0.5, 0.5, 0.5, 0.2, 0.1))
    with pytest.raises(ValueError):
        quad.ratios(1.0, 2.0, 1.0, 1.0, 0.1, 0.5, 5.0)
    with pytest.raises(ValueError):
        quad.ratios(2.0, 1.0, 1.0, 1.0, 0.5, 0.1, 5.0)


def test_the_2p2_ratios_of_masses_and_axes():
    assert quad.ratios_2p2(2.0, 1.0, 1.5, 0.75, 0.1, 0.2, 1.0) == pytest.approx((0.5, 0.5, 0.75, 0.1, 0.2))
    with pytest.raises(ValueError):
        quad.ratios_2p2(1.0, 2.0, 1.0, 0.5, 0.1, 0.2, 1.0)
    with pytest.raises(ValueError):
        quad.ratios_2p2(1.0, 0.5, 2.0, 1.0, 0.1, 0.2, 1.0)
    with pytest.raises(ValueError):
        quad.ratios_2p2(2.0, 1.0, 1.0, 0.5, 0.1, 1.2, 1.0)


def test_the_2p2_range_check_and_hierarchy():
    assert quad.outside_range_2p2(0.5, 0.5, 0.8, 0.2, 0.2, 0.1, 0.1, 0.1) == ()
    assert quad.outside_range_2p2(0.05, 1.2, 1.5, 0.2, 1.2, 0.1, 0.99, 0.1) == ("q_in1", "q_in2", "q_out", "alpha_in2_out", "e_in2")
    assert quad.hierarchy_holds_2p2(0.2, 0.2, 0.3, 0.3, 0.3)
    assert not quad.hierarchy_holds_2p2(0.6, 0.2, 0.5, 0.0, 0.5)
    assert not quad.hierarchy_holds_2p2(0.2, 0.6, 0.0, 0.5, 0.5)


def test_the_range_check_names_what_the_authors_did_not_document():
    assert quad.outside_range(0.5, 1.0, 1.0, 0.2, 0.2, 0.1, 0.1, 0.1) == ()
    assert quad.outside_range(0.05, 4.0, 2.5, 0.2, 1.2, 0.1, 0.99, 0.1) == ("q_in", "q_mid", "q_out", "alpha_mid_out", "e_mid")


def test_the_hierarchy_test_is_apoapsis_inside_the_next_periapsis():
    assert quad.hierarchy_holds(0.2, 0.2, 0.3, 0.3, 0.3)
    assert not quad.hierarchy_holds(0.6, 0.2, 0.5, 0.5, 0.0)
    assert not quad.hierarchy_holds(0.2, 0.6, 0.0, 0.5, 0.5)


def test_every_number_in_the_rule_file_carries_a_tier_a_tweak_level_and_a_basis():
    for name, entry in {**quad.NUMBERS, **{"2p2 " + k: v for k, v in quad.NUMBERS_2P2.items()}}.items():
        assert entry["tier"] in {f"T{n}" for n in range(1, 10)} and entry["tweak"] in {f"L{n}" for n in range(5)} and entry["basis"], name


class _Boom:
    def __reduce__(self):
        return (os.system, ("echo should-not-run > /dev/null",))


def test_the_builder_reads_a_pickle_without_running_its_code(tmp_path):
    path = tmp_path / "evil.pkl"
    path.write_bytes(pickle.dumps(_Boom()))
    with pytest.raises(pickle.UnpicklingError):
        builder.read_state(str(path))


def test_the_builder_refuses_a_network_of_another_shape():
    state = {"coefs_": [np.zeros((11, 50)), np.zeros((50, 1))], "intercepts_": [np.zeros(50), np.zeros(1)], "activation": "logistic",
             "out_activation_": "logistic", "classes_": np.array([0, 1])}
    with pytest.raises(SystemExit):
        builder.network(state)
