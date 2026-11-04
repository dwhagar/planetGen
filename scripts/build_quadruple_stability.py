#!/usr/bin/env python3
"""
Builds `src/planetgen/physics/data/quadruple_stability_3p1_weights.json` (GEN.1029) and `quadruple_stability_2p2_weights.json` (GEN.1031): the weights
of the 3+1 and 2+2 quadruple stability classifiers of [R4575] as plain numbers, so the
generator needs no machine-learning library.

    python scripts/build_quadruple_stability.py 3p1 /path/to/mlp_model_3p1_ghost_v1.2.2.pkl
    python scripts/build_quadruple_stability.py 2p2 /path/to/mlp_model_2p2_ghost_v1.2.2.pkl

The pickles are https://github.com/pavanvyn/quadruple-stability/blob/main/mlp_model_3p1_ghost_v1.2.2.pkl and the same name with 2p2 (209 KB each, never committed). A pickle can run code
when it is loaded, so it is read with an unpickler that refuses every global except the array constructors of numpy and replaces every scikit-learn
class with a bag that only keeps its state; nothing else is executed. The state of a scikit-learn `MLPClassifier` holds the weight matrices
(`coefs_`), the biases (`intercepts_`) and the activation names; the script refuses a model that is not 11 inputs, four hidden layers of 50 and one
output, all logistic. Numbers are written to nine significant digits (a relative change under 1e-9, a change of the probability under 1e-8).

The validity ranges and the tier and tweak level of every number are in `quadruple_stability_3p1.json` and `quadruple_stability_2p2.json`, which this script
does not touch. Both networks have the same shape (11 inputs, four hidden layers of 50, one output); only the meaning of the inputs differs.
"""

import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "src", "planetgen", "physics", "data")
KINDS = ("3p1", "2p2")
LAYERS = (11, 50, 50, 50, 50, 1)


class _Bag:
    """Stands in for any scikit-learn class: keeps the state it is given and does nothing else."""

    def __init__(self, *args, **kwargs):
        pass

    def __setstate__(self, state):
        self.state = state

    def __call__(self, *args, **kwargs):
        return _Bag()


_CORE = getattr(np, "_core", None) or np.core
_ALLOWED = {
    ("numpy", "dtype"): np.dtype,
    ("numpy", "ndarray"): np.ndarray,
    ("numpy.core.multiarray", "_reconstruct"): _CORE.multiarray._reconstruct,
    ("numpy.core.multiarray", "scalar"): _CORE.multiarray.scalar,
    ("numpy._core.multiarray", "_reconstruct"): _CORE.multiarray._reconstruct,
    ("numpy._core.multiarray", "scalar"): _CORE.multiarray.scalar,
}


class _SafeUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) in _ALLOWED:
            return _ALLOWED[(module, name)]
        if module.startswith("sklearn") or module.startswith("numpy.random"):
            return _Bag
        raise pickle.UnpicklingError(f"refused {module}.{name}")


def read_state(path):
    """The state dictionary of the scikit-learn model in the pickle at `path`, read without running its code."""
    with open(path, "rb") as handle:
        return _SafeUnpickler(handle).load().state


def network(state):
    """`{"activation", "coefs", "intercepts"}` from a model state; refuses any other shape of network."""
    coefs = [np.asarray(c, dtype=float) for c in state["coefs_"]]
    intercepts = [np.asarray(b, dtype=float) for b in state["intercepts_"]]
    shapes = tuple([coefs[0].shape[0]] + [c.shape[1] for c in coefs])
    if shapes != LAYERS or state["activation"] != "logistic" or state["out_activation_"] != "logistic" or list(state["classes_"]) != [0, 1]:
        raise SystemExit(f"expected a logistic network of layers {LAYERS} with classes [0, 1], got {shapes}, {state['activation']}, {state['out_activation_']}")
    rounded = lambda a: [[float(f"{x:.9g}") for x in row] for row in np.atleast_2d(a)] if a.ndim == 2 else [float(f"{x:.9g}") for x in a]
    return {"activation": "logistic", "layers": list(LAYERS), "coefs": [rounded(c) for c in coefs], "intercepts": [rounded(b) for b in intercepts]}


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in KINDS:
        sys.exit(f"usage: build_quadruple_stability.py {'|'.join(KINDS)} /path/to/model.pkl")
    output = os.path.join(DATA, f"quadruple_stability_{sys.argv[1]}_weights.json")
    with open(output, "w", encoding="utf-8") as handle:
        json.dump(network(read_state(sys.argv[2])), handle, separators=(",", ":"))
        handle.write("\n")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
