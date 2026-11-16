# tests/test_pair_overlap.py

"""
OPS.242: the pairs of commands that may overlap (`queue/independence.json`), each proved by running them overlapped (`tests/pair_overlap_support.py`, `tests/pair_overlap_recipes.py`).

* `test_pair_is_independent` is the proof: the two stages overlapped, in both orders, at 1 and 15 workers, store the same content as one after the other, and each command touches only the
  tables `queue/touches.py` declares for it. It fails with the rows that differ.
* The file and the recipes are held together: every listed pair has a recipe, every recipe is listed, and the tables the pair is declared to share (`touches.shared_between`) are the ones
  its entry names, so a declaration that changes (a new conflict) fails here until the pair is proved again.
"""

import pytest

from planetgen.queue import touches
from tests import pair_overlap_recipes as recipes
from tests.cut_neutral_support import Context
from tests.pair_overlap_support import prove


def _listed():
    return {tuple(sorted(pair)): entry for pair, entry in touches.load_pairs().items()}


def _key(pair):
    return tuple(sorted(pair))


def test_every_listed_pair_has_a_recipe_and_every_recipe_is_listed():
    assert set(_listed()) == {_key(pair) for pair in recipes.PAIRS}


@pytest.mark.parametrize("pair", sorted(recipes.PAIRS))
def test_a_listed_pair_names_what_it_shares_and_its_proof(pair):
    entry = _listed()[_key(pair)]
    assert entry["grain"] == "table"
    assert entry["shared"] == touches.shared_between(*pair), f"{pair}: the tables the commands share changed; prove the pair again and list what it shares"
    assert entry["proof"] == f"tests/test_pair_overlap.py::test_pair_is_independent[{pair[0]}-{pair[1]}]"


@pytest.mark.db
@pytest.mark.parametrize("pair", sorted(recipes.PAIRS), ids=lambda pair: f"{pair[0]}-{pair[1]}")
def test_pair_is_independent(mysql_config, monkeypatch, redis_server, pair):
    verdict = prove(recipes.PAIRS[pair], Context(mysql_config, monkeypatch))
    print(verdict.line(), "|", verdict.stored, flush=True)
    assert verdict.passed, verdict.line() + "".join(f"\n    {line[:400]}" for line in [*verdict.undeclared, *verdict.differences])


@pytest.mark.db
def test_the_proof_fails_for_two_stages_that_are_not_independent(mysql_config, monkeypatch, redis_server):
    """The luminosity pass reads the marks the mass pass leaves when it is cut; overlapped, it is cut before them and stores other stars. The proof must see that."""
    verdict = prove(recipes.CONTROL, Context(mysql_config, monkeypatch))
    print(verdict.line(), flush=True)
    assert not verdict.passed and verdict.differences
