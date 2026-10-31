# planetgen/util/request_odds.py

"""
What a request does to a draw (GEN.1130)
=========================================

`+habitable_world` and `+intelligent_life` raise the odds of the wanted outcome and never force it. Each draw that leads to the outcome is
made as it would be without the request; the request then keeps a draw that gave the outcome always, and one that did not only one time in
`tuning.REQUEST_ODDS`, so the odds of the outcome (its chance over the chance it does not happen) are multiplied by `tuning.REQUEST_ODDS` and
no more. A draw that is not kept is made again, as many times as the draw's own limit allows, and the last one stands: the request may come
out unmet, and the system says so. Nothing is placed, aimed or adapted to make the outcome.
"""

from planetgen import tuning
from planetgen.util import draw


def keeps(wanted, rng=draw):
    """Whether a draw stands under a request: always when it gave the wanted outcome, else one time in `tuning.REQUEST_ODDS`. One draw, always."""
    pick = rng.random()
    return bool(wanted) or pick < 1.0 / tuning.REQUEST_ODDS


def raised_chance(chance):
    """The chance `chance` becomes when its odds are multiplied by `tuning.REQUEST_ODDS` (0 and 1 stay put)."""
    if chance <= 0.0 or chance >= 1.0:
        return min(max(chance, 0.0), 1.0)
    odds = tuning.REQUEST_ODDS * chance
    return odds / (1.0 - chance + odds)
