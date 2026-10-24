# tests/first_fill.py

"""A test that wants the sectors it fills to hold exactly what it asked for turns off the first fill's draw down to the
backfill's deepest floor (GEN.633), the way a test of the backfill turns the backfill off."""

from planetgen.generation import run_galaxy


def without_the_deepest_floor(monkeypatch):
    """Makes `run_galaxy._fill_context` skip `_take_to_the_deepest_floor` for the test."""
    monkeypatch.setattr(run_galaxy, "_take_to_the_deepest_floor", lambda *_args, **_kwargs: None)
