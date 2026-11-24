# tests/air_history_off.py

"""A test of what a build makes besides the air (a moon's origin, orbit and makeup) builds the systems without the record of the air history (GEN.1745), which was most of a system's cost (TEST.831): about 60 percent of the
build of a system with moons. The record is made at the end of a body's build from the finished body and draws only on streams of its own (`air_record`), so nothing else a body holds moves and the draw stream
is the same after it (the moons of the 60 seeded systems of `test_moon_composition`: one digest of every moon's fields and the next draw, with and without it). A body gets no `air_history`; its life ages are still made."""

import contextlib

from planetgen.physics import air_record


def _no_record(*_args, **_kwargs):
    """`air_record.record_of` for a build with no record: the body gets none (`life._history_of` stores `None`)."""
    return None


def without_the_air_history(monkeypatch):
    """Makes `air_record.record_of` give no record for the test."""
    monkeypatch.setattr(air_record, "record_of", _no_record)


@contextlib.contextmanager
def air_history_off():
    """The same for a fixture that has no `monkeypatch` of its own (a module-scoped one): the record is back after the `with`."""
    import pytest
    with pytest.MonkeyPatch.context() as patch:
        without_the_air_history(patch)
        yield
