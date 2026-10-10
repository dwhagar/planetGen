"""UX.84: every long step is registered, and a registered step that is predicted to pass 15
seconds (or has no history) draws a bar while a short one draws nothing."""

import ast
import pathlib

import pytest

from planetgen import tuning
from planetgen.generation import steps

PACKAGE = pathlib.Path(steps.__file__).resolve().parents[1]

# Where each call keeps its kind: a position, or a keyword.
KIND_ARGUMENT = {"Step": 1, "hooked": 2, "worker_step": 3}


def _kinds_used():
    found = {}
    for path in PACKAGE.rglob("*.py"):
        if path == pathlib.Path(steps.__file__).resolve():
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = node.func.attr if isinstance(node.func, ast.Attribute) else getattr(node.func, "id", None)
            kind = None
            if name in KIND_ARGUMENT and len(node.args) > KIND_ARGUMENT[name]:
                kind = node.args[KIND_ARGUMENT[name]]
            elif name == "StageProgress":
                kind = next((kw.value for kw in node.keywords if kw.arg == "kind"), None)
            if isinstance(kind, ast.Constant) and isinstance(kind.value, str):
                found.setdefault(kind.value.split(":")[0], []).append(f"{path.relative_to(PACKAGE)}:{node.lineno}")
    return found


def test_every_step_kind_in_the_code_is_registered():
    unregistered = {kind: where for kind, where in _kinds_used().items()
                    if kind not in steps.STEP_KINDS and kind != "stages"}
    assert not unregistered, f"add these to steps.STEP_KINDS (so they are timed and get a bar): {unregistered}"


def test_the_registry_holds_no_kind_nothing_uses():
    used = set(_kinds_used()) | {"sector", "stages"}      # `sector_bar` fixes "sector"; StageProgress defaults to "stages"
    assert set(steps.STEP_KINDS) <= used, f"unused kinds: {sorted(set(steps.STEP_KINDS) - used)}"


class _Stats:
    """Recorded speeds: `seconds` for any amount of work, or none at all."""

    def __init__(self, seconds):
        self.seconds = seconds

    def pool_rate(self, kind, workers):
        """Units a second, for a 100 unit step to take `seconds`."""
        return None if self.seconds is None else 100.0 / self.seconds

    def record(self, *args, **kwargs):
        pass

    def flush(self):
        pass


class _Display:
    def __init__(self):
        self.tasks = []

    def add_task(self, description, total=None, **kwargs):
        self.tasks.append(description)
        return len(self.tasks)

    def update(self, *_a, **_k):
        pass

    def remove_task(self, *_a):
        pass


@pytest.mark.parametrize("kind", sorted(steps.STEP_KINDS))
def test_a_step_predicted_over_the_limit_draws_its_bar_at_once(kind):
    display = _Display()
    with steps.Step("long", kind, 100, stats=_Stats(tuning.PROGRESS_BAR_SECONDS * 4), progress=display):
        assert display.tasks == ["long"]


@pytest.mark.parametrize("kind", sorted(steps.STEP_KINDS))
def test_a_step_with_no_recorded_speed_counts_as_long_and_draws(kind):
    display = _Display()
    with steps.Step("unknown", kind, 100, stats=_Stats(None), progress=display):
        assert display.tasks == ["unknown"]


@pytest.mark.parametrize("kind", sorted(steps.STEP_KINDS))
def test_a_step_predicted_short_draws_nothing(kind):
    display = _Display()
    with steps.Step("short", kind, 100, stats=_Stats(1.0), progress=display):
        assert display.tasks == []


def test_hooked_hands_a_callers_own_callback_through_and_draws_nothing_without_a_display():
    seen = []
    with steps.hooked(lambda *a: seen.append(a), "x", "containment", 3) as report:
        report("Sectors", 1, 3)
    assert seen == [("Sectors", 1, 3)]
    with steps.hooked(None, "x", "containment", 3) as report:      # no open display: nothing to draw on
        report("Sectors", 1, 3)


def test_hooked_draws_a_step_on_the_display_it_is_given(monkeypatch):
    display = _Display()
    with steps.hooked(None, "Containment: sectors", "containment", 3, progress=display, force=True) as report:
        report("Sectors", 1, 3)
    assert display.tasks == ["Containment: sectors"]
