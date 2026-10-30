"""OPS.73: every constant in `tuning.py` has a reader, so a knob nobody reads cannot sit there looking like it does something."""

import ast
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TUNING = ROOT / "planetgen" / "tuning.py"

WAITING_FOR_ITS_READER = {
    "STELLAR_NS_MASS_GAUSSIANS": "GEN.930 (the neutron star masses follow from the star model; the sourced Gaussians are its input)",
}
"""dict: The constants with no reader yet and the item that wires each; an entry goes when its item merges."""


def _constants():
    tree = ast.parse(TUNING.read_text(encoding="utf-8"))
    names = []
    for node in tree.body:
        targets = [node.target] if isinstance(node, ast.AnnAssign) else node.targets if isinstance(node, ast.Assign) else []
        names += [target.id for target in targets if isinstance(target, ast.Name) and target.id.isupper()]
    return names


def _names_read():
    """Every name read, attribute taken, name imported or name given as a string in the package and the tests."""
    read = set()
    for path in [*(ROOT / "planetgen").rglob("*.py"), *(ROOT / "tests").glob("*.py")]:
        if path == pathlib.Path(__file__).resolve():
            continue
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load):
                read.add(node.id)
            elif isinstance(node, ast.Attribute):
                read.add(node.attr)
            elif isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value.isidentifier():
                read.add(node.value)
            elif isinstance(node, ast.ImportFrom):
                read.update(alias.name for alias in node.names)
    return read


def test_no_tuning_constant_lacks_a_reader():
    read = _names_read()
    unread = [name for name in _constants() if name not in read and name not in WAITING_FOR_ITS_READER]
    assert not unread, f"tuning constants nothing reads (delete them, or wire them): {unread}"


def test_the_waiting_list_holds_only_constants_that_still_wait():
    read = _names_read()
    names = set(_constants())
    stale = [name for name in WAITING_FOR_ITS_READER if name not in names or name in read]
    assert not stale, f"WAITING_FOR_ITS_READER names constants that are gone or now read: {stale}"
