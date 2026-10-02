# stellarObjects/versionKey.py

"""
The version key (DB.6; docs/design/reproducible-galaxies.md, section 4):
22 uppercase hex digits naming exactly what a galaxy seed needs to make
the same galaxy again -- the PlanetGen release and the environment it ran
in. One helper (`version_key`) computes it from the running code, so the
galaxy, each run's history row and the log line agree.

| Part                              | Digits    |
|-----------------------------------|-----------|
| PlanetGen MAJOR, REVISION, BUILD  | 4 + 4 + 6 |
| Python major, minor, micro        | 2 + 2 + 2 |
| OS                                | 1         |
| Architecture                      | 1         |

PlanetGen 7.127.352 on Python 3.12.3, Linux, x86-64 is
`0007007F000160030C0300`. Each part is packed in its own digits, not
added up: a plain sum collides (7.127.352 and 7.128.351 both sum to 486).
"""

import platform
import sys

from ._version import __version__

KEY_DIGITS = 22
"""int: A version key's length in hex digits."""

OS_DIGITS = {"linux": 0, "windows": 1, "darwin": 3}
"""dict: `platform.system()` (lower case) -> OS digit (Boss, 2026-10-02:
0 Linux, 1 Windows, 3 macOS)."""

OTHER_UNIX_DIGIT = 2
"""int: Any other Unix or BSD (`OTHER_UNIX_SYSTEMS`)."""

OTHER_UNIX_SYSTEMS = ("freebsd", "openbsd", "netbsd", "dragonfly", "sunos", "aix", "cygwin")

ARCH_DIGITS = {
    "x86_64": 0, "amd64": 0, "x64": 0,
    "arm64": 1, "aarch64": 1, "arm64e": 1,
    "i386": 2, "i486": 2, "i586": 2, "i686": 2, "x86": 2,
    "armv6l": 3, "armv7l": 3, "armv7": 3, "arm": 3,
    "riscv64": 4,
}
"""dict: `platform.machine()` (lower case) -> architecture digit: 0
x86-64 (the default), 1 ARM64 (Apple Silicon included), 2 32-bit x86, 3
32-bit ARM, 4 RISC-V 64."""

UNKNOWN_DIGIT = 0xF
"""int: An OS or architecture not listed."""


def os_digit(system=None):
    """The OS digit for `system` (`platform.system()` by default)."""
    name = (platform.system() if system is None else system).strip().lower()
    if name in OS_DIGITS:
        return OS_DIGITS[name]
    if name.startswith(OTHER_UNIX_SYSTEMS):
        return OTHER_UNIX_DIGIT
    return UNKNOWN_DIGIT


def arch_digit(machine=None):
    """The architecture digit for `machine` (`platform.machine()` by
    default)."""
    name = (platform.machine() if machine is None else machine).strip().lower()
    return ARCH_DIGITS.get(name, UNKNOWN_DIGIT)


def _release_parts(version):
    parts = str(version).split(".")
    if len(parts) != 3 or not all(part.isdigit() for part in parts):
        raise ValueError(f"a PlanetGen version is MAJOR.REVISION.BUILD, not {version!r}")
    return tuple(int(part) for part in parts)


def version_key(version=None, python=None, system=None, machine=None):
    """
    The 22-hex-digit version key.

    Args:
        version (str, optional): `MAJOR.REVISION.BUILD`; this release's by
            default.
        python (tuple, optional): `(major, minor, micro)`; the running
            Python's by default.
        system (str, optional): `platform.system()`'s answer.
        machine (str, optional): `platform.machine()`'s answer.

    Returns:
        str: 22 uppercase hex digits.

    Raises:
        ValueError: If a part doesn't fit its digits.
    """
    major, revision, build = _release_parts(__version__ if version is None else version)
    py_major, py_minor, py_micro = tuple(sys.version_info[:3] if python is None else python)
    widths = ((major, 4), (revision, 4), (build, 6), (py_major, 2), (py_minor, 2), (py_micro, 2),
              (os_digit(system), 1), (arch_digit(machine), 1))
    for value, digits in widths:
        if not 0 <= value < 16 ** digits:
            raise ValueError(f"{value} doesn't fit in {digits} hex digits")
    return "".join(f"{value:0{digits}X}" for value, digits in widths)


def python_version():
    """The running Python's version, `3.12.3`."""
    return platform.python_version()


def platform_name():
    """The OS and architecture as `platform` reports them, `Linux x86_64`."""
    return f"{platform.system()} {platform.machine()}".strip()


def current():
    """`{"version_key", "planetgen_version", "python_version", "platform"}`
    of the running code: what the galaxy and each run's history row
    store."""
    return {
        "version_key": version_key(),
        "planetgen_version": __version__,
        "python_version": python_version(),
        "platform": platform_name(),
    }


def run_line(galaxy_seed, run):
    """
    The line every generation run writes first (OPS.10): `Galaxy seed
    <32 hex>, PlanetGen <version> (<version key>), run: <run>`.

    Args:
        galaxy_seed (bytes | None): The galaxy's 16-byte seed; `None`
            when it has none yet (never planned, or a run that doesn't
            touch the galaxy).
        run (str): What runs, such as `galaxy --ring 3`.
    """
    seed = bytes(galaxy_seed).hex().upper() if galaxy_seed is not None else "none yet"
    return f"Galaxy seed {seed}, PlanetGen {__version__} ({version_key()}), run: {run}"
