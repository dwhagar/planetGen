# planetgen/db/render.py

"""
On-demand rendering of a stored star system's text.

No page text is stored (schema v29): `describe.describe` builds the system's Document from the
database rows (`load_star_system`) and `describe.render` writes it as text or Markdown
when someone asks, for the code views, the downloads and the command line.
`describe_system` is the Document itself, for the web page, which writes it as HTML.
"""

from planetgen import describe
from planetgen.db.store import load_star_system

FORMATS = tuple(describe.RENDERERS)
"""tuple: The formats a system can be rendered in (`describe.RENDERERS`)."""


def _check_format(fmt):
    if fmt not in FORMATS:
        raise ValueError(f"unknown format {fmt!r} (expected one of {', '.join(FORMATS)})")


def render_star_system(system, fmt):
    """
    Renders an already-loaded (or freshly generated) `StarSystem` in `fmt`.

    Args:
        system (StarSystem): The system to render.
        fmt (str): One of `FORMATS`.

    Returns:
        str: The text.

    Raises:
        ValueError: If `fmt` isn't one of `FORMATS`.
    """
    _check_format(fmt)
    return describe.render(describe.describe(system), fmt)


def render_system_text(conn, star_system_id, fmt):
    """
    Renders one stored system from its database rows.

    Args:
        conn (Connection): An open connection.
        star_system_id (int): The `star_systems.id` to render.
        fmt (str): One of `FORMATS`.

    Returns:
        str: The text.

    Raises:
        ValueError: If no such system exists, or `fmt` is unknown.
    """
    _check_format(fmt)
    return render_star_system(load_star_system(conn, star_system_id), fmt)


def describe_system(conn, star_system_id):
    """
    One stored system's Document (`describe.describe`), JSON-safe.

    Raises:
        ValueError: If no such system exists.
    """
    return describe.describe(load_star_system(conn, star_system_id))
