# planetgen/db/render.py

"""
On-demand rendering of a stored star system's text.

No page text is stored (schema v29): `describe.describe` builds the system's Document from the
database rows (`load_star_system`) and `describe.render` writes it as text, Markdown or wikitext
when someone asks, for the code views, the downloads, the command line and the wiki upload.
`render_system_sections` gives the page's Markdown split per star/planet/moon/belt/comet for
`web/lib/systempage.py`'s expandable system list until the page renders the Document itself
(UX.105).
"""

from contextlib import contextmanager

from planetgen import describe
from planetgen.db.store import load_star_system

FORMATS = tuple(describe.RENDERERS)
"""tuple: The formats a system can be rendered in (`describe.RENDERERS`)."""


def _check_format(fmt):
    if fmt not in FORMATS:
        raise ValueError(f"unknown format {fmt!r} (expected one of {', '.join(FORMATS)})")


@contextmanager
def _format(system, fmt):
    """Sets the system's shared `SystemConfig.MARKDOWN` flag for the
    duration of a render, restoring it afterward. Every reconstructed
    object shares that one config, so this switches the whole page."""
    _check_format(fmt)
    config = system.system_config
    original = config.MARKDOWN
    config.MARKDOWN = fmt == "markdown"
    try:
        yield
    finally:
        config.MARKDOWN = original


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


def _without_header(paragraphs):
    """Drops a leading `#`-style heading paragraph -- the system list
    already shows each body's name on its own row."""
    if paragraphs and paragraphs[0].lstrip().startswith("#"):
        return paragraphs[1:]
    return paragraphs


def _body_markdown(body):
    """One planet's, moon's, belt's or comet's own Markdown, without its
    heading and (for a planet) without its moons' sections, which
    `Planet.to_paragraph_list` appends after its own."""
    paragraphs = body.to_paragraph_list()
    moon_paragraph_count = sum(len(moon.to_paragraph_list()) for moon in getattr(body, "moons", []))
    if moon_paragraph_count:
        paragraphs = paragraphs[:-moon_paragraph_count]
    return "\n\n".join(_without_header(paragraphs))


def render_system_sections(conn, star_system_id):
    """
    Renders one stored system's page as Markdown, split into the pieces
    `html/system.py`'s system list shows under each row.

    Args:
        conn (Connection): An open connection.
        star_system_id (int): The `star_systems.id` to render.

    Returns:
        dict: `overview` (the pair's own data table for a binary, the
            system summary, a close pair's age sentence and any system
            flavor text, in page order), plus `stars`, `planets`, `moons`,
            `belts` and `comets`, each mapping a row `id` (as a string, so
            the dict survives a JSON round trip unchanged) to that row's
            own Markdown.

    Raises:
        ValueError: If no such system exists.
    """
    system = load_star_system(conn, star_system_id)
    sections = {"overview": "", "stars": {}, "planets": {}, "moons": {}, "belts": {}, "comets": {}}

    with _format(system, "markdown"):
        overview = []
        if system.binary_type == "close":
            data_block, age_sentence = system.star.to_paragraph_list()[:2]
            overview += [data_block, system.summary_paragraph(), age_sentence]
        elif system.binary_type == "wide":
            overview += [system.wide_binary.to_paragraph_list()[0], system.summary_paragraph()]
        else:
            overview.append(system.summary_paragraph())
        if system.system_flavor_text:
            overview.append(f"Sensors show {system.system_flavor_text}")
        sections["overview"] = "\n\n".join(overview)

        for star in system.stars:
            sections["stars"][str(star.db_id)] = "\n\n".join(star.to_paragraph_list())

        for body in system.planets + system.secondary_planets:
            if body.body_type == "a":
                sections["belts"][str(body.db_id)] = _body_markdown(body)
                continue
            sections["planets"][str(body.db_id)] = _body_markdown(body)
            for moon in body.moons:
                sections["moons"][str(moon.db_id)] = _body_markdown(moon)

        for comet in system.comets + system.secondary_comets:
            sections["comets"][str(comet.db_id)] = _body_markdown(comet)

    return sections
