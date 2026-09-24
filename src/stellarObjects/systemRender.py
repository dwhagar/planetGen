# stellarObjects/systemRender.py

"""
On-demand rendering of a stored star system's wiki page text.

Since schema v28 no page text is stored (see `schema.sql`'s "v28" header
note): both formats are rebuilt from the database rows whenever they are
needed -- `_db.load_star_system` reconstructs the generation object graph
and `StarSystem.__str__` renders it, exactly as generation itself used to
right before saving. `render_system_text` gives the whole page, for the
Wikitext/Markdown code view and for wiki upload; `render_system_sections`
gives the same Markdown split per star/planet/moon/belt/comet, for
`html/system.py`'s expandable system list.
"""

from contextlib import contextmanager

from ._db import load_star_system

FORMATS = ("wikitext", "markdown")
"""tuple: The page formats `render_system_text` accepts."""


@contextmanager
def _format(system, fmt):
    """Sets the system's shared `SystemConfig.MARKDOWN` flag for the
    duration of a render, restoring it afterward. Every reconstructed
    object shares that one config, so this switches the whole page."""
    if fmt not in FORMATS:
        raise ValueError(f"unknown format {fmt!r} (expected one of {', '.join(FORMATS)})")
    config = system.system_config
    original = config.MARKDOWN
    config.MARKDOWN = fmt == "markdown"
    try:
        yield
    finally:
        config.MARKDOWN = original


def render_star_system(system, fmt):
    """
    Renders an already-loaded (or freshly generated) `StarSystem` as a full
    wiki page.

    Args:
        system (StarSystem): The system to render.
        fmt (str): `"wikitext"` or `"markdown"`.

    Returns:
        str: The page text.

    Raises:
        ValueError: If `fmt` isn't one of `FORMATS`.
    """
    with _format(system, fmt):
        return str(system)


def render_system_text(conn, star_system_id, fmt):
    """
    Renders one stored system's full wiki page from its database rows.

    Args:
        conn (Connection): An open connection.
        star_system_id (int): The `star_systems.id` to render.
        fmt (str): `"wikitext"` or `"markdown"`.

    Returns:
        str: The page text.

    Raises:
        ValueError: If no such system exists, or `fmt` is unknown.
    """
    if fmt not in FORMATS:
        raise ValueError(f"unknown format {fmt!r} (expected one of {', '.join(FORMATS)})")
    return render_star_system(load_star_system(conn, star_system_id), fmt)


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
