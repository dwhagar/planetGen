# planetgen/util/star_count.py

"""
The stars filter of a search (UX.239): how many stars a system holds, as the search page, the search box and `GET /api/search?stars=` write it.
The web layer reads the database only through the API, so the parser lives here, where both the API (`db.query.search`) and the page can use it.
"""

STAR_COUNT_NAMES = (("single", 1), ("binary", 2), ("trinary", 3), ("quaternary", 4), ("quintenary", 5))
"""tuple: The named star counts of the stars filter (UX.239): `(name, stars)`. `quintenary` and the digits above 5 stand for five or more stars as a
choice, see `STAR_COUNT_CHOICES`."""

STAR_COUNT_CHOICES = (("1", "Single"), ("2", "Binary"), ("3", "Trinary"), ("4", "Quaternary"), ("5+", "Quintenary and beyond"))
"""tuple: The choices the search form offers for the stars filter, `(value, label)`: exactly one to four stars, five or more."""


def parse_star_count(text):
    """
    The stars filter of a search (UX.239) as `(least, most)`, `most` `None` for "or more", or `None` for no filter (an empty `text`). `text` is a count
    (`3`: exactly three stars), a count and a plus (`3+`: three or more), or a name (`single`, `binary`, `trinary`, `quaternary` for exactly that many,
    `quintenary` for five or more, `quintenary+` likewise); the case is not significant.

    Raises:
        ValueError: Text that is none of these, or a count under 1.
    """
    word = (text or "").strip().lower()
    if not word:
        return None
    named = dict(STAR_COUNT_NAMES)
    plus = word.endswith("+")
    base = word[:-1].strip() if plus else word
    if base in named:
        count = named[base]
        return (count, None) if plus or base == "quintenary" else (count, count)
    if base.isdigit() and int(base) >= 1:
        return (int(base), None) if plus else (int(base), int(base))
    raise ValueError(f"stars must be a count such as 3, 3+ for three or more, or one of {', '.join(name for name, _ in STAR_COUNT_NAMES)}; got {text!r}")


def star_count_label(stars):
    """The words for a stars filter `(least, most)`: `Binary`, `Quintenary and beyond`, `7 stars or more`, `6 stars`."""
    if stars is None:
        return ""
    least, most = stars
    names = {count: name.capitalize() for name, count in STAR_COUNT_NAMES}
    if most == least:
        return names[least] if least in names and least < 5 else f"{least} stars"
    if most is None:
        return "Quintenary and beyond" if least == 5 else f"{least} stars or more"
    return f"{least} to {most} stars"
