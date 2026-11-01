# planetgen/web/lib/autolink.py

"""
The auto-linker: makes sure a page's text links to the help pages that explain the words it uses.

The terms come from `help_registry.terms()`, so each part of the program says which of its words have a help page. The rules (Boss,
2026-10-31):

* The count is per page, and per term: the first `LIMIT` (2) mentions of a term are linked.
* A link already on the page counts, whatever its text and whether or not its address carries a `#fragment`, a `?query` or a closing
  slash. A term whose help page the page already links `LIMIT` times or more is left alone; with one link already, one more mention is linked.
* Text inside an existing link is never changed, and neither is heading text, a button, a form control, a script, a style, an
  attribute, a comment or the site's own header and menus. A table cell is linked like any other text.
* A term matches whole words only, with the case exactly as registered. It is not a word when a letter, digit, underscore or hyphen touches
  it on either side (`14c`, `4cm`, `a-4c`), when it is the tail of a decimal (`1.4c`, `1,4c`), or when it is the head of one (`Class 2.5`).
  Those neighbors count through entities (`&eacute;4c`) and through inline tags (`x<b>4c</b>`). A block tag (a paragraph, a cell) is a gap.
* A bare id followed by a unit symbol is a measurement, not a class (`Class 2 m`, `Class 3a kg`); the symbols are the unit registry's.
* The longer of two overlapping terms wins.
* A term with an `after` prefix matches only after that exact text, as a word of its own (`Subclass 2` and `XClass 2` do not count), or inside
  a list that starts with it ("Classes 3a, 4b and 5"); the prefix stays outside the link. Every life class id has the prefix "Class ", so a
  measurement ("4c" as a speed, "4 C", "3 g", "2 m", "15 Mpc") is never taken for a class.
* A page does not link to itself.
* Markup the reader cannot give back exactly as it came in (an unclosed comment, an entity with no semicolon) is left as written, with no links.
* A page's markup can switch the linker off for a stretch with the comments `<!--autolink:off-->` and `<!--autolink:on-->`; the pairs nest, an
  `on` with no `off` does nothing, and an `off` never closed switches it off to the end of the page.

`link_page(html, terms, url_of, own_url=None)` changes only what is between `<main` and `</main>`. `link(fragment, ...)` works on a
fragment. Pure functions, no Flask.
"""

import re
from html import escape, unescape
from html.parser import HTMLParser
from urllib.parse import urlsplit

LIMIT = 2
"""int: How many links to a term's help page a page should hold."""

SKIP = frozenset({
    "a", "h1", "h2", "h3", "h4", "h5", "h6", "button", "textarea", "select", "option", "label", "script", "style", "title",
    "svg", "math", "head", "header", "nav", "footer", "sl-menu", "sl-menu-item", "sl-dropdown", "sl-button", "sl-select", "sl-input",
    "code", "pre", "noscript", "template", "iframe", "object",
})
"""frozenset: The tags whose text is never linked (and the content of which is not scanned)."""

INLINE = frozenset({
    "a", "abbr", "b", "bdi", "bdo", "cite", "code", "data", "dfn", "em", "i", "kbd", "mark", "q", "s", "samp", "small", "span", "strong",
    "sub", "sup", "time", "u", "var", "wbr", "sl-tooltip", "sl-icon",
})
"""frozenset: The tags that do not break a word: `x<b>4c</b>` is the word `x4c`."""

OFF = "autolink:off"
ON = "autolink:on"
"""str: The HTML comments a page's markup can use to keep the linker out of a stretch of text (`<!--autolink:off-->` ... `<!--autolink:on-->`), for
text that is a name or a quotation and must stay as written."""

CONTEXT = 24
"""int: How many characters of the text before and after a stretch the matcher can see (for a prefix or a unit split by an inline tag)."""

_VOID = frozenset({"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "source", "track", "wbr"})


class _Token:
    """One piece of the fragment: a tag or comment (`raw` as written, `tag` its lower-case name, `kind` "start", "end" or "other"),
    or a stretch of text (`segments`: `("d", raw)` for characters and `("e", raw)` for an entity)."""

    __slots__ = ("kind", "raw", "tag", "segments", "skip")

    def __init__(self, kind, raw="", tag=""):
        self.kind, self.raw, self.tag, self.segments, self.skip = kind, raw, tag, [], False


class _Reader(HTMLParser):
    """Splits a fragment into tokens and marks the text that must not be linked."""

    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.tokens = []
        self.stack = []
        self.off = 0

    def _text(self):
        if not self.tokens or self.tokens[-1].kind != "text":
            token = _Token("text")
            token.skip = bool(self.off or any(tag in SKIP for tag in self.stack))
            self.tokens.append(token)
        return self.tokens[-1]

    def handle_starttag(self, tag, attrs):
        self.tokens.append(_Token("start", self.get_starttag_text(), tag))
        if tag not in _VOID:
            self.stack.append(tag)

    def handle_startendtag(self, tag, attrs):
        self.tokens.append(_Token("other", self.get_starttag_text(), tag))

    def handle_endtag(self, tag):
        self.tokens.append(_Token("end", f"</{tag}>", tag))
        if tag in self.stack:
            while self.stack and self.stack.pop() != tag:
                pass

    def handle_data(self, data):
        self._text().segments.append(("d", data))

    def handle_entityref(self, name):
        self._text().segments.append(("e", f"&{name};"))

    def handle_charref(self, name):
        self._text().segments.append(("e", f"&#{name};"))

    def handle_comment(self, data):
        self.tokens.append(_Token("other", f"<!--{data}-->", "!--"))
        mark = data.strip()
        if mark == OFF:
            self.off += 1
        elif mark == ON and self.off:
            self.off -= 1

    def handle_decl(self, decl):
        self.tokens.append(_Token("other", f"<!{decl}>", "!"))

    def handle_pi(self, data):
        self.tokens.append(_Token("other", f"<?{data}>", "?"))


def _address(href):
    """An address as the linker compares it: its path without a closing slash (a `#fragment` and a `?query` are dropped)."""
    path = urlsplit(href.strip()).path
    return path.rstrip("/") or "/"


def _existing(tokens):
    """The links the fragment already has, counted by address."""
    count = {}
    for token in tokens:
        if token.kind == "start" and token.tag == "a":
            match = re.search(r"""href\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s>]+))""", token.raw, re.I)
            if match:
                href = unescape(next(group for group in match.groups() if group is not None))
                count[_address(href)] = count.get(_address(href), 0) + 1
    return count


def _plain(token):
    """The text of a text token as a reader sees it (entities decoded) and, for each segment, where it lies in that text and in the raw one."""
    text, raw, parts = [], 0, []
    pos = 0
    for kind, source in token.segments:
        shown = source if kind == "d" else unescape(source)
        parts.append((kind, pos, pos + len(shown), raw, raw + len(source)))
        text.append(shown)
        pos += len(shown)
        raw += len(source)
    return "".join(text), parts


def _raw_offset(parts, index, end):
    """The offset in the raw text of the plain offset `index`, or `None` when it falls inside an entity (`end` says it ends a match)."""
    for kind, a, b, ra, rb in parts:
        if kind == "d" and (a <= index < b or (end and index == b)):
            return ra + (index - a)
        if kind == "e":
            if index == a:
                return ra
            if end and index == b:
                return rb
    if parts and index == parts[-1][2]:
        return parts[-1][4]
    return None


def _neighbor(tokens, at, step, plain_of, want):
    """Up to `want` characters of the text beside token `at`, looking through inline tags and comments (`step` -1 before, +1 after); a block tag, or the
    end of the fragment, gives nothing."""
    gathered = ""
    i = at + step
    while 0 <= i < len(tokens) and len(gathered) < want:
        token = tokens[i]
        if token.kind == "text":
            shown = plain_of(i)
            gathered = shown[-(want - len(gathered)):] + gathered if step < 0 else gathered + shown[:want - len(gathered)]
        elif token.tag not in INLINE and token.tag != "!--":
            break
        i += step
    return gathered


ITEM = r"\d+[a-z]?"
SEPARATOR = r"(?:, and |, or |, | and/or | and | or | to |/)"
"""str: What may sit between the items of a list that follows a prefix ("Class 3a, 4b and 5"): a term in such a list is still a class."""


NOT_UNITS = frozenset({"in", "turn", "min", "(fraction)", "mi", "ms (Unix)", "s (Unix)", "yr (Unix)"})
"""frozenset: Unit symbols that are also ordinary words, so "Class 2 in the" is not read as a measurement."""

EXTRA_UNITS = ("mg", "µg", "μg", "Sv", "mSv", "Gy", "mol", "Msun", "Mjup", "Rjup", "Re", "Rj", "Gya", "Mya", "dex", "mag", "ppb", "°", "×")
"""tuple: Symbols written in the pages that the unit registry does not list."""


def _unit_pattern():
    """A pattern for a unit symbol after a number: every symbol of the unit registry (`physics/unit_registry.py`) that is not an ordinary word,
    and `EXTRA_UNITS`, longest first. Read when a term is first compiled."""
    try:
        from planetgen.physics import unit_registry
        symbols = {unit.symbol for unit in unit_registry.UNITS.values()}
    except Exception:    # the linker works without the registry, with the extra symbols only
        symbols = set()
    symbols = sorted((symbols - NOT_UNITS) | set(EXTRA_UNITS), key=len, reverse=True)
    return "|".join(re.escape(symbol) for symbol in symbols if symbol)


_units = []


def _matcher(term):
    """The pattern of a term: whole word only, not the end of a decimal before it or the start of one after it, and, for a term that is a bare id,
    not a number with a unit (`2 m`, `4c` as a speed is already out for want of the prefix, `3a` in `3a kg`)."""
    lead = r"(?<![\w-])(?<!\d[.,])"
    tail = r"(?![\w-])(?![.,]\d)"
    if term.after:
        if not _units:
            _units.append(_unit_pattern())
        tail += rf"(?![ \u00a0\u202f\u2009]?(?:{_units[0]})(?![\w]))"
    return re.compile(rf"{lead}(?P<term>{re.escape(term.text)}){tail}")


def _prefix(term):
    """For a term with an `after` prefix, the pattern the text before it must end with: the prefix as a word of its own (also in the plural,
    "Classes "), then any list of items and separators ("Class 3a, 4b and "). `None` for a term with no prefix."""
    if not term.after:
        return None
    base = term.after.rstrip()
    return re.compile(rf"(?<![\w-]){re.escape(base)}(?:es)?{re.escape(term.after[len(base):])}(?:{ITEM}{SEPARATOR})*$")


def link(fragment, terms, url_of, own_url=None):
    """
    `fragment` with help links added for `terms` (`help_registry.Term`s), `url_of(term)` giving a term's address. Returns the fragment
    unchanged when no term can appear in it.
    """
    if "&" not in fragment:     # an entity can spell a term (`&#52;c`), so only a fragment without one is screened by plain search
        terms = [t for t in terms if t.text in fragment and (not t.after or t.after.strip() in fragment)]
    if not terms:
        return fragment
    reader = _Reader()
    reader.feed(fragment)
    reader.close()
    tokens = reader.tokens
    if "".join(token.raw or "".join(source for _kind, source in token.segments) for token in tokens) != fragment:
        return fragment       # markup the reader cannot give back exactly (an unclosed comment, a bare `&amp`, an upper-case end tag): leave it as written
    existing = _existing(tokens)
    own = _address(own_url) if own_url else None
    urls, matchers, left = {}, [], {}
    for term in terms:
        if term.text in urls:
            continue
        urls[term.text] = url_of(term)
        address = _address(urls[term.text])
        left[term.text] = 0 if address == own else max(0, LIMIT - existing.get(address, 0))
        matchers.append((term, _matcher(term), _prefix(term)))
    cache = {}

    def plain_of(i):
        if i not in cache:
            cache[i] = _plain(tokens[i])
        return cache[i][0]

    out = []
    for i, token in enumerate(tokens):
        if token.kind != "text":
            out.append(token.raw)
            continue
        raw = "".join(source for _kind, source in token.segments)
        if token.skip or not any(left.values()):
            out.append(raw)
            continue
        text, parts = _plain(token)
        before = _neighbor(tokens, i, -1, plain_of, CONTEXT)
        after = _neighbor(tokens, i, +1, plain_of, CONTEXT)
        whole = before + text + after
        found = []
        for term, rx, prefix in matchers:
            if not left[term.text]:
                continue
            for match in rx.finditer(whole):
                start, end = match.start("term") - len(before), match.end("term") - len(before)
                if start < 0 or end > len(text):
                    continue
                if prefix is not None and not prefix.search(whole[:match.start("term")]):
                    continue
                found.append((start, end, term))
        found.sort(key=lambda f: (-(f[1] - f[0]), f[0]))      # the longer term wins an overlap; of equal length the earlier
        taken, chosen = [], []
        for start, end, term in found:
            if any(start < b and a < end for a, b in taken):
                continue
            taken.append((start, end))
            chosen.append((start, end, term))
        pieces, at = [], 0
        for start, end, term in sorted(chosen):
            if not left[term.text]:
                continue
            raw_start, raw_end = _raw_offset(parts, start, False), _raw_offset(parts, end, True)
            if raw_start is None or raw_end is None or raw_start < at:
                continue
            left[term.text] -= 1
            pieces.append(raw[at:raw_start])
            pieces.append(f'<a href="{escape(urls[term.text], quote=True)}">{raw[raw_start:raw_end]}</a>')
            at = raw_end
        pieces.append(raw[at:])
        out.append("".join(pieces))
    return "".join(out)


def link_page(html, terms, url_of, own_url=None):
    """`html` with the part between `<main` and `</main>` linked by `link`; a page with no `<main>` is returned as it came."""
    start = html.find("<main")
    end = html.rfind("</main>")
    if start < 0 or end < start:
        return html
    return html[:start] + link(html[start:end], terms, url_of, own_url) + html[end:]
