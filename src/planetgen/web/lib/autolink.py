# planetgen/web/lib/autolink.py

"""
The auto-linker: makes sure a page's text links to the help pages that explain the words it uses.

The terms come from `help_registry.terms()`, so each part of the program says which of its words have a help page. The rules (Boss,
2026-10-31):

* The count is per page, and per term: the first `LIMIT` (2) mentions of a term are linked.
* A link already on the page counts. A term whose help page the page already links `LIMIT` times or more is left alone; with one link
  already, one more mention is linked.
* Text inside an existing link is never changed, and neither is heading text, a button, a form control, a script, a style or the site's
  own header and menus. A table cell is linked like any other text.
* A term matches whole words only, case exactly, and the longer of two overlapping terms wins. A term with an `after` prefix ("Class "
  before a bare life class number) matches only when that text comes right before it, and the prefix stays outside the link.
* A page does not link to itself.
* A page's markup can switch the linker off for a stretch with the comments `<!--autolink:off-->` and `<!--autolink:on-->`.

`link_page(html, terms, url_of, own_url=None)` changes only what is between `<main` and `</main>`. `link(fragment, ...)` works on a
fragment. Pure functions, no Flask.
"""

import re
from html import escape
from html.parser import HTMLParser

LIMIT = 2
"""int: How many links to a term's help page a page should hold."""

SKIP = frozenset({
    "a", "h1", "h2", "h3", "h4", "h5", "h6", "button", "textarea", "select", "option", "label", "script", "style", "title",
    "svg", "head", "header", "nav", "footer", "sl-menu", "sl-menu-item", "sl-dropdown", "sl-button", "sl-select", "sl-input", "code", "pre",
})
"""frozenset: The tags whose text is never linked (and the content of which is not scanned)."""

OFF = "autolink:off"
ON = "autolink:on"
"""str: The HTML comments a page's markup can use to keep the linker out of a stretch of text (`<!--autolink:off-->` ... `<!--autolink:on-->`), for
text that is a name or a quotation and must stay as written. They nest."""

_VOID = frozenset({"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "source", "track", "wbr"})


class _Linker(HTMLParser):
    def __init__(self, matchers, urls, own_url, existing):
        super().__init__(convert_charrefs=False)
        self.out = []
        self.matchers = matchers          # [(term, compiled regex)]
        self.urls = urls                  # {term text: url}
        self.own_url = own_url
        self.left = {}                    # {term text: links still to add}
        for term, _rx in matchers:
            url = urls[term.text]
            self.left[term.text] = 0 if url == own_url else max(0, LIMIT - existing.get(url, 0))
        self.skip = 0
        self.stack = []

    def handle_starttag(self, tag, attrs):
        self.out.append(self.get_starttag_text())
        if tag not in _VOID:
            self.stack.append(tag)
            if tag in SKIP:
                self.skip += 1

    def handle_startendtag(self, tag, attrs):
        self.out.append(self.get_starttag_text())

    def handle_endtag(self, tag):
        self.out.append(f"</{tag}>")
        if tag in self.stack:
            while self.stack:
                top = self.stack.pop()
                if top in SKIP:
                    self.skip -= 1
                if top == tag:
                    break

    def handle_data(self, data):
        self.out.append(data if self.skip or not any(self.left.values()) else self._link(data))

    def handle_entityref(self, name):
        self.out.append(f"&{name};")

    def handle_charref(self, name):
        self.out.append(f"&#{name};")

    def handle_comment(self, data):
        self.out.append(f"<!--{data}-->")
        if data.strip() == OFF:
            self.skip += 1
        elif data.strip() == ON and self.skip:
            self.skip -= 1

    def handle_decl(self, decl):
        self.out.append(f"<!{decl}>")

    def handle_pi(self, data):
        self.out.append(f"<?{data}>")

    def _link(self, text):
        found = []
        for term, rx in self.matchers:
            if not self.left[term.text]:
                continue
            for match in rx.finditer(text):
                found.append((match.start("term"), match.end("term"), term))
        # The longer term wins an overlap; of equal length the earlier one.
        found.sort(key=lambda f: (-(f[1] - f[0]), f[0]))
        taken, chosen = [], []
        for start, end, term in found:
            if any(start < b and a < end for a, b in taken):
                continue
            taken.append((start, end))
            chosen.append((start, end, term))
        chosen.sort()
        pieces, at = [], 0
        for start, end, term in chosen:
            if not self.left[term.text]:
                continue
            self.left[term.text] -= 1
            pieces.append(text[at:start])
            pieces.append(f'<a href="{escape(self.urls[term.text], quote=True)}">{text[start:end]}</a>')
            at = end
        pieces.append(text[at:])
        return "".join(pieces)


class _Counter(HTMLParser):
    """Counts the links the fragment already has, by address."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.count = {}

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.count[href] = self.count.get(href, 0) + 1


def _matcher(term):
    after = f"(?<={re.escape(term.after)})" if term.after else ""
    return re.compile(rf"{after}(?<![\w-])(?P<term>{re.escape(term.text)})(?![\w-])")


def link(fragment, terms, url_of, own_url=None):
    """
    `fragment` with help links added for `terms` (`help_registry.Term`s), `url_of(term)` giving a term's address. Returns the fragment
    unchanged when no term appears in it.
    """
    terms = [t for t in terms if t.text in fragment]
    if not terms:
        return fragment
    urls, seen = {}, set()
    matchers = []
    for term in terms:
        if term.text in seen:
            continue
        seen.add(term.text)
        urls[term.text] = url_of(term)
        matchers.append((term, _matcher(term)))
    counter = _Counter()
    counter.feed(fragment)
    counter.close()
    linker = _Linker(matchers, urls, own_url, counter.count)
    linker.feed(fragment)
    linker.close()
    return "".join(linker.out)


def link_page(html, terms, url_of, own_url=None):
    """`html` with the part between `<main` and `</main>` linked by `link`; a page with no `<main>` is returned as it came."""
    start = html.find("<main")
    end = html.rfind("</main>")
    if start < 0 or end < start:
        return html
    return html[:start] + link(html[start:end], terms, url_of, own_url) + html[end:]
