#!/usr/bin/env python3
"""
Builds the documentation website, `docs/site/`: every Markdown file in the repository's documentation as an HTML page,
so GitHub Pages (Settings > Pages > Deploy from a branch > `main` / `docs`) can serve the docs as a browsable site.

Boss (2026-11-08, 6:13 PM EST): "a script that will build HTML pages that github will serve us up for us from the
docs directory."

How it works
============
* Every `docs/**/*.md` file, and the top-level `README.md`, `INSTALL.md`, `CHANGELOG.md`, `DATA-LICENSE.md` and
  `LICENSE.md` (under `docs/site/repo/`), becomes `docs/site/<same path>.html`.
* Markdown is read with `markdown-it-py` (CommonMark plus tables and strikethrough; already a dependency of `rich`).
  Each heading gets a GitHub-style anchor, so every `#anchor` in the sources still works.
* Links: a link to another Markdown file becomes a link to its page; a link to a generated page or other file under
  `docs/` becomes a relative link to that file; a link to a file elsewhere in the repository (`src/...`) becomes a
  link to it on GitHub; a link on the web stays. Every link is relative, so the site works at any address.
* Formulas: chemical formulas in prose (`H2O`, `CO2`, `CH4`, `H2SO4`) get subscripts (H<sub>2</sub>O) outside code;
  `$...$` and `$$...$$` LaTeX is rendered by KaTeX (loaded from cdnjs, pinned to a version); the TeX stays readable
  as text when KaTeX cannot load.
* Every page has a side menu of the whole site and a search box, both read from `docs/site/nav.json` (built here), and
  a table of contents of its own headings; `docs/site/index.html` is the front page.
* Only this script writes `docs/site/`. It writes nothing outside it, and a file whose text did not change is left
  alone. `--check` fails if any page is out of date or a stale one is left behind. Like the other global documents
  (DOC.69), the TODO thread rebuilds it through `scripts/build_global_docs.py`.

Usage
=====
    python scripts/build_docs_site.py            # write docs/site/
    python scripts/build_docs_site.py --check    # fail if docs/site/ is out of date
    python scripts/build_docs_site.py --chemicals  # list every chemical formula the subscripting touched (to review)
"""

import argparse
import html
import json
import os
import posixpath
import re
import sys
from urllib.parse import unquote, urlsplit

from markdown_it import MarkdownIt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")
SITE = os.path.join(DOCS, "site")
GITHUB = "https://github.com/dwhagar/planetGen/blob/main/"
KATEX = "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9"
ROOT_FILES = ("README.md", "INSTALL.md", "CHANGELOG.md", "DATA-LICENSE.md", "LICENSE.md")
"""tuple: the repository's own top-level Markdown files, published under `site/repo/`."""

KEEP_NAMES = ("docs/REFERENCES.", "docs/references/", "DATA-LICENSE.md", "LICENSE.md", "CHANGELOG.md")
"""tuple: path fragments of the files that keep author names (`rewrite_citations.KEEP`, DOC.99); no page is built from them."""

SECTION_TITLES = {
    "": "Start here", "reference": "Reference: using and configuring the program", "deployment": "Deployment",
    "design": "Design notes", "plan": "The plan", "references": "Reference papers and datasets",
    "archive": "Archive (replaced designs)", "repo": "The repository's own files",
}
"""dict: first folder of a page to its heading in the menu."""

# --- chemical formulas -----------------------------------------------------------------------------------------

ELEMENTS = set(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb "
    "Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au "
    "Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm".split())
MOLECULES = {"H2", "N2", "O2", "F2", "Cl2", "Br2", "I2", "S2", "S8", "P4", "O3", "H2O", "D2O"}
"""set: formulas made of one element with a digit, which are only subscripted when named here (`C4`, `P100` are not formulas)."""
UPPERCASE_FORMULAS = set(
    "CO2 N2 CH4 O2 NH3 H2 H2O H2S SO2 O3 N2O C2H6 H2SO4 NO2 PH3 CH3OH S2 H2CO C2H2 S8 F2 K2O N2O4 C2H4 P2O5 HNO3 H2O2 SO3 "
    "CS2 CF4 C6H6 HO2 NH2 N2H4 C3H8 C4H10 C2N2 C3O2 C4H2 CO3 SO4 NO3 HC3N CH2O CH3SH CH3CN NH4 H2CO3 HCO3 CH3 C2H C2H3 "
    "C2H5 CH3O2 P4 S2O S2O2 H2COH C4N2 C6H5 PH2 SF6 NH4SH HSO3 HSO4 HNO2 CH2 CF2 CF3 H3NO4 OF2".split())
"""set: formulas of capital letters and digits only that are subscripted; every other all-capital string with digits is
a label (`BC13`, `F5V`, `IF6`, `HP22`) or a name, and a rule cannot tell them apart. Add a species here when a page has one."""
NOT_FORMULAS = {"Tc(H2)", "F(CO2)", "F(N2)", "NH32H2O", "FeO2He", "Na2He", "FeO1", "H13CO", "SO42", "Bi2S3", "BrONO2", "MnAu2",
                "La2CuO4", "CaAl12O19", "WFC3", "CF3H"}
"""set: strings that pass the element check but are not formulas worth subscripting."""
_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_FORMULA = re.compile(r"(?<![\w.\-/#$&;])((?:[A-Z][a-z]?\d*){1,10}(?:\((?:[A-Z][a-z]?\d*){1,6}\)\d*)?)(?![\w$&])")
CHEMICALS_SEEN = {}
"""dict: formula to count, filled while pages are built (`--chemicals`)."""


def is_formula(text):
    """True if `text` reads as a chemical formula with at least one count: `H2O`, `CO2`, `Fe2O3`, not `GEN1377`."""
    if text in NOT_FORMULAS or not any(c.isdigit() for c in text):
        return False
    flat = text.replace("(", "").replace(")", "")
    pieces = _TOKEN.findall(flat)
    if "".join(sym + num for sym, num in pieces) != flat or any(sym not in ELEMENTS or len(num) > 2 for sym, num in pieces):
        return False
    if len(pieces) == 1:
        return text in MOLECULES
    return text in UPPERCASE_FORMULAS or any(c.islower() for c in text)


def subscript(text):
    """`text` (plain prose, already escaped) with each chemical formula's counts in `<sub>`."""
    def one(match):
        formula = match.group(1)
        if not is_formula(formula):
            return formula
        CHEMICALS_SEEN[formula] = CHEMICALS_SEEN.get(formula, 0) + 1
        return re.sub(r"(\d+)", r"<sub>\1</sub>", formula)
    return _FORMULA.sub(one, text)


SKIP_TAGS = {"code", "pre", "script", "style", "kbd", "samp", "textarea"}


def decorate(body):
    """Subscripts chemical formulas in the text of rendered HTML `body`, leaving code, links' addresses and math alone."""
    out, depth, math = [], 0, 0
    for part in re.split(r"(<[^>]+>)", body):
        if part.startswith("<") and part.endswith(">"):
            name = re.match(r"</?\s*([a-zA-Z0-9]+)", part)
            name = name.group(1).lower() if name else ""
            closing = part.startswith("</")
            if name in SKIP_TAGS and not part.endswith("/>"):
                depth += -1 if closing else 1
            if name == "span" and 'class="math' in part:
                math += 1
            elif name == "span" and math and closing:
                math -= 1
            out.append(part)
        else:
            out.append(part if depth > 0 or math > 0 else subscript(part))
    return "".join(out)


# --- math ------------------------------------------------------------------------------------------------------

_FENCE = re.compile(r"^(```|~~~)")
_DISPLAY = re.compile(r"\$\$(.+?)\$\$", re.S)
_INLINE = re.compile(r"(?<![\\$\w])\$(?!\s)([^$\n]{1,240}?)(?<![\s\\])\$(?![\w$])")
_CODE_SPAN = re.compile(r"(`+)(?:.|\n)*?\1")
_MATHY = re.compile(r"[\\^_{}]|^[A-Za-z]\w{0,2}$|^[\d.]+$")


def protect_math(text):
    """`(text with each formula replaced by a private marker, [(display?, tex)])`; code blocks and spans are untouched."""
    found = []

    def keep(display, tex):
        found.append((display, tex.strip()))
        return f"{len(found) - 1}"

    out, buffer, fenced = [], [], False
    for line in text.split("\n"):
        if _FENCE.match(line.lstrip()):
            if not fenced:
                out.append(_process_prose("\n".join(buffer), keep))
                buffer = []
            fenced = not fenced
            out.append(line)
            continue
        (out if fenced else buffer).append(line)
    out.append(_process_prose("\n".join(buffer), keep))
    return "\n".join(part for part in out if part is not None), found


def _process_prose(prose, keep):
    """Replaces the math in `prose` (a stretch of Markdown with no fenced code), skipping inline code spans."""
    spans, last, pieces = [m.span() for m in _CODE_SPAN.finditer(prose)], 0, []

    def convert(chunk):
        chunk = _DISPLAY.sub(lambda m: keep(True, m.group(1)), chunk)
        return _INLINE.sub(lambda m: keep(False, m.group(1)) if _MATHY.search(m.group(1)) else m.group(0), chunk)

    for start, end in spans:
        pieces.append(convert(prose[last:start]))
        pieces.append(prose[start:end])
        last = end
    pieces.append(convert(prose[last:]))
    return "".join(pieces)


def restore_math(body, found):
    """Puts the formulas back as KaTeX targets."""
    def one(match):
        display, tex = found[int(match.group(1))]
        escaped = html.escape(tex, quote=False)
        return (f'<span class="math display">\\[{escaped}\\]</span>' if display
                else f'<span class="math">\\({escaped}\\)</span>')
    return re.sub("(\\d+)", one, body)


# --- markdown ----------------------------------------------------------------------------------------------------

def slug(text, seen):
    """GitHub's heading anchor for `text`: lower case, no punctuation, spaces to hyphens, `-1` for repeats."""
    base = re.sub(r"[^\w\- ]", "", text.strip().lower()).replace(" ", "-")
    n = seen.get(base, 0)
    seen[base] = n + 1
    return base if n == 0 else f"{base}-{n}"


def make_parser():
    return MarkdownIt("commonmark", {"html": True}).enable(["table", "strikethrough"])


def inline_text(token):
    return "".join(child.content for child in (token.children or []) if child.type in ("text", "code_inline"))


class Page:
    """One Markdown file and where its page goes."""

    def __init__(self, source, out):
        self.source = source  # path relative to the repository root, with slashes
        self.out = out  # path relative to docs/site/, with slashes
        self.title = ""
        self.summary = ""
        self.headings = []  # (level, text, anchor)
        self.body = ""
        self.text = ""


def discover():
    """Every `Page` the site holds, ordered by path."""
    pages = []
    for folder, dirs, files in os.walk(DOCS):
        dirs[:] = sorted(d for d in dirs if d not in ("__pycache__", "site"))
        for name in sorted(files):
            if name.endswith(".md"):
                rel = os.path.relpath(os.path.join(folder, name), ROOT).replace(os.sep, "/")
                if any(part in rel for part in KEEP_NAMES):
                    continue  # these keep sources by name (DOC.99); links to them go to GitHub, the reference list to docs/REFERENCES.html
                pages.append(Page(rel, rel[len("docs/"):][:-3] + ".html"))
    for name in ROOT_FILES:
        if os.path.exists(os.path.join(ROOT, name)) and not any(part in name for part in KEEP_NAMES):
            pages.append(Page(name, "repo/" + name[:-3] + ".html"))
    return sorted(pages, key=lambda p: p.out)


def link_target(page, href, by_source):
    """Where `href`, written in `page`'s Markdown, points on the site."""
    parts = urlsplit(href)
    if parts.scheme or parts.netloc or href.startswith(("#", "mailto:")):
        return href
    path = unquote(parts.path)
    base = posixpath.dirname(page.source)
    target = posixpath.normpath(posixpath.join(base, path)) if path else page.source
    fragment = f"#{parts.fragment}" if parts.fragment else ""
    here = posixpath.dirname(page.out) or "."
    if target in by_source:
        return posixpath.relpath(by_source[target].out, here) + fragment
    if target.startswith("docs/") and os.path.isfile(os.path.join(ROOT, target)):
        return posixpath.relpath(target, posixpath.join("docs/site", here)) + fragment
    if os.path.isdir(os.path.join(ROOT, target)):
        return f"{GITHUB.replace('/blob/', '/tree/')}{target.rstrip('/')}"
    if os.path.isfile(os.path.join(ROOT, target)):
        return GITHUB + target + fragment
    return href


def render(page, by_source):
    """Fills in `page.title`, `.summary`, `.headings`, `.body` and `.text` from its Markdown file."""
    with open(os.path.join(ROOT, page.source), encoding="utf-8") as handle:
        source = handle.read()
    protected, math = protect_math(source)
    parser = make_parser()
    tokens = parser.parse(protected)
    seen, summary = {}, ""
    for i, token in enumerate(tokens):
        if token.type == "heading_open":
            text = inline_text(tokens[i + 1])
            anchor = slug(text, seen)
            token.attrSet("id", anchor)
            level = int(token.tag[1:])
            page.headings.append((level, text, anchor))
            if level == 1 and not page.title:
                page.title = text
        elif token.type == "paragraph_open" and not summary and page.title and not tokens[i + 1].content.startswith(("|", "<")):
            summary = tokens[i + 1].content
        elif token.type == "inline":
            for child in token.children or []:
                if child.type == "link_open":
                    child.attrSet("href", link_target(page, child.attrGet("href") or "", by_source))
                elif child.type == "image":
                    child.attrSet("src", link_target(page, child.attrGet("src") or "", by_source))
    body = restore_math(parser.renderer.render(tokens, parser.options, {}), math)
    body = body.replace("<table>", '<div class="table"><table>').replace("</table>", "</table></div>")
    page.body = decorate(body)
    page.title = page.title or os.path.basename(page.source)[:-3].replace("-", " ").capitalize()
    plain = re.sub(r"\s+", " ", re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", re.sub(r"[`*_]", "", summary))).strip()
    sentence = re.match(r"(.+?[.!?])(\s|$)", plain)
    first = sentence.group(1) if sentence else plain
    page.summary = first if len(first) <= 240 else first[:237] + "..."
    page.text = re.sub(r"\s+", " ", re.sub(r"[#*`|>\-]|\[([^\]]*)\]\([^)]*\)", " ", source)).strip()[:600]


# --- pages -------------------------------------------------------------------------------------------------------

STYLE = """
:root{--bg:#F3F5F8;--surface:#fff;--sunk:#E9EDF3;--ink:#152030;--muted:#566275;--line:#D8DEE7;--link:#1F5FA8;--code:#EEF1F6}
@media (prefers-color-scheme: dark){:root{--bg:#0C1118;--surface:#141B25;--sunk:#1A2330;--ink:#E5EAF1;--muted:#97A3B5;--line:#273244;--link:#7DB4F0;--code:#1B2431;color-scheme:dark}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
a{color:var(--link);text-decoration:none} a:hover{text-decoration:underline}
header.top{position:sticky;top:0;z-index:5;background:var(--surface);border-bottom:1px solid var(--line);display:flex;gap:12px;align-items:center;padding:8px 16px}
header.top .brand{font-weight:700;color:var(--ink)} header.top input{flex:1;max-width:420px;padding:6px 10px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--ink);font:inherit}
header.top .src{margin-left:auto;font-size:14px;color:var(--muted)} #menu-toggle{display:none}
.layout{display:flex;align-items:flex-start}
nav#nav{width:290px;flex:none;position:sticky;top:49px;max-height:calc(100vh - 49px);overflow:auto;padding:12px 10px 40px 16px;font-size:14px;border-right:1px solid var(--line)}
nav#nav details{margin:2px 0} nav#nav summary{cursor:pointer;font-weight:600;padding:2px 0} nav#nav ul{list-style:none;margin:0;padding-left:12px}
nav#nav li{padding:1px 0} nav#nav a.here{font-weight:700;color:var(--ink)} #results{list-style:none;padding:0;margin:0}
#results li{padding:3px 0}
main{flex:1;min-width:0;max-width:980px;padding:22px 28px 80px}
main h1{font-size:1.9rem;line-height:1.2;margin:0 0 .4em} main h2{font-size:1.4rem;margin:1.6em 0 .4em;border-bottom:1px solid var(--line);padding-bottom:3px}
main h3{font-size:1.15rem;margin:1.4em 0 .3em} main h4{margin:1.2em 0 .2em}
main :is(h1,h2,h3,h4,h5,h6){scroll-margin-top:60px}
code{background:var(--code);padding:1px 5px;border-radius:4px;font:.88em ui-monospace,Menlo,Consolas,monospace}
pre{background:var(--code);padding:12px 14px;border-radius:8px;overflow:auto} pre code{background:none;padding:0}
blockquote{margin:1em 0;padding:2px 14px;border-left:4px solid var(--line);color:var(--muted)}
.table{overflow-x:auto;margin:1em 0} table{border-collapse:collapse;font-size:14.5px} th,td{border:1px solid var(--line);padding:5px 9px;vertical-align:top;text-align:left}
th{background:var(--sunk)} img{max-width:100%} sub{font-size:.75em;line-height:0}
details.toc{background:var(--surface);border:1px solid var(--line);border-radius:8px;padding:6px 14px;margin:0 0 18px;font-size:14.5px}
details.toc ul{list-style:none;padding-left:14px;margin:4px 0} details.toc>ul{padding-left:0}
.math.display{display:block;overflow-x:auto;margin:.8em 0}
.meta{color:var(--muted);font-size:13px;margin:0 0 14px}
.cards{list-style:none;padding:0} .cards li{background:var(--surface);border:1px solid var(--line);border-radius:8px;padding:8px 13px;margin:7px 0}
.cards li p{margin:2px 0 0;color:var(--muted);font-size:14.5px}
@media (max-width:860px){nav#nav{display:none;position:fixed;top:49px;left:0;z-index:4;background:var(--surface);width:84vw;max-height:calc(100vh - 49px)}
nav#nav.open{display:block} #menu-toggle{display:inline-block;background:none;border:1px solid var(--line);border-radius:6px;color:var(--ink);padding:4px 9px}
main{padding:16px 16px 60px}}
"""

SCRIPT = """
(function(){
var nav=document.getElementById('nav'),root=nav.dataset.root,here=nav.dataset.here;
function tree(entries){var html='',groups={};entries.forEach(function(e){(groups[e.section]=groups[e.section]||[]).push(e);});
 Object.keys(groups).forEach(function(s){var open=groups[s].some(function(e){return e.path===here;});
  html+='<details'+(open?' open':'')+'><summary>'+s+'</summary><ul>'+groups[s].map(function(e){
   return '<li><a href="'+root+e.path+'"'+(e.path===here?' class="here"':'')+'>'+e.title.replace(/&/g,'&amp;').replace(/</g,'&lt;')+'</a></li>';}).join('')+'</ul></details>';});
 return html;}
var data=null,box=document.getElementById('q');
fetch(root+'nav.json').then(function(r){return r.json();}).then(function(d){data=d;nav.innerHTML=tree(d);}).catch(function(){nav.innerHTML='<a href="'+root+'index.html">Site index</a>';});
box.addEventListener('input',function(){var q=box.value.trim().toLowerCase();if(!data)return;
 if(!q){nav.innerHTML=tree(data);return;}
 var words=q.split(/\\s+/),hits=data.filter(function(e){var h=(e.title+' '+e.path+' '+e.text).toLowerCase();return words.every(function(w){return h.indexOf(w)>=0;});}).slice(0,60);
 nav.innerHTML='<ul id="results">'+hits.map(function(e){return '<li><a href="'+root+e.path+'">'+e.title.replace(/&/g,'&amp;').replace(/</g,'&lt;')+'</a><br><small>'+e.path+'</small></li>';}).join('')+'</ul>';});
var t=document.getElementById('menu-toggle');t.addEventListener('click',function(){nav.classList.toggle('open');});
})();
"""

KATEX_TAGS = f"""<link rel="stylesheet" href="{KATEX}/katex.min.css">
<script defer src="{KATEX}/katex.min.js"></script>
<script defer src="{KATEX}/contrib/auto-render.min.js" onload="renderMathInElement(document.body,{{delimiters:[{{left:'\\\\[',right:'\\\\]',display:true}},{{left:'\\\\(',right:'\\\\)',display:false}}],throwOnError:false}})"></script>"""


def shell(title, description, root, here, body, source=None):
    """A whole page. `root` is the relative path from the page to `docs/site/`; `here` its path inside the site."""
    edit = (f'<a class="src" href="{GITHUB}{html.escape(source)}">Source on GitHub</a>' if source else "")
    uses_math = 'class="math' in body
    return f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title, quote=False)} - planetGen docs</title>
<meta name="description" content="{html.escape(description, quote=True)}">
<style>{STYLE}</style>
{KATEX_TAGS if uses_math else ""}
<header class="top"><button id="menu-toggle" aria-label="Menu">Menu</button><a class="brand" href="{root}index.html">planetGen docs</a>
<input id="q" type="search" placeholder="Search the docs" aria-label="Search the docs">{edit}</header>
<div class="layout">
<nav id="nav" data-root="{root}" data-here="{html.escape(here)}"><a href="{root}index.html">Site index</a></nav>
<main>
{body}
</main>
</div>
<script>{SCRIPT}</script>
"""


def toc(headings):
    """A collapsed table of contents of the `##` and `###` headings, or empty when there are fewer than four."""
    items = [(lvl, text, anchor) for lvl, text, anchor in headings if lvl in (2, 3)]
    if len(items) < 4:
        return ""
    lines = [f'<li style="margin-left:{(lvl - 2) * 14}px"><a href="#{anchor}">{html.escape(text, quote=False)}</a></li>'
             for lvl, text, anchor in items[:200]]
    return f'<details class="toc"><summary>On this page ({len(items)} headings)</summary><ul>{"".join(lines)}</ul></details>'


def section_of(page):
    first = page.out.split("/")[0] if "/" in page.out else ""
    return SECTION_TITLES.get(first, first.capitalize())


def build_files():
    """`{path under docs/site/: text}` of the whole site."""
    CHEMICALS_SEEN.clear()
    pages = discover()
    by_source = {p.source: p for p in pages}
    for page in pages:
        render(page, by_source)
    files = {}
    for page in pages:
        depth = page.out.count("/")
        root = "../" * depth
        meta = f'<p class="meta">{html.escape(page.source)}</p>'
        body = meta + toc(page.headings) + page.body
        # the first heading is the page's own; the meta line goes under it
        if body.find("</h1>") >= 0:
            head, _, rest = page.body.partition("</h1>")
            body = head + "</h1>" + meta + toc(page.headings) + rest
        files[page.out] = shell(page.title, page.summary or page.title, root or "./", page.out, body, page.source)
    entries = [{"path": p.out, "title": p.title, "section": section_of(p), "text": p.summary + " " + p.text[:200]} for p in pages]
    for name in sorted(os.listdir(os.path.join(DOCS, "generated"))):
        if name.endswith(".html"):
            entries.append({"path": "../generated/" + name, "title": name[:-5].replace("-", " ").capitalize() + " (generated page)",
                            "section": "Generated pages", "text": ""})
    entries.append({"path": "../REFERENCES.html", "title": "References (papers, datasets, catalogs)", "section": "Reference papers and datasets", "text": ""})
    order = list(SECTION_TITLES.values()) + ["Generated pages"]
    entries.sort(key=lambda e: (order.index(e["section"]) if e["section"] in order else 99, e["path"]))
    files["nav.json"] = json.dumps(entries, ensure_ascii=False, separators=(",", ":")) + "\n"
    files["index.html"] = front_page(entries)
    return files


def front_page(entries):
    sections = {}
    for e in entries:
        sections.setdefault(e["section"], []).append(e)
    parts = ["<h1>planetGen documentation</h1>",
             '<p class="meta">Built by scripts/build_docs_site.py from the Markdown files in the repository. Search with the box at the top.</p>']
    for section, items in sections.items():
        parts.append(f"<h2>{html.escape(section)}</h2><ul class=\"cards\">")
        for e in items:
            parts.append(f'<li><a href="{html.escape(e["path"])}">{html.escape(e["title"], quote=False)}</a> '
                         f'<small class="meta">{html.escape(e["path"])}</small></li>')
        parts.append("</ul>")
    return shell("Documentation site", "Every Markdown document of the planetGen repository as a web page, with a menu and search.",
                 "./", "index.html", "\n".join(parts))


def existing():
    found = set()
    for folder, _dirs, files in os.walk(SITE):
        for name in files:
            found.add(os.path.relpath(os.path.join(folder, name), SITE).replace(os.sep, "/"))
    return found


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="fail if docs/site/ is out of date")
    parser.add_argument("--chemicals", action="store_true", help="list the chemical formulas that got subscripts")
    args = parser.parse_args(argv)
    files = build_files()
    if args.chemicals:
        for formula, count in sorted(CHEMICALS_SEEN.items()):
            print(f"{count:6d}  {formula}")
        return 0
    have = existing() if os.path.isdir(SITE) else set()
    stale = sorted(have - set(files))
    changed = []
    for path, text in sorted(files.items()):
        full = os.path.join(SITE, path)
        current = open(full, encoding="utf-8").read() if os.path.exists(full) else None
        if current != text:
            changed.append(path)
            if not args.check:
                os.makedirs(os.path.dirname(full), exist_ok=True)
                with open(full, "w", encoding="utf-8", newline="\n") as handle:
                    handle.write(text)
    if not args.check:
        for path in stale:
            os.remove(os.path.join(SITE, path))
        print(f"docs/site: {len(files)} pages, {len(changed)} written, {len(stale)} removed")
        return 0
    if changed or stale:
        print(f"docs/site is out of date ({len(changed)} changed, {len(stale)} stale): run python scripts/build_docs_site.py")
        return 1
    print("docs site is up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
