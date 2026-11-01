#!/usr/bin/env python3
"""
Moves every citation out of the repository's text into the central registry (DOC.98): an author-year mention, a DOI and an arXiv number
become the source's citation ID in square brackets, `[R0123]` (the numbered style of NASA documents; the ID is the number).

    python scripts/rewrite_citations.py report [--scope SCOPE ...]   # what would change, by folder and kind; writes nothing
    python scripts/rewrite_citations.py apply  [--scope SCOPE ...]   # rewrite the files in place
    python scripts/rewrite_citations.py check  [--scope SCOPE ...]   # exit 1 when a citation is still written out

A scope is a git path prefix (`docs/design`, `src/planetgen/physics`). With no scope the whole repository is read. What is rewritten:

* A text file (`.md`, `.txt`, `.rst`, `.yml`, `.tsv`): the whole file.
* A Python file: its comments and docstrings, always; its other string literals only with `--strings` (a string can reach a page or a
  stored record, so the output is checked first). After rewriting, the syntax tree of the file, docstrings left out, must equal the one
  before, or the file is not written.
* A JavaScript file: its `//` and `/* */` comments, always; strings only with `--strings`.
* A data or markup file (`.json`, `.csv`, `.html`, `.sh`): only with `--strings`.

A mention is replaced only when the registry ties it to exactly one source: the last settle step (`scripts/settle_references.py`) recorded the
paper, or the volume and page, the identifier or the letter beside the mention (`2023a`) pick one of several. A mention the registry calls
"not a paper" (a time stamp, an organization, a standard) or "for review" is left alone, and so is one that fits two sources and nothing
beside it chooses: those are listed by `report` as ambiguous. Nothing is guessed.

The files that must keep names are in `KEEP` (the registry itself and the build tools, the data license, the changelog history, third-party source
texts). Standard library only, apart from `build_references`.
"""

import argparse
import ast
import io
import json
import os
import re
import subprocess
import sys
import tokenize
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import build_references as br  # noqa: E402

TEXT_FILES = (".md", ".txt", ".rst", ".yml", ".yaml", ".tsv", ".toml")
DATA_FILES = (".json", ".csv", ".html", ".sh")
PY_FILES = (".py",)
JS_FILES = (".js", ".mjs")

KEEP = (
    "docs/REFERENCES.", "docs/reference/references-data.json", "docs/references/",
    "scripts/build_references.py", "scripts/settle_references.py", "scripts/references_curated.py", "scripts/rewrite_citations.py",
    "DATA-LICENSE.md", "LICENSE.md", "CHANGELOG.md",
    "src/tests/test_references.py", "src/tests/test_citations_links.py", "src/tests/test_references_crawl.py", "src/tests/test_references_pages.py", "src/tests/test_citations_outside_the_registry.py",
    "src/planetgen/datasets/catalogue.py",
    "node_modules/", "vendor/", "three.min", ".min.js", "package-lock",
)
"""tuple: Path fragments of the files that keep author names, DOIs and arXiv numbers (the registry and its builders, the data license, the changelog
history, the download catalogue that needs its addresses)."""

NOT_FILES = ("src/.test-reports/",)
"""tuple: Path fragments of generated reports no one edits."""

JOURNALS = ("ApJ", "ApJS", "ApJL", "AJ", "A&A", "A&AS", "MNRAS", "Icarus", "Nature", "Science", "PASP", "ARA&A", "GeCoA", "GCA", "JGR", "GRL", "EPSL", "PSJ", "P&SS",
            "Astrobiology", "Space Sci. Rev.", "SSRv", "Phys. Rev. D", "PhRvD", "PhRvL", "Rev. Mod. Phys.", "RvMP", "A&ARv", "Sol. Phys.", "CeMDA", "NatAs", "NatGe",
            "Natur", "Sci", "Meteoritics", "Meteorit. Planet. Sci.", "MAPS", "JCP", "J. Chem. Phys.", "J. Phys. Chem. Ref. Data", "JPCRD")
"""tuple: Journal abbreviations whose volume and page that follow a replaced citation (`ApJ 758, 29`) are dropped with it."""

NAME = br.NAME
PART = r"(?:(?:Di|De|Del|Della|Van|Von|van|von|de|der|den|du|La|Le|Mac|Mc|St\.|ter|ten|el|al)\s+)"
SEP = r"(?:,\s+and\s+|\s+and\s+|,\s+&\s+|\s+&\s+|,\s+)"
BEFORE = re.compile(r"(?:(?:" + PART + r"?" + NAME + r")" + SEP + r"){1,3}\Z|" + PART + r"\Z|(?:(?:" + PART + r"?" + NAME + r")" + SEP + r"){1,3}" + PART + r"\Z")
"""re.Pattern: The names that come before the author `br.AUTHYEAR` found, which belong to the same author list: `Moe, Kratter` before `and Badenes 2019`
(the harvest reads only the last two names), `Moe and Di` before `Stefano 2017`."""

LETTERS = re.compile(r"(?:\s*[,/&]\s*[a-z]\b)*")
JOURNAL = r"(?:" + "|".join(re.escape(j) for j in sorted(JOURNALS, key=len, reverse=True)) + r")"
TAIL = re.compile(r"(?P<sep>,?\s+)" + JOURNAL + r"\s+\d{1,4}(?:,\s*[A-Za-z]?\d{1,6})?\b")
JOURNAL_PART = re.compile(JOURNAL + r"\s+\d{1,4}(?:,\s*[A-Za-z]?\d{1,6})?\b")
PLACE = re.compile(r"(?:Table|Tables|Fig\.|Figs\.|Figure|Section|Sections|Sec\.|Eq\.|Appendix|Annex|p\.|pp\.)\s+[A-Z]?\d+(?:\.\d+)*[a-z]?(?:\s*(?:to|and|-)\s*\d+)?")

DOI_LINK = re.compile(r"\[[^\]\n]*\]\((?:https?://(?:dx\.)?doi\.org/(?P<d>[^)\s]+)|https?://arxiv\.org/(?:abs|pdf)/(?P<x>[^)\s]+?)(?:v\d+)?(?:\.pdf)?)\)", re.I)
DOI_BARE = re.compile(r"(?:https?://(?:dx\.)?doi\.org/|\bdoi:\s*|\bDOI:?\s+)(?P<d>10\.\d{4,9}/(?:[^\s\"'`<>,;|\\\[\]{}()]|\([^\s()]*\))+)", re.I)
DOI_NAKED = re.compile(r"(?<![\w/.])(?P<d>10\.\d{4,9}/(?:[^\s\"'`<>,;|\\\[\]{}()]|\([^\s()]*\))+)")
ARXIV_ID = re.compile(
    r"https?://(?:www\.|export\.|ar5iv\.(?:labs\.)?)?arxiv\.org/(?:abs|pdf|html)/(?P<x>" + br.RE_ARXIV_NEW + "|" + br.RE_ARXIV_OLD + r")(?:v\d+)?(?:\.pdf)?"
    r"|\barXiv[:\s.]*(?:abs/)?(?P<x2>" + br.RE_ARXIV_NEW + "|" + br.RE_ARXIV_OLD + r")(?:v\d+)?", re.I)

ID = r"R\d{4,}"
STOP = set("""about their these those which where while there other using their under between through could would should around within without
across study results first second third stars star planet planets using based model models""".split())
BRACKET = re.compile(r"\[(" + ID + r")\]")
GROUP = re.compile(r"\[(" + ID + r")\](?:\s*(?:;|,|\band\b|&)\s*\[(" + ID + r")\])+")


class Registry:
    """The registry as the rewriter needs it: which source an author-year mention names, and the ID of a DOI or an arXiv number."""

    def __init__(self, cache=None, all_candidates=False):
        self.all_candidates = all_candidates
        self.c = cache or br.load_cache()
        br.assign_ids(self.c)
        self.ids = self.c["ids"]
        self.status = self.c["authoryear_status"]
        self.ents = self.c["entries"]
        papers = {k: e for k, e in self.ents.items() if (e.get("authors") or e.get("title")) and e.get("type") != "dataset"}
        datas = {k: e for k, e in self.ents.items() if (e.get("authors") or e.get("title")) and e.get("type") == "dataset"}
        self.suffix = {}
        self.suffix.update(br.sorted_entries(papers)[1])
        self.suffix.update(br.sorted_entries(datas)[1])

    def id_of(self, key):
        return self.ids[key.lower()]

    def key_of(self, m):
        a, b = m.group("a"), m.group("b")
        kind = "pair" if b else ("etal" if "et al" in m.group(0) else "one")
        return "|".join((a, b or "", kind, m.group("y")))

    def pick(self, key, letters, after, before=""):
        """The one source ID a mention names, or `(None, why)`: tied to one source; or to several and the letter, a DOI or arXiv number, or the
        volume and page after it pick one."""
        st = self.status.get(key, {})
        matched = [m for m in st.get("matched", []) if m in self.ents]
        if not matched:
            return None, "not a paper" if st.get("not_a_paper") else ("for review" if st.get("review") else "unknown")
        ids = sorted({self.id_of(m) for m in matched})
        if len(ids) == 1:
            return ids[0], ""
        if letters:
            want = [x.strip() for x in re.split(r"[,/&]", letters) if x.strip()]
            hit = sorted({self.id_of(m) for m in matched if self.suffix.get(m) in want})
            if len(hit) == 1:
                return hit[0], ""
        head = after[:260]
        found = []
        for m in matched:
            e = self.ents[m]
            vol, page = str(e.get("volume") or ""), str(e.get("page") or "")
            if (e.get("doi") and e["doi"].lower() in head.lower()) or (e.get("arxiv") and e["arxiv"].lower() in head.lower()):
                found.append(m)
            elif vol and page and re.search(r"\b" + re.escape(vol) + r",?\s+" + re.escape(page) + r"\b", head):
                found.append(m)
        hit = sorted({self.id_of(m) for m in found})
        if len(hit) == 1:
            return hit[0], ""
        # title words of one candidate (not of the others) in the sentence around the mention
        window = (before + " " + after[:400]).lower()
        score = {}
        for m in matched:
            words = {w for w in re.findall(r"[a-z]{5,}", (self.ents[m].get("title") or "").lower()) if w not in STOP}
            score[self.id_of(m)] = len({w for w in words if w in window})
        ranked = sorted(score.items(), key=lambda kv: -kv[1])
        if ranked[0][1] >= 1 and ranked[0][1] > ranked[1][1]:
            return ranked[0][0], ""
        return None, "ambiguous (%s)" % ", ".join(ids)


def join_group(ids):
    return "[" + ", ".join(sorted(set(ids))) + "]"


def tidy(text):
    """Merges neighbouring brackets, unwraps a parenthesis that holds only a citation, and moves a place inside the source into the bracket."""
    text = GROUP.sub(lambda m: join_group(re.findall(ID, m.group(0))), text)
    text = re.sub(r"\(\s*(\[" + ID + r"(?:, " + ID + r")*\])\s*\)", r"\1", text)
    text = re.sub(r"\[(" + ID + r")\],?\s+(Table|Tables|Fig\.|Figs\.|Figure|Section|Sec\.|Eq\.|Appendix|Annex|p\.|pp\.)\s+([A-Z]?\d+(?:\.\d+)*[a-z]?(?:\s*(?:to|and|-)\s*\d+)?)(?![\w.]\w)", r"[\1: \2 \3]", text)
    # the same source twice in a row
    text = re.sub(r"\[(" + ID + r")\](\s*[,;]\s*|\s+)\[\1\]", r"[\1]", text)
    return text


def biblio_only(seg):
    """`(ok, extra IDs, place)` for the inside of a parenthesis: ok when it holds nothing but journal, volume, page, a place in the source, other
    citation IDs and words like Letters, so the whole parenthesis goes with the citation it follows."""
    rest = seg
    extra = re.findall(r"\[(" + ID + r")\]", rest)
    rest = re.sub(r"\[" + ID + r"\]", "", rest)
    rest = JOURNAL_PART.sub("", rest)
    rest = re.sub(JOURNAL + r"(?:\s+Letters?)?\b", "", rest)
    places = PLACE.findall(rest)
    rest = PLACE.sub("", rest)
    rest = re.sub(r"\b(?:Letters?|Erratum|in press|submitted|preprint|arXiv|doi|published|accepted)\b", "", rest, flags=re.I)
    return bool(re.fullmatch(r"[\s,;.:]*", rest)), extra, (places[0] if len(places) == 1 else "")


def after_year(text, pos, opened, sid):
    """
    What follows the year of a mention. Returns `(end, place, extra IDs, reopen)`: where the mention ends, a place inside the source (`Table 4`), the other
    IDs named in the same parenthesis, and whether the parenthesis stays open because words follow that are not bibliography ("Name (2019, found that ...)").
    """
    if opened and text[pos:pos + 1] == ")":
        return pos + 1, "", [], False
    start = pos
    if not opened:
        lead = re.match(r"\s*\(", text[pos:])
        if not lead:
            return pos, "", [], False
        start = pos + lead.end()
    close = text.find(")", start)
    seg = text[start:close] if 0 < close - start <= 200 else None
    if seg is not None and "(" not in seg and "\n\n" not in seg:
        ok, extra, place = biblio_only(seg)
        if ok:
            return close + 1, place, extra, False
    if not opened:
        return pos, "", [], False
    skip = re.match(r"[,;]?\s*", text[pos:])
    return pos + (skip.end() if skip else 0), "", [], True


def rewrite_text(text, reg, stats=None, where=""):
    """`text` with each tied citation replaced; `stats` (a Counter) counts the outcomes; `reg` is a `Registry`."""
    stats = stats if stats is not None else Counter()

    def ref_of(doi=None, arxiv=None):
        names = []
        if doi:
            d = doi.rstrip(".,;:").lower()
            while d and d[-1] in ")]" and d.count("(") < d.count(")"):
                d = d[:-1]
            names.append(d)
        if arxiv:
            names.append("arxiv:" + re.sub(r"v\d+$", "", arxiv.lower()))
        for n in names:
            if n and n in reg.ids:
                return reg.ids[n]
        return None

    def links(m):
        sid = ref_of(m.group("d"), m.group("x"))
        if not sid:
            stats["identifier not in the registry"] += 1
            return m.group(0)
        stats["identifier"] += 1
        return "[%s]" % sid

    def bare(m):
        sid = ref_of(doi=m.group("d"))
        if not sid:
            stats["identifier not in the registry"] += 1
            return m.group(0)
        stats["identifier"] += 1
        return "[%s]" % sid

    def arx(m):
        sid = ref_of(arxiv=m.group("x") or m.group("x2"))
        if not sid:
            stats["identifier not in the registry"] += 1
            return m.group(0)
        stats["identifier"] += 1
        return "[%s]" % sid

    text = DOI_LINK.sub(links, text)
    text = DOI_BARE.sub(bare, text)
    text = DOI_NAKED.sub(bare, text)
    text = ARXIV_ID.sub(arx, text)

    out, pos = [], 0
    for m in br.AUTHYEAR.finditer(text):
        a, b, y = m.group(1), m.group(2), m.group(3)
        if m.start() < pos or a in br.NOT_NAMES or (b and b in br.NOT_NAMES) or not (1900 <= int(y) <= 2026):
            continue
        kind = "pair" if b else ("etal" if "et al" in m.group(0) else "one")
        key = "|".join((a, b or "", kind, y))
        end = m.end()
        letters = m.group(4) or ""
        more = LETTERS.match(text, end)
        if more and more.end() > end:
            letters += more.group(0)
            end = more.end()
        opened = "(" in m.group(0)[:m.group(0).rindex(y)]
        sid, why = reg.pick(key, letters, text[end:end + 300], text[max(0, m.start() - 200):m.start()])
        group = []
        if not sid and why.startswith("ambiguous") and reg.all_candidates:
            group = re.findall(ID, why)
            sid = group[0]
            stats["cited every candidate (same author and year, nothing beside it chooses)"] += 1
            stats["ambiguous: " + key] += 1
        if not sid:
            stats[why.split(" (")[0]] += 1
            if why.startswith("ambiguous"):
                stats["ambiguous: " + key] += 1
            continue
        start = m.start()
        before = BEFORE.search(text, max(pos, start - 160), start)
        if before and before.start() >= pos:
            start = before.start()
        out.append(text[pos:start])
        new_end, place, extra, reopen = after_year(text, end, opened, sid)
        ids = sorted({sid, *extra, *group})
        out.append("[%s%s]" % (", ".join(ids), (": " + place) if place and len(ids) == 1 else "") + (" (" if reopen else ""))
        pos = new_end
        stats["author-year"] += 1
        if not opened:
            tail = TAIL.match(text, pos)
            if tail:
                pos = tail.end()
                stats["journal tail"] += 1
    out.append(text[pos:])
    text = tidy("".join(out))
    return tidy(text)


# ---------------------------------------------------------------------------------------------------------- per file kind

def python_spans(text):
    """`[(start, end, kind)]` for each comment ("comment"), docstring ("doc") and other string ("string") of a Python file, by character offset."""
    lines = text.split("\n")
    offs = [0]
    for ln in lines:
        offs.append(offs[-1] + len(ln) + 1)
    docs = set()
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return None
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                docs.add((first.value.lineno, first.value.col_offset))
    spans = []
    for tk in tokenize.generate_tokens(io.StringIO(text).readline):
        if tk.type == tokenize.COMMENT:
            kind = "comment"
        elif tk.type == tokenize.STRING:
            kind = "doc" if tk.start in docs else "string"
        else:
            continue
        spans.append((offs[tk.start[0] - 1] + tk.start[1], offs[tk.end[0] - 1] + tk.end[1], kind))
    return spans


def doc_less_dump(text):
    tree = ast.parse(text)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                node.body[0] = ast.Pass()
    return ast.dump(tree)


def rewrite_python(text, reg, stats, strings=False):
    spans = python_spans(text)
    if spans is None:
        stats["python that does not parse"] += 1
        return text
    out, pos = [], 0
    for s, e, kind in spans:
        if kind == "string" and not strings:
            continue
        piece = text[s:e]
        new = rewrite_text(piece, reg, stats) if kind != "string" else rewrite_string_literal(piece, reg, stats)
        if new != piece:
            out.append(text[pos:s])
            out.append(new)
            pos = e
            stats["python " + kind] += 1
    out.append(text[pos:])
    new = "".join(out)
    if new != text and not strings and doc_less_dump(new) != doc_less_dump(text):
        stats["python changed in more than its comments and docstrings"] += 1
        return text
    if new != text:
        try:
            ast.parse(new)
        except SyntaxError:
            stats["python that no longer parses"] += 1
            return text
    return new


def rewrite_string_literal(piece, reg, stats):
    return rewrite_text(piece, reg, stats)


JS_COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.S)


def rewrite_js(text, reg, stats, strings=False):
    if strings:
        return rewrite_text(text, reg, stats)
    return JS_COMMENT.sub(lambda m: rewrite_text(m.group(0), reg, stats), text)


def rewrite_file(path, text, reg, stats, strings=False):
    if path.endswith(PY_FILES):
        return rewrite_python(text, reg, stats, strings)
    if path.endswith(JS_FILES):
        return rewrite_js(text, reg, stats, strings)
    if path.endswith(TEXT_FILES):
        return rewrite_text(text, reg, stats)
    if path.endswith(DATA_FILES) and strings:
        return rewrite_text(text, reg, stats)
    return text


def tracked(scopes):
    names = subprocess.check_output(["git", "ls-files"], cwd=ROOT, text=True).split("\n")
    for n in names:
        if not n or not n.endswith(TEXT_FILES + DATA_FILES + PY_FILES + JS_FILES):
            continue
        if any(k in n for k in KEEP + NOT_FILES) or (scopes and not any(n.startswith(s) for s in scopes)):
            continue
        yield n


def cmd(args):
    reg = Registry(all_candidates=args.all_candidates)
    total = Counter()
    by_dir = defaultdict(Counter)
    ambiguous = Counter()
    changed = 0
    for n in tracked(args.scope):
        path = os.path.join(ROOT, n)
        try:
            with open(path, encoding="utf-8") as fh:
                text = fh.read()
        except (OSError, UnicodeDecodeError):
            continue
        stats = Counter()
        new = rewrite_file(n, text, reg, stats, args.strings)
        for k, v in stats.items():
            if k.startswith("ambiguous: "):
                ambiguous[k[11:]] += v
            else:
                total[k] += v
                by_dir["/".join(n.split("/")[:3] if n.startswith("src/planetgen") else n.split("/")[:2])][k] += v
        if new != text:
            changed += 1
            if args.cmd == "apply":
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(new)
    print("files %s: %d" % ("rewritten" if args.cmd == "apply" else "that would change", changed))
    for k, v in sorted(total.items()):
        print("  %-60s %6d" % (k, v))
    if args.verbose:
        for d, c in sorted(by_dir.items()):
            print(d, dict(c))
        print("ambiguous (key: places):")
        for k, v in ambiguous.most_common(args.verbose):
            print("  ", k, v)
    if args.cmd == "check" and (total["author-year"] or total["identifier"] or ambiguous):
        sys.exit(1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("report", "apply", "check"))
    ap.add_argument("--scope", nargs="*", default=[])
    ap.add_argument("--all-candidates", action="store_true", help="where nothing chooses between two papers of one author and year, cite both")
    ap.add_argument("--strings", action="store_true", help="also rewrite string literals and data files")
    ap.add_argument("-v", "--verbose", nargs="?", const=30, type=int, default=0)
    cmd(ap.parse_args())
