#!/usr/bin/env python3
"""DOC.95 (References R3b): tie every author-year citation in the repository's text to its paper, or label it as not a paper.

Works on docs/reference/references-data.json after `build_references.py harvest` has run. For each author-year mention that is not yet
tied to a paper, in this order:

  1. triage    - a time stamp, a month, a sentence word or an organisation, standard, web page or program is labelled "not a paper"
                 with what it is and a link (the labels are in references_curated.py);
  2. local     - a full citation elsewhere in the text (an arXiv number or DOI beside the mention, or an entry already resolved) is
                 matched on first author, year and (for pairs) second author, at any author position;
  3. crossref  - the rest are searched by author, year and the citing sentences; a record is accepted only when the first author and
                 year agree AND either the volume and first page appear beside the mention, or two title words appear in the citing
                 sentences, or (pairs) both authors agree with the year exact and one title word or the journal appears. A tie or a
                 near miss goes to the review list with the best candidates and why it was not taken. Nothing is guessed.

Every decision is written to the cache (authoryear_status[key]) so a later harvest does not lose it; the REFERENCES pages and the test
in src/tests/test_references.py read it from there.

    python scripts/settle_references.py triage
    python scripts/settle_references.py local
    python scripts/settle_references.py search [--limit N] [--min-places N]
    python scripts/settle_references.py report
"""
import argparse
import json
import math
import os
import re
import sys
import urllib.parse
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_references as br  # noqa: E402
import references_curated as cur  # noqa: E402


def split_key(key):
    a, b, kind, y = key.split("|")
    return a, b, kind, int(y)


def fam_match(family, name):
    """A cited name matches an author's family name when they agree after folding, also for 'Di Stefano' cited as 'Stefano'."""
    f, n = (br.fold(re.sub(r"(?i)ae|oe|ue", lambda m: m.group(0)[0], t)) for t in (family, name))
    return bool(n) and len(n) >= 3 and (f == n or f.endswith(n) or (len(f) >= 5 and n.startswith(f)) or (len(n) >= 5 and f.startswith(n)))


def snippets(c, key):
    v = c["candidates"]["authoryear"].get(key) or {}
    s = list(v.get("ctxs") or [])
    if not s and v.get("ctx"):
        s = [v["ctx"]]
    st = c.get("authoryear_status", {}).get(key, {})
    return s or st.get("ctxs") or []


def status_of(c, key):
    return c.setdefault("authoryear_status", {}).setdefault(key, {"n": c["candidates"]["authoryear"].get(key, {}).get("n", 0), "matched": []})


def tied(v):
    return bool(v.get("matched"))


def settled(v):
    return tied(v) or bool(v.get("not_a_paper"))


# ---------------------------------------------------------------------------------------------------------------- triage
def cmd_triage(c):
    n_noise = n_np = 0
    for key in list(c["candidates"]["authoryear"]):
        v = status_of(c, key)
        if settled(v):
            continue
        a, b, kind, y = split_key(key)
        first = a
        if first in cur.NOISE:
            v["not_a_paper"] = {"what": cur.NOISE[first], "url": "", "class": "not a citation"}
            n_noise += 1
        elif first in cur.NOT_PAPERS:
            what, url = cur.NOT_PAPERS[first]
            v["not_a_paper"] = {"what": what, "url": url, "class": "organisation, standard, web page or program"}
            n_np += 1
    br.log("labelled not a citation:", n_noise, "labelled organisation, standard, page or program:", n_np)


# ---------------------------------------------------------------------------------------------------------------- local
def entry_index(c):
    ids = {}
    for k, e in c["entries"].items():
        if e.get("doi"):
            ids[e["doi"].lower()] = k
        if e.get("arxiv"):
            ids[e["arxiv"].lower()] = k
        ids[k.lower()] = k
    return ids


def authors_fit(e, a, b, kind, year):
    """True when entry e carries author a (and, for a pair, b right after it) with a year within one of the cited year."""
    ey = e.get("year")
    if not ey or abs(int(ey) - year) > 1:
        return False
    au = e.get("authors") or []
    if kind == "pair":
        return any(fam_match(au[i][0], a) and fam_match(au[i + 1][0], b) for i in range(min(len(au) - 1, 4)))
    return bool(au) and fam_match(au[0][0], a)


def ids_in(text):
    out = []
    for m in br.ARXIV.finditer(text):
        out.append(re.sub(r"v\d+$", "", m.group(1)).lower())
    for m in br.ARXIV_LIST.finditer(text):
        out += [i.lower() for i in re.findall(br.RE_ARXIV_NEW, m.group(1))]
    for m in br.DOI.finditer(text):
        d = br.clean_doi(m.group(1))
        if d:
            out.append(d.lower())
    return out


BARE_ARXIV = re.compile(r"(?<![\d.])((?:0[7-9]|1\d|2[0-6])(?:0[1-9]|1[0-2])\.\d{4,5}|[a-z]+(?:-[a-z]+)?/\d{7})(?:v\d+)?(?![\d])")


def cmd_bare(c):
    """An arXiv number written beside a mention without the word 'arXiv' (for example 'Wipperfurth 2020 (1907.12184)'): fetch the paper
    and tie the mention when the cited name is one of its first authors and the year is within two."""
    want = {}
    for key in c["candidates"]["authoryear"]:
        v = status_of(c, key)
        if tied(v) or v.get("not_a_paper"):
            continue
        for sn in snippets(c, key):
            for m in BARE_ARXIV.finditer(sn):
                want.setdefault(key, set()).add(m.group(1))
    ids = sorted({i for s_ in want.values() for i in s_})
    have = {e["arxiv"]: k for k, e in c["entries"].items() if e.get("arxiv")}
    got = br.fetch_arxiv([i for i in ids if i not in have])
    br.log("bare arXiv numbers beside unsettled mentions:", len(ids), "fetched:", len(got))
    n = 0
    for key, idset in want.items():
        a, b, kind, y = split_key(key)
        for i in sorted(idset):
            e = None
            if i in have:
                e = c["entries"][have[i]]
                k = have[i]
            elif i in got:
                g = got[i]
                e = {"authors": g["authors"], "year": g["year"], "title": g["title"], "arxiv": i, "doi": g["doi"], "journal": "", "volume": "", "page": "",
                     "type": "preprint", "src": "arxiv", "journal_ref": g["journal_ref"]}
                k = "arXiv:" + i
            if not e or not e.get("year") or not (y - 2 <= int(e["year"]) <= y + 1):
                continue
            au = (e.get("authors") or [])[:4]
            if not any(fam_match(x[0], a) for x in au) or (kind == "pair" and not any(fam_match(x[0], b) for x in au)):
                continue
            c["entries"].setdefault(k, e)
            v = status_of(c, key)
            v["matched"] = [k]
            v["evidence"] = "an arXiv number beside the mention (%s) belongs to a paper whose authors and year fit" % i
            v.pop("review", None)
            n += 1
            break
    br.log("tied by a bare arXiv number:", n)


FULLCITE = re.compile(r"^\W{0,6}([A-Z][^\W\d_][\w'\u2019\-]*),\s+[A-Z]\..{0,160}?\(((?:19|20)\d\d)[a-z]?\)\.?\s*[_*\"\u201c]*([^.?!\n]{12,220})", re.M)


def title_sim(a, b):
    import difflib
    na, nb = (re.sub(r"[^a-z0-9 ]", "", t.lower()) for t in (a, b))
    return difflib.SequenceMatcher(None, na, nb).ratio()


def cmd_fulltext(c, extra, workers):
    """Reference lists in the text write full citations ('Cockell, C. S. (2006). Title. Journal'). For a mention that is still open, take
    the title written there for the same first author and a year within one, look the title up, and tie the mention only when the
    record's title is a near copy (similarity 0.9 or more) and its first author and year agree."""
    titles = {}
    for pth, label in br.iter_texts(extra):
        t = br.read(pth)
        for m in FULLCITE.finditer(t):
            title = re.sub(r"[_*`]+", "", m.group(3)).strip(" ._")
            if len(title.split()) >= 3:
                titles.setdefault((br.fold(m.group(1)), int(m.group(2))), set()).add(title)
    todo = []
    for key in c["candidates"]["authoryear"]:
        v = status_of(c, key)
        if tied(v) or v.get("not_a_paper"):
            continue
        a, b, kind, y = split_key(key)
        for yy in (y, y - 1, y + 1):
            for t in titles.get((br.fold(a), yy), ()):
                todo.append((key, t, yy))
    br.log("open mentions with a full citation written somewhere:", len({k for k, _, _ in todo}), "titles:", len(todo))

    def one(item):
        key, title, yy = item
        a, b, kind, y = split_key(key)
        try:
            items = cr_items({"query.title": title, "query.author": a, "rows": "3", "mailto": br.MAILTO})
        except SearchFailed:
            return key, None
        for it in items:
            au = it.get("author") or []
            e = br.crossref_entry(it)
            if au and any(fam_match(x.get("family", ""), a) for x in au[:2]) and e["year"] and abs(int(e["year"]) - yy) <= 1 and title_sim(title, e["title"]) >= 0.9 \
                    and (kind != "pair" or any(fam_match(x.get("family", ""), b) for x in au[:4])):
                return key, (e, title)
        return key, None
    n = 0
    for key, res in br.parallel(one, todo, workers):
        v = status_of(c, key)
        if res and not tied(v):
            e, title = res
            e["via"] = "a full citation written in the text, title matched at Crossref"
            e["src"] = "crossref-search"
            c["entries"].setdefault(e["doi"], e)
            v["matched"] = [e["doi"]]
            v["evidence"] = "a full citation in the text gives the title (%s) and the record's first author, year and title agree" % title[:80]
            v.pop("review", None)
            n += 1
    br.log("tied by a full citation in the text:", n)


def cmd_hints(c, workers):
    """HINTS in references_curated.py give, for a key, the title (or the start of it) that a person identified from the citing sentences.
    The key is tied only when Crossref returns a record whose title is a near copy of the hint (similarity 0.75), whose author list
    contains the cited first author (and second author for a pair), and whose year is within two of the cited year."""
    items = [(k, t) for k, t in cur.HINTS.items() if k in c["candidates"]["authoryear"] and not tied(status_of(c, k))]

    def one(item):
        key, title = item
        a, b, kind, y = split_key(key)
        try:
            recs = cr_items({"query.bibliographic": title + " " + a + " " + str(y), "query.author": (a + " " + b).strip(), "rows": "8", "mailto": br.MAILTO})
        except SearchFailed:
            return key, None, "search failed"
        best = None
        for it in recs:
            au = it.get("author") or []
            e = br.crossref_entry(it)
            if not au or not e["year"] or abs(int(e["year"]) - y) > 2:
                continue
            if not any(fam_match(x.get("family", ""), a) for x in au[:3]):
                continue
            if kind == "pair" and not any(fam_match(x.get("family", ""), b) for x in au[:4]):
                continue
            if re.match(r"\s*(erratum|correction|corrigendum|comment on|reply)", e["title"], re.I):
                continue
            sim = max(title_sim(title, e["title"]), title_sim(title, e["title"][:len(title)]))
            if not best or sim > best[0]:
                best = (sim, e)
        if best and best[0] >= 0.75:
            return key, best[1], "sim %.2f" % best[0]
        return key, None, "no record close enough (%s)" % ("best %.2f: %s" % (best[0], best[1]["title"][:50]) if best else "none by this author")
    n = 0
    for key, e, why in br.parallel(one, items, workers):
        if e:
            v = status_of(c, key)
            e["via"] = "title identified from the citing sentences, matched at Crossref"
            e["src"] = "crossref-search"
            c["entries"].setdefault(e["doi"], e)
            v["matched"] = [e["doi"]]
            v["hand"] = True
            v["evidence"] = "the title was identified by hand from the citing sentences and the record agrees (%s)" % why
            v.pop("review", None)
            n += 1
        else:
            print("  hint failed:", key, "|", cur.HINTS[key][:50], "|", why)
    br.log("tied by a hand-identified title:", n, "of", len(items))


def cmd_local(c):
    ents = c["entries"]
    ids = entry_index(c)
    byfam = defaultdict(list)  # last token of a family name -> entry keys (first four authors)
    for k, e in ents.items():
        for au in (e.get("authors") or [])[:4]:
            byfam[br.fold(au[0])].append(k)
            toks = re.split(r"[\s\-]+", au[0])
            if len(toks) > 1:
                byfam[br.fold(toks[-1])].append(k)
    n_id = n_name = 0
    for key in list(c["candidates"]["authoryear"]):
        v = status_of(c, key)
        if settled(v):
            continue
        a, b, kind, y = split_key(key)
        hit = []
        # (a) an identifier beside the mention whose record fits
        for sn in snippets(c, key):
            for i in ids_in(sn):
                k = ids.get(i)
                if k and k not in hit and authors_fit(ents[k], a, b, kind, y):
                    hit.append(k)
        if hit:
            v["matched"] = hit
            v["evidence"] = "an identifier beside the mention (arXiv number or DOI) belongs to a paper whose authors and year fit"
            n_id += 1
            continue
        # (b) a record already in the list: author at any early position, year within one (pairs: both names, side by side)
        cand = [k for k in set(byfam.get(br.fold(a), [])) if authors_fit(ents[k], a, b, kind, y)]
        if kind == "one":
            cand = [k for k in cand if len(ents[k]["authors"]) == 1 or int(ents[k].get("year") or 0) == y]
        exact = [k for k in cand if int(ents[k].get("year") or 0) == y]
        pick = exact or cand
        if pick and (kind != "one" or len(pick) == 1 or exact):
            v["matched"] = sorted(pick)[:3]
            v["evidence"] = "a record already in the list has this author, year and (for a pair) second author"
            n_name += 1
    br.log("tied by an identifier beside the mention:", n_id, "tied to an existing record by name and year:", n_name)


# ---------------------------------------------------------------------------------------------------------------- search
def content_words(texts):
    w = set()
    for t in texts:
        w |= br.words(t)
    return w


def nums_after(texts):
    out = set()
    for t in texts:
        out |= set(re.findall(r"(?<![\d.])\d{1,6}(?![\d.]\d)", t))
    return out


DOMAIN_CT = ("Icarus Astrophysical Journal Astronomy Astrophysics Geophysical Research Planetary Space Science Nature Astrobiology Geochimica "
             "Monthly Notices Royal Astronomical Society Earth Letters Physical Review Chemistry Physics")
_VOCAB = {}


def vocab(c):
    """Words that occur in the titles and journals of at least three papers already in the list: the project's subject area, used to
    tell an astronomy, planetary or geoscience paper from a namesake in medicine or economics."""
    if not _VOCAB:
        from collections import Counter
        cnt = Counter()
        for e in c["entries"].values():
            cnt.update(br.words((e.get("title") or "") + " " + (e.get("journal") or "")))
        _VOCAB.update({w: n for w, n in cnt.items() if n >= 3})
    return _VOCAB


_JOURNALS = set()


def domain_journals(c):
    """Journals that at least two papers in the list appeared in: the project's own journals (Icarus, ApJ, GeoRL ...)."""
    if not _JOURNALS:
        from collections import Counter
        cnt = Counter((e.get("journal") or "").strip().lower() for e in c["entries"].values() if e.get("journal"))
        _JOURNALS.update(j for j, n in cnt.items() if n >= 2 and j)
    return _JOURNALS


def relevance(c, e):
    v = vocab(c)
    return len({w for w in br.words((e.get("title") or "") + " " + (e.get("journal") or "")) if w in v})


class SearchFailed(Exception):
    pass


def cr_items(params):
    """One Crossref query. A rate limit (HTTP 429) is waited out; a query that still fails raises SearchFailed so the key is retried
    on a later run instead of being recorded as 'nothing found'."""
    import time
    import urllib.error
    import urllib.request
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(params)
    for i in range(7):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": br.UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())["message"]["items"]
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return []
            time.sleep(min(60, 3 * 2 ** i))
        except (OSError, ValueError, KeyError):
            time.sleep(min(60, 3 * 2 ** i))
    raise SearchFailed(url[:200])


def score(c, item, a, b, kind, y, words, nums, text_fold):
    au = item.get("author") or []
    if not au:
        return None
    ok = False
    if kind == "pair":
        ok = any(fam_match(au[i].get("family", ""), a) and fam_match(au[i + 1].get("family", ""), b) for i in range(min(len(au) - 1, 4)))
    else:
        ok = fam_match(au[0].get("family", ""), a)
    if not ok:
        return None
    e = br.crossref_entry(item)
    if not e["year"] or abs(int(e["year"]) - y) > 1:
        return None
    if re.match(r"\s*(erratum|errata|correction|corrigendum|retraction|addendum|comment on|reply to|reply:|response to)", e["title"], re.I):
        return None  # a correction notice is never the paper that is cited
    title_hit = sorted(br.words(e["title"]) & words)
    vp = bool(e["volume"] and e["page"] and str(e["volume"]) in nums and str(e["page"]) in nums)
    jr = False
    cont = (item.get("container-title") or [""])[0]
    if cont:
        cw = [w for w in re.findall(r"[a-z]{4,}", cont.lower()) if w not in br.STOP]
        jr = bool(cw) and all(w in text_fold for w in cw[:2])
        if e["journal"] and re.search(r"\b" + re.escape(e["journal"]) + r"\b", text_fold, re.I):
            jr = True
    cites = item.get("is-referenced-by-count", 0) or 0
    rel = relevance(c, e)
    dom = (e.get("journal") or "").strip().lower() in domain_journals(c)
    s = (5 if vp else 0) + 2 * len(title_hit) + (2 if jr else 0) + (3 if int(e["year"]) == y else 0) + (1 if kind != "one" or len(au) == 1 else 0) \
        + min(rel, 6) + (4 if dom else 0) + math.log10(1 + cites) / 2
    return {"entry": e, "s": s, "dom": dom, "vp": vp, "title": title_hit, "journal": jr, "cites": cites, "exact": int(e["year"]) == y, "rel": rel}


def accept(ranked, kind):
    """The acceptance rule on the scored records (best first); returns (record or None, why, probable). A record is taken only when it
    has the cited first author (and second author, for a pair) and a year within one, AND the text supports it: the volume and first
    page, a title word, or the journal appear in the citing sentences; or it is the only such record in one of the project's journals
    and the citation is a pair. Anything weaker is 'probable' (kept for a person to check) and is not tied."""
    if not ranked:
        return None, "no record with this first author and year", None
    best = ranked[0]
    runner = ranked[1] if len(ranked) > 1 else None
    if best["vp"] and (not runner or best["s"] - runner["s"] >= 1):
        return best, "first author and year agree, and the volume and first page appear beside the mention", None
    dom = [r for r in ranked if r["dom"]]
    if not dom:
        return None, "no record in the project's journals (best: %s)" % (best["entry"]["journal"] or "no journal"), None
    top = dom[0]
    second = dom[1] if len(dom) > 1 else None
    if second and top["s"] - second["s"] < 2:
        return None, "two records in the project's journals fit equally (scores %.1f and %.1f)" % (top["s"], second["s"]), top
    if len(top["title"]) >= 2 and top["exact"]:
        return top, "first author and year agree, and title words appear in the citing sentences: " + ", ".join(top["title"][:5]), None
    if top["title"] and top["journal"]:
        return top, "first author and year agree, and a title word (%s) and the journal appear in the citing sentences" % top["title"][0], None
    if kind == "pair" and top["exact"] and (top["title"] or top["journal"]):
        return top, "both authors and the year agree, and %s appears in the citing sentences" % ("a title word (%s)" % top["title"][0] if top["title"] else "the journal"), None
    if kind == "pair" and top["exact"] and not second:
        return top, "the only record in the project's journals with both authors and this year", None
    return None, "one record in the project's journals fits the author and year, but nothing in the citing sentences confirms it", top


def search_key(c, key, override=None):
    a, b, kind, y = split_key(key)
    sn = snippets(c, key)
    if override:
        first, second, oy, tw = override
        a, b = first, second
        kind = "pair" if second else kind
        if oy:
            y = oy
        sn = sn + [tw]
    words = content_words(sn)
    nums = nums_after(sn)
    text_fold = " ".join(sn).lower()
    au_q = (a + " " + b).strip() if kind == "pair" else a
    flt = "from-pub-date:%d,until-pub-date:%d" % (y - 1, y + 1)
    trail = " ".join(re.findall(r"[A-Za-z0-9&.]+", (sn[0].split(str(y), 1)[1][:60] if sn and str(y) in sn[0] else "")))
    queries = [
        {"query.author": au_q, "query.container-title": DOMAIN_CT, "filter": flt, "rows": "60"},
        {"query.author": au_q, "query.bibliographic": " ".join(sorted(words))[:300], "filter": flt, "rows": "30"},
        {"query.bibliographic": ("%s %s %d %s" % (a, b, y, trail)).strip()[:300], "filter": flt, "rows": "10"},
    ]
    cands = {}
    for q in queries:
        for it in cr_items(dict(q, mailto=br.MAILTO)):
            d = (it.get("DOI") or "").lower()
            if d and d not in cands:
                sc = score(c, it, a, b, kind, y, words, nums, text_fold)
                if sc:
                    cands[d] = sc
    ranked = sorted(cands.values(), key=lambda r: -r["s"])
    got, why, probable = accept(ranked, kind)
    if got:
        return {"accepted": got, "why": why}
    near = [{"doi": r["entry"]["doi"], "title": r["entry"]["title"], "year": r["entry"]["year"], "journal": r["entry"]["journal"],
             "volume": r["entry"]["volume"], "page": r["entry"]["page"], "cites": r["cites"], "dom": r["dom"], "score": round(r["s"], 1)}
            for r in sorted(ranked, key=lambda r: (not r["dom"], -r["s"]))[:3]]
    return {"accepted": None, "why": why, "near": near, "probable": probable["entry"]["doi"] if probable else None}


def apply_result(c, key, res):
    v = status_of(c, key)
    tried = v.setdefault("tried", [])
    if "crossref" not in tried:
        tried.append("crossref")
    if res["accepted"]:
        e = res["accepted"]["entry"]
        e["via"] = "crossref search on author, year and the citing sentences"
        e["src"] = "crossref-search"
        c["entries"].setdefault(e["doi"], e)
        v["matched"] = [e["doi"]]
        v["evidence"] = res["why"]
        v.pop("review", None)
    else:
        v["review"] = {"why": res["why"], "near": res.get("near", []), "probable": res.get("probable")}


def cmd_search(c, limit, min_places, workers):
    aliases = {k: v for k, v in cur.ALIASES.items()}
    todo = []
    for key in c["candidates"]["authoryear"]:
        v = status_of(c, key)
        a, b, kind, y = split_key(key)
        if settled(v) or "crossref" in v.get("tried", []):
            continue
        if v["n"] < min_places or y < 1900:
            continue
        todo.append(key)
    todo.sort(key=lambda k: -status_of(c, k)["n"])
    if limit:
        todo = todo[:limit]
    br.log("keys to search:", len(todo))
    vocab(c), domain_journals(c)  # fill the shared word and journal sets before the threads start

    def one(key):
        try:
            return key, search_key(c, key, aliases.get(split_key(key)[0]))
        except SearchFailed as ex:  # a network hiccup must not lose the batch; the key stays untried
            return key, {"accepted": None, "why": "search failed: %s" % ex}
    done = found = 0
    for i in range(0, len(todo), 200):
        for key, res in br.parallel(one, todo[i:i + 200], workers):
            if res["why"].startswith("search failed"):
                continue
            apply_result(c, key, res)
            found += bool(res["accepted"])
            done += 1
        br.save_cache(c)
        br.log("searched", done, "tied", found)


# ---------------------------------------------------------------------------------------------------------------- report
def counts(c):
    out = {"tied": [0, 0], "not": [0, 0], "review": [0, 0]}
    for k, v in c.get("authoryear_status", {}).items():
        n = v.get("n", 0)
        g = "tied" if tied(v) else "not" if v.get("not_a_paper") else "review"
        out[g][0] += 1
        out[g][1] += n
    return out


def cmd_review(c, top, offset, min_places, probable_only):
    """Print the review list heaviest first for a person (or a careful reader) to adjudicate: the citing sentence and the nearest records."""
    rows = sorted(((v["n"], k, v) for k, v in c["authoryear_status"].items() if v.get("review") and k in c["candidates"]["authoryear"] and v["n"] >= min_places),
                  key=lambda r: (-r[0], r[1]))
    if probable_only:
        rows = [r for r in rows if r[2]["review"].get("probable")]
    for n, k, v in rows[offset:offset + top]:
        sn = snippets(c, k)
        a, b, kind, y = split_key(k)
        first = sn[0] if sn else ""
        i = first.find(str(y))
        print("%d  %s   [%s]" % (n, k, v["review"]["why"][:70]))
        print("   > " + first[max(0, i - 70):i + 130])
        for x in v["review"].get("near", [])[:3]:
            print("   - %s | %s | %s %s %s | c%s%s | %s" % (x["year"], x["title"][:70], x["journal"], x["volume"], x["page"], x["cites"], " dom" if x["dom"] else "", x["doi"]))


def probables(c, min_places):
    """The review keys that have one likely record (the only record in the project's journals with this author and year), heaviest first."""
    rows = [(v["n"], k, v) for k, v in c["authoryear_status"].items()
            if k in c["candidates"]["authoryear"] and v.get("review") and v["review"].get("probable") and not tied(v) and v["n"] >= min_places]
    return sorted(rows, key=lambda r: (-r[0], r[1]))


def cmd_probables(c, min_places, offset, top):
    for i, (n, k, v) in list(enumerate(probables(c, min_places)))[offset:offset + top]:
        sn = snippets(c, k)
        y = k.split("|")[3]
        t = sn[0] if sn else ""
        j = t.find(y)
        x = next((x for x in v["review"]["near"] if x["doi"] == v["review"]["probable"]), v["review"]["near"][0])
        print("%d %d %s :: %s :: %s (%s %s %s)" % (i, n, k.rstrip("|"), t[max(0, j - 55):j + 50].replace("\n", " "), x["title"][:70], x["journal"], x["volume"], x["page"]))


def cmd_accept(c, spec, min_places):
    """Tie the chosen probables (indices as `probables` prints them, '3,5,10-14') to their likely record, after a person checked each."""
    rows = probables(c, min_places)
    idx = set()
    for part in spec.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            idx |= set(range(int(lo), int(hi) + 1))
        elif part.strip():
            idx.add(int(part))
    done = 0
    for i in sorted(idx):
        n, k, v = rows[i]
        doi = v["review"]["probable"]
        res = None
        for it in cr_items({"filter": "doi:" + doi, "rows": "1", "mailto": br.MAILTO}):
            res = it
        if not res:
            continue
        e = br.crossref_entry(res)
        e["via"] = "chosen by hand from the citing sentences; author and year checked against the record"
        e["src"] = "crossref-search"
        c["entries"].setdefault(doi, e)
        v["matched"] = [doi]
        v["hand"] = True
        v["evidence"] = "the only record in the project's journals with this author and year, and a person read the citing sentence and agreed"
        v.pop("review", None)
        done += 1
    br.log("tied", done)


def cmd_lookup(c):
    """Read lines 'KEY<TAB>words from the title' on stdin and print the best Crossref records for each, for a person choosing a PICK."""
    for line in sys.stdin:
        if not line.strip():
            continue
        key, words_ = line.rstrip("\n").split("\t")[:2]
        a, b, kind, y = split_key(key)
        q = {"query.bibliographic": ("%s %s %d %s" % (a, b, y, words_)).strip(), "query.author": (a + " " + b).strip(), "rows": "5",
             "filter": "from-pub-date:%d,until-pub-date:%d" % (y - 1, y + 1), "mailto": br.MAILTO}
        try:
            items = cr_items(q)
        except SearchFailed:
            print(key, "SEARCH FAILED")
            continue
        print(key, "|", words_)
        for it in items[:4]:
            au = it.get("author") or []
            e = br.crossref_entry(it)
            ok = bool(au) and any(fam_match(x.get("family", ""), a) for x in au[:3])
            print("    %s %s | %s | %s %s %s | %s" % ("OK " if ok else "no ", e["year"], e["title"][:70], e["journal"], e["volume"], e["page"], e["doi"]))


def apply_curated(c):
    """Apply the hand decisions in references_curated.py: PICKS (key -> DOI or arXiv id, checked against the record), RECORDS (a
    pre-DOI paper or book written out in full, with a search link) and DROPS (key -> why it is not a paper after all)."""
    ents = c["entries"]
    for key, ref in cur.PICKS.items():
        if key not in c["candidates"]["authoryear"]:
            continue
        v = status_of(c, key)
        a, b, kind, y = split_key(key)
        if tied(v) and v.get("hand"):
            continue
        if ref.startswith("arXiv:"):
            got = br.fetch_arxiv([ref[6:]]).get(ref[6:])
            e = None if not got else {"authors": got["authors"], "year": got["year"], "title": got["title"], "arxiv": ref[6:], "doi": got["doi"],
                                     "journal": "", "volume": "", "page": "", "type": "preprint", "src": "arxiv", "journal_ref": got["journal_ref"]}
        else:
            e = br.fetch_crossref(ref) or br.fetch_datacite(ref)
        if not e:
            br.log("PICK not found:", key, ref)
            continue
        br.log("PICK", key, "->", (e.get("title") or "")[:60], e.get("year"))
        au = e.get("authors") or []
        fam_ok = (bool(au) and any(fam_match(x[0], a) for x in au[:3])) or key.split("|")[0] in cur.ORG_PICKS
        yr_ok = e.get("year") and abs(int(e["year"]) - y) <= 5
        if not (fam_ok and yr_ok):
            br.log("PICK refused (author or year disagree):", key, ref, au[:2], e.get("year"))
            continue
        k = ref if not ref.startswith("arXiv:") else ref
        e.setdefault("src", "crossref")
        e["via"] = "chosen by hand from the citing sentences; author and year checked against the record"
        ents.setdefault(k, e)
        v["matched"] = [k]
        v["hand"] = True
        v["evidence"] = "chosen by hand from the citing sentences; the record's author and year agree with the citation"
        v.pop("review", None)
    for key, rec in cur.RECORDS.items():
        if key not in c["candidates"]["authoryear"]:
            continue
        v = status_of(c, key)
        k = "manual:" + re.sub(r"[^a-z0-9]+", "-", key.lower()).strip("-")
        e = dict(rec)
        e.setdefault("type", "book" if rec.get("publisher") and not rec.get("journal") else "journal-article")
        e.setdefault("doi", "")
        e.setdefault("editors", [])
        e.setdefault("book", "")
        e["src"] = "hand"
        e["url"] = rec.get("url") or br.ads_title_url(rec["authors"][0][0], rec["year"], rec["title"])
        ents[k] = e
        v["matched"] = [k]
        v["hand"] = True
        v["evidence"] = "no DOI exists or was found; the record is written out from the citing sentences (journal, volume, page or publisher) with a search link"
        v.pop("review", None)
    for key, why in cur.DEMOTE.items():
        if key in c["candidates"]["authoryear"]:
            v = status_of(c, key)
            if tied(v) and not v.get("hand"):
                v["matched"] = []
                v.pop("evidence", None)
                v["review"] = {"why": why, "near": [], "probable": None}
                v.setdefault("tried", ["crossref"])
    for key, why in cur.DROPS.items():
        if key not in c["candidates"]["authoryear"]:
            continue
        v = status_of(c, key)
        if not tied(v):
            v["not_a_paper"] = {"what": why[0], "url": why[1] if len(why) > 1 else "", "class": "not a paper (decided by hand)"}
            v.pop("review", None)


def cmd_reset(c):
    """Forget every search result (not the triage or the local ties) so the search can be run again with a changed rule."""
    gone = [k for k, e in c["entries"].items() if e.get("via") == "crossref search on author, year and the citing sentences"]
    for k in gone:
        del c["entries"][k]
    for v in c.get("authoryear_status", {}).values():
        if "crossref" in v.get("tried", []):
            v.pop("tried", None)
            v.pop("review", None)
            if v.get("evidence", "").startswith(("first author", "both authors", "the only record")):
                v["matched"] = []
                v.pop("evidence", None)
    br.log("removed", len(gone), "searched entries")


def cmd_report(c):
    for g, (k, n) in counts(c).items():
        print("%-7s %5d keys %6d places" % (g, k, n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["triage", "local", "search", "report", "reset", "review", "curated", "lookup", "probables", "accept", "bare", "fulltext", "hints"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-places", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--probable", action="store_true")
    ap.add_argument("--idx", default="")
    ap.add_argument("--extra", nargs="*", default=[])
    a = ap.parse_args()
    c = br.load_cache()
    if a.cmd == "triage":
        cmd_triage(c)
    elif a.cmd == "local":
        cmd_local(c)
    elif a.cmd == "search":
        cmd_search(c, a.limit, a.min_places, a.workers)
    elif a.cmd == "reset":
        cmd_reset(c)
    elif a.cmd == "review":
        cmd_review(c, a.top, a.offset, a.min_places, a.probable)
        return
    elif a.cmd == "probables":
        cmd_probables(c, a.min_places, a.offset, a.top)
        return
    elif a.cmd == "hints":
        cmd_hints(c, a.workers)
    elif a.cmd == "fulltext":
        cmd_fulltext(c, a.extra, a.workers)
    elif a.cmd == "bare":
        cmd_bare(c)
    elif a.cmd == "accept":
        cmd_accept(c, a.idx, a.min_places)
    elif a.cmd == "lookup":
        cmd_lookup(c)
        return
    elif a.cmd == "curated":
        apply_curated(c)
    if a.cmd in ("triage", "local", "reset", "curated", "accept", "bare", "fulltext", "hints"):
        br.save_cache(c)
    cmd_report(c)


if __name__ == "__main__":
    main()
