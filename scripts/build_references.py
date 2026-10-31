#!/usr/bin/env python3
"""
Builds docs/REFERENCES.md and docs/REFERENCES.html: every paper, dataset and catalog the project cites, in the AAS/ADS style
that NASA journals use (author-year, alphabetical, journal, volume, page, DOI or arXiv or ADS link).

    python scripts/build_references.py harvest [--extra DIR ...]   # scan the repo (and extra folders) for citations -> candidates file
    python scripts/build_references.py resolve [--offline]         # look up each candidate (arXiv API, Crossref), update the cache
    python scripts/build_references.py render                      # write the .md and .html from the cache
    python scripts/build_references.py match                       # tie author-year mentions to papers by first author, year and second author
    python scripts/settle_references.py triage|local|search|report # label non-papers, tie the rest by identifier and Crossref (DOC.95)
    python scripts/build_references.py render --check              # fail if the committed pages differ from the cache

The cache `docs/reference/references-data.json` is the one record both pages are rendered from, so the two always match. A citation is
never guessed: an author-year mention that no identifier or search result can be tied to a paper with matching first author and year
is labelled "not a paper" (with what it is and a link) or listed "for review" with what was tried, and a test fails if any mention has none
of those three outcomes.
"""

import argparse
import html
import json
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
CACHE = os.path.join(ROOT, "docs", "reference", "references-data.json")
OUT_MD = os.path.join(ROOT, "docs", "REFERENCES.md")
OUT_HTML = os.path.join(ROOT, "docs", "REFERENCES.html")
MAILTO = "david.hagar@gmail.com"

TEXT_EXT = (".md", ".py", ".js", ".mjs", ".html", ".json", ".sh", ".txt", ".csv", ".tsv", ".bib", ".yml", ".yaml", ".rst", ".toml")
SKIP_PARTS = ("docs/REFERENCES", "docs/reference/references-data", "scripts/build_references.py", "scripts/settle_references.py", "scripts/references_curated.py", "/node_modules/", "/vendor/",
              "three.min", ".min.js", "package-lock")
MAX_BYTES = 4_000_000

RE_ARXIV_NEW = r"\d{4}\.\d{4,5}"
RE_ARXIV_OLD = r"[a-z]+(?:-[a-z]+)?(?:\.[A-Z]{2})?/\d{7}"
ARXIV = re.compile(r"(?:arXiv[:\s.]*(?:abs/)?|arxiv\.org/(?:abs|pdf)/)(" + RE_ARXIV_NEW + "|" + RE_ARXIV_OLD + r")(?:v\d+)?", re.I)
ARXIV_LIST = re.compile(r"arXiv[:\s]*((?:" + RE_ARXIV_NEW + r"(?:v\d+)?(?:\s*(?:,|and|;|&)\s*)?)+)", re.I)
DOI = re.compile(r"\b(10\.\d{4,9}/(?:[^\s\"'`<>,;|\\\[\]{}()]|\([^\s()]*\))+)")
BIBCODE = re.compile(r"\b((?:19|20)\d\d[A-Za-z][A-Za-z&.]{3}[.\d]{4}[A-Za-z.]\d{4}[A-Z])\b")
NAME = r"(?:[A-Z][a-z]?['’]?[A-Z]?[A-Za-zÀ-ſ-]{2,}|van [A-Z][a-z]+|de [A-Z][a-z]+)"
AUTHYEAR = re.compile(r"\b(" + NAME + r")(?:,? (?:and|&) (" + NAME + r")| et al\.?)?,? \(?((?:19|20)\d\d)([a-z])?\b")
NOT_NAMES = set("""Boss Research Status Quality Split Epoch Version Phase Part Table Figure Section Release Stage Item Items Round
Train Lane Report Reports Plan Study Census Sweep Draft Check Audit Test Tests Rule Rules Step Day Week Month January February March April
May June July August September October November December Mon Tue Wed Thu Fri Sat Sun Gen Core Aux Bug Star Tides Air Sources Planets Map
Db Sec Todo Main Merge Notes Note Data Phase Year Years Since Before After During Until From Between Around Roughly About Over Under Near
Late Early Circa Approx Total Count Counts Mark Release Released Filed Added Fixed Closed Done Open Order Version Versions Update Updated
Decided Decision Decisions Built Merged Opened Written Started Ended Dated Report Epochs Stage Era Eras Unit Units Mass Radius Age Time
Date Dates Number Numbers Size Sizes Level Levels Class Classes Type Types Model Models Method Methods Result Results Source Sources
Gate Gates Pass Passes Run Runs Row Rows Entry Entries List Lists Case Cases Set Sets Test Group Groups Task Tasks Page Pages Line Lines
Index Code Tier Tiers Tweak Tweaks Value Values Range Ranges Fill Fills Gap Gaps Win Wins Fix Fixes Bugs Wave Waves Batch Batches Slice""".split())

NOT_CITATIONS = {"EDT", "EST", "UTC", "Wikipedia", "OSHA", "Orbital", "TwinMind", "NASA", "Python", "GitHub", "Claude", "Anthropic", "Alembic", "Redis", "MySQL"}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def iter_texts(extra):
    """Yield (label, text) for every text file of the repo and of the extra folders."""
    files = subprocess.check_output(["git", "ls-files"], cwd=ROOT, text=True).split("\n")
    for f in files:
        if f and f.endswith(TEXT_EXT) and not any(s in "/" + f for s in SKIP_PARTS):
            yield os.path.join(ROOT, f), f
    for d in extra:
        for dp, dns, fns in os.walk(d):
            for fn in fns:
                if fn.endswith(TEXT_EXT):
                    p = os.path.join(dp, fn)
                    yield p, os.path.relpath(p, d)


def read(p):
    try:
        if os.path.getsize(p) > MAX_BYTES:
            return ""
        with open(p, errors="ignore") as fh:
            return fh.read()
    except OSError:
        return ""


# Publisher stems that a wrapped line leaves behind, and DOIs mistyped in the text with the right one (checked against Crossref).
DOI_STEMS = {"10.1051/0004-6361", "10.1051/0004-"}
DOI_FIXES = {
    "10.1007/978-1-4020-6933-8-7": "10.1007/978-1-4020-6933-8_7",
    "10.1146/annurev-astro081915-023441": "10.1146/annurev-astro-081915-023441",
    "10.1111/cgf13202": "10.1111/cgf.13202",
    "10.1016/0019-1035(63)90074-6": "10.1016/0019-1035(63)90072-1",
    "10.1038/s41550-020-1069-3": "10.1038/s41550-020-1069-4",
    "10.3847/1538-4357/abf22d": "10.3847/1538-4357/abfd9c",
    "10.1007/s11214-015-0158-y": "10.1146/annurev-earth-060313-054834",
}


def clean_doi(d):
    d = d.rstrip(".:")
    while d and d[-1] in ")]" and d.count("(") < d.count(")"):
        d = d[:-1]
    if d.count("(") > d.count(")") or re.search(r"/(?:\d+)?$", d):
        return None
    d = re.sub(r"(?i)\.(?:pdf|html?|xml)$", "", d)
    d = re.sub(r"(?i)/(?:pdf|full|epdf|abstract)$", "", d).replace("\u2013", "-").lower()
    # a DOI cut at a line break in the text ends in a hyphen or is only a publisher stem; it names no paper
    if d.endswith("-") or d in DOI_STEMS:
        return None
    d = DOI_FIXES.get(d, d)
    return d if len(d) > 9 else None


def harvest(extra):
    arx, dois, bibs, ays = Counter(), Counter(), Counter(), Counter()
    ay_ctx = {}
    ay_more = defaultdict(list)
    for p, label in iter_texts(extra):
        t = read(p)
        if not t:
            continue
        for m in ARXIV.finditer(t):
            arx[re.sub(r"v\d+$", "", m.group(1).lower() if "/" in m.group(1) else m.group(1))] += 1
        for m in ARXIV_LIST.finditer(t):
            for i in re.findall(RE_ARXIV_NEW, m.group(1)):
                arx[i] += 1
        for m in DOI.finditer(t):
            d = clean_doi(m.group(1))
            if d:
                dois[d] += 1
        for m in BIBCODE.finditer(t):
            bibs[m.group(1)] += 1
        for m in AUTHYEAR.finditer(t):
            a, b, y, s = m.group(1), m.group(2), m.group(3), m.group(4)
            if a in NOT_NAMES or (b and b in NOT_NAMES):
                continue
            if not (1900 <= int(y) <= 2026):
                continue
            kind = "pair" if b else ("etal" if "et al" in m.group(0) else "one")
            key = (a, b or "", kind, y)
            ays[key] += 1
            if key not in ay_ctx:
                ay_ctx[key] = t[max(0, m.start() - 90):m.end() + 90].replace("\n", " ")
            # up to six different citing places per key (60 characters before, 220 after), kept so the settle step can read a
            # journal, volume and page or title words that follow the author and year
            if len(ay_more[key]) < 6:
                snip = t[max(0, m.start() - 60):m.end() + 220].replace("\n", " ")
                if all(snip[:80] != o[:80] for o in ay_more[key]):
                    ay_more[key].append(snip)
    return arx, dois, bibs, ays, ay_ctx, ay_more


def load_cache():
    if os.path.exists(CACHE):
        with open(CACHE) as fh:
            return json.load(fh)
    return {"entries": {}, "candidates": {}, "authoryear": {}, "unresolved_ids": {}}


def save_cache(c):
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    with open(CACHE, "w") as fh:
        json.dump(c, fh, indent=1, sort_keys=True, ensure_ascii=False)
        fh.write("\n")


def cmd_harvest(args):
    arx, dois, bibs, ays, ctx, more = harvest(args.extra)
    c = load_cache()
    c["candidates"] = {
        "arxiv": dict(arx), "doi": dict(dois), "bibcode": dict(bibs),
        "authoryear": {"|".join(k): {"n": n, "ctx": ctx[k], "ctxs": more[k]} for k, n in ays.items()},
    }
    counts = find_archives(args.extra)
    c["archives"] = [{"name": ARCHIVES[i][0], "provider": ARCHIVES[i][1], "url": ARCHIVES[i][2], "mentions": n} for i, n in sorted(counts.items())]
    c["setup_datasets"] = catalogue_datasets()
    save_cache(c)
    log("archives", len(c["archives"]), "setup datasets", len(c["setup_datasets"]))
    log("arXiv", len(arx), "DOI", len(dois), "bibcode", len(bibs), "author-year", len(ays))


UA = "planetGen-references/1.0 (mailto:%s)" % MAILTO


def http_get(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json, application/atom+xml, */*"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (404, 400, 410):
                return None
            last = e
        except Exception as e:  # network trouble: wait and retry
            last = e
        time.sleep(2 * (i + 1))
    log("giving up on", url, last)
    return None


ATOM = "{http://www.w3.org/2005/Atom}"
ARX = "{http://arxiv.org/schemas/atom}"


def split_name(full):
    full = full.strip()
    if "," in full:
        fam, giv = full.split(",", 1)
        return fam.strip(), giv.strip()
    parts = full.split()
    if len(parts) == 1:
        return parts[0], ""
    i = len(parts) - 1
    while i > 0 and parts[i - 1].lower() in ("van", "von", "de", "der", "den", "del", "di", "da", "la", "le", "ten", "ter", "du", "dos", "das"):
        i -= 1
    return " ".join(parts[i:]), " ".join(parts[:i])


def initials(giv):
    out = []
    for w in re.split(r"[\s]+", giv.strip()):
        if not w:
            continue
        for part in w.split("-"):
            c = part.strip(".")[:1]
            if c:
                out.append(c.upper() + ".")
    return " ".join(out)


def fetch_arxiv(ids):
    out = {}
    for i in range(0, len(ids), 60):
        chunk = ids[i:i + 60]
        data = http_get("http://export.arxiv.org/api/query?max_results=%d&id_list=%s" % (len(chunk), ",".join(chunk)))
        time.sleep(3)
        if not data:
            continue
        try:
            root = ET.fromstring(data)
        except ET.ParseError:
            continue
        for e in root.findall(ATOM + "entry"):
            eid = e.findtext(ATOM + "id", "")
            m = re.search(r"abs/(.+?)(?:v\d+)?$", eid)
            title = e.findtext(ATOM + "title", "")
            if not m or title.strip().lower() == "error":
                continue
            aid = m.group(1)
            authors = [split_name(a.findtext(ATOM + "name", "")) for a in e.findall(ATOM + "author")]
            out[aid] = {
                "authors": [[f, initials(g)] for f, g in authors],
                "title": re.sub(r"\s+", " ", title).strip(),
                "year": int((e.findtext(ATOM + "published", "") or "0")[:4]) or None,
                "doi": (e.findtext(ARX + "doi") or "").lower() or None,
                "journal_ref": e.findtext(ARX + "journal_ref"),
            }
    return out


JOURNALS = {
    "the astrophysical journal": "ApJ", "the astrophysical journal letters": "ApJL", "the astrophysical journal supplement series": "ApJS",
    "the astronomical journal": "AJ", "astronomy & astrophysics": "A&A", "astronomy and astrophysics": "A&A",
    "monthly notices of the royal astronomical society": "MNRAS", "monthly notices of the royal astronomical society: letters": "MNRAS",
    "icarus": "Icarus", "nature": "Nature", "nature astronomy": "NatAs", "nature geoscience": "NatGe", "nature communications": "NatCo",
    "science": "Science", "science advances": "SciA", "planetary and space science": "P&SS", "geophysical research letters": "GeoRL",
    "journal of geophysical research: planets": "JGRE", "journal of geophysical research: space physics": "JGRA",
    "journal of geophysical research": "JGR", "journal of geophysical research: solid earth": "JGRB",
    "journal of geophysical research: atmospheres": "JGRD", "geochimica et cosmochimica acta": "GeCoA",
    "earth and planetary science letters": "E&PSL", "annual review of astronomy and astrophysics": "ARA&A",
    "annual review of earth and planetary sciences": "AREPS", "publications of the astronomical society of the pacific": "PASP",
    "publications of the astronomical society of japan": "PASJ", "publications of the astronomical society of australia": "PASA",
    "space science reviews": "SSRv", "reviews of modern physics": "RvMP", "physical review letters": "PhRvL", "physical review d": "PhRvD",
    "physical review a": "PhRvA", "physical review e": "PhRvE", "physical review c": "PhRvC", "physical review b": "PhRvB",
    "physical review": "PhRv", "astrobiology": "AsBio", "international journal of astrobiology": "IJAsB",
    "journal of quantitative spectroscopy and radiative transfer": "JQSRT", "the journal of chemical physics": "JChPh",
    "journal of physical and chemical reference data": "JPCRD", "atomic data and nuclear data tables": "ADNDT",
    "astronomy and computing": "A&C", "astronomische nachrichten": "AN", "celestial mechanics and dynamical astronomy": "CeMDA",
    "reports on progress in physics": "RPPh", "reviews in mineralogy and geochemistry": "RvMG", "chemical reviews": "ChRv",
    "the planetary science journal": "PSJ", "research notes of the aas": "RNAAS", "the open journal of astrophysics": "OJAp",
    "astrophysics and space science": "Ap&SS", "solar physics": "SoPh", "living reviews in solar physics": "LRSP",
    "geophysical journal international": "GeoJI", "physics of the earth and planetary interiors": "PEPI",
    "meteoritics & planetary science": "M&PS", "meteoritics and planetary science": "M&PS", "pure and applied geophysics": "PApGe",
    "journal of atmospheric sciences": "JAtS", "journal of the atmospheric sciences": "JAtS", "applied optics": "ApOpt",
    "journal of geophysical research: oceans": "JGRC", "geochemistry, geophysics, geosystems": "GGG", "science of the total environment": "STOTE",
    "radiation research": "RadR", "annals of the icrp": "AnICRP", "journal of the american chemical society": "JAChS",
    "the journal of physical chemistry a": "JPCA", "chemical physics letters": "CPL", "chemical physics": "CP",
    "the journal of physical chemistry": "JPhCh", "journal of molecular spectroscopy": "JMoSp", "the astronomical journal supplement": "AJS",
    "advances in space research": "AdSpR", "acta astronautica": "AcAau", "earth, moon, and planets": "EM&P", "experimental astronomy": "ExA",
    "proceedings of the national academy of sciences": "PNAS", "nuclear physics a": "NuPhA", "astroparticle physics": "APh",
    "new astronomy reviews": "NewAR", "new astronomy": "NewA", "the astronomy and astrophysics review": "A&ARv",
    "the observatory": "Obs", "physics reports": "PhR", "physics letters b": "PhLB", "journal of fluid mechanics": "JFM",
    "journal of the royal astronomical society of canada": "JRASC", "bulletin of the american astronomical society": "BAAS",
    "classical and quantum gravity": "CQGra", "journal of cosmology and astroparticle physics": "JCAP", "universe": "Univ",
    "galaxies": "Galax", "remote sensing": "RemS", "atmosphere": "Atmos", "life": "Life", "geosciences": "Geosci",
}


def abbrev_journal(long, short):
    k = re.sub(r"\s+", " ", (long or "").strip().lower().replace("&amp;", "&"))
    if k in JOURNALS:
        return JOURNALS[k]
    return (short or long or "").strip()


def crossref_entry(d):
    msg = d
    authors = [[a.get("family") or a.get("name") or "", initials(a.get("given", ""))] for a in msg.get("author", []) if a.get("family") or a.get("name")]
    issued = msg.get("issued", {}).get("date-parts", [[None]])[0][0] or msg.get("published", {}).get("date-parts", [[None]])[0][0] \
        or msg.get("published-online", {}).get("date-parts", [[None]])[0][0]
    ct = (msg.get("container-title") or [""])[0]
    sct = (msg.get("short-container-title") or [""])[0]
    page = msg.get("page") or msg.get("article-number") or ""
    first_page = re.split(r"[-\u2013]", page)[0] if page else ""
    typ = msg.get("type", "")
    return {
        "authors": authors, "year": issued, "title": re.sub(r"<[^>]+>", "", re.sub(r"\s+", " ", (msg.get("title") or [""])[0])).strip(),
        "journal": abbrev_journal(ct, sct), "volume": msg.get("volume", ""), "page": first_page, "type": typ,
        "publisher": msg.get("publisher", ""), "doi": msg.get("DOI", "").lower(),
        "editors": [[e.get("family", ""), initials(e.get("given", ""))] for e in msg.get("editor", [])] if typ.startswith("book") or typ == "dataset" else [],
        "book": ct if typ in ("book-chapter", "proceedings-article") else "",
    }


def fetch_crossref(doi):
    data = http_get("https://api.crossref.org/works/%s?mailto=%s" % (urllib.parse.quote(doi, safe="/()"), MAILTO))
    if not data:
        return None
    try:
        return crossref_entry(json.loads(data)["message"])
    except (ValueError, KeyError):
        return None


def fetch_datacite(doi):
    data = http_get("https://api.datacite.org/dois/%s" % urllib.parse.quote(doi, safe="/()"))
    if not data:
        return None
    try:
        a = json.loads(data)["data"]["attributes"]
    except (ValueError, KeyError):
        return None
    au = []
    for c in a.get("creators", []):
        fam = c.get("familyName") or c.get("name", "")
        au.append([fam, initials(c.get("givenName", ""))])
    return {
        "authors": au, "year": a.get("publicationYear"), "title": (a.get("titles") or [{}])[0].get("title", ""), "journal": "", "volume": "",
        "page": "", "type": "dataset", "publisher": a.get("publisher", ""), "doi": doi, "editors": [], "book": "",
    }


def parallel(fn, items, workers=8):
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(fn, items))


def cmd_resolve(args):
    c = load_cache()
    cand = c["candidates"]
    ents = c["entries"]
    have_doi = {e["doi"]: k for k, e in ents.items() if e.get("doi")}
    have_arx = {e["arxiv"]: k for k, e in ents.items() if e.get("arxiv")}
    bad = c.setdefault("unresolved_ids", {})
    for k in [k for k in bad if k not in cand["doi"] and k not in cand["arxiv"]]:
        del bad[k]  # an identifier the text no longer carries (or that the harvest now repairs) is not a standing failure
    todo_arx = [a for a in cand["arxiv"] if a not in have_arx and a not in bad]
    log("arXiv to look up:", len(todo_arx))
    got = fetch_arxiv(todo_arx)
    for a in todo_arx:
        g = got.get(a)
        if not g:
            bad[a] = "arXiv id not found by the arXiv API"
            continue
        e = {"authors": g["authors"], "year": g["year"], "title": g["title"], "arxiv": a, "doi": g["doi"], "journal": "", "volume": "", "page": "",
             "type": "preprint", "src": "arxiv", "journal_ref": g["journal_ref"]}
        ents["arXiv:" + a] = e
        have_arx[a] = "arXiv:" + a
    dois = set(cand["doi"]) | {e["doi"] for e in ents.values() if e.get("doi")}
    todo = [d for d in dois if d not in bad and not any(k == d for k in ents) and not (d in have_doi and ents[have_doi[d]].get("src") == "crossref")]
    log("DOIs to look up:", len(todo))

    def one(d):
        r = fetch_crossref(d) or fetch_datacite(d)
        return d, r
    for d, r in parallel(one, todo):
        if not r:
            bad[d] = "DOI not found at Crossref or DataCite"
            continue
        r["src"] = "crossref"
        # merge with an arXiv entry that has this DOI
        k = have_doi.get(d)
        if k and k in ents:
            r["arxiv"] = ents[k].get("arxiv")
            del ents[k]
        ents[d] = r
        have_doi[d] = d
        if r.get("arxiv"):
            have_arx[r["arxiv"]] = d
    # link arXiv ids named by the arXiv entries that now share a DOI
    for k in list(ents):
        e = ents[k]
        if k.startswith("arXiv:") and e.get("doi") and e["doi"] in ents:
            ents[e["doi"]]["arxiv"] = e["arxiv"]
            del ents[k]
    save_cache(c)
    log("entries:", len(ents), "unresolved ids:", len(bad))


def fold(s):
    import unicodedata
    return re.sub(r"[^a-z]", "", unicodedata.normalize("NFKD", s.replace("\u00df", "ss")).encode("ascii", "ignore").decode().lower())


STOP = set("the of and for with from that this are was were its their into over under between about using use based new first one two".split())


def words(s):
    return {w for w in re.findall(r"[a-z]{4,}", s.lower()) if w not in STOP}


def match_authoryear(c):
    """Tie each author-year mention to a resolved entry by first author, year and (for pairs) second author. Never a guess:
    a mention that fits more than one distinct paper is kept apart as ambiguous."""
    ents = c["entries"]
    idx = defaultdict(list)
    for k, e in ents.items():
        if e.get("authors") and e.get("year"):
            idx[(fold(e["authors"][0][0]), int(e["year"]))].append(k)
    out = {}
    for key, v in c["candidates"]["authoryear"].items():
        a, b, kind, y = key.split("|")
        y = int(y)
        hits = []
        for yy in (y, y - 1, y + 1):
            for k in idx.get((fold(a), yy), []):
                e = ents[k]
                if kind == "pair" and not (len(e["authors"]) > 1 and fold(e["authors"][1][0]) == fold(b)):
                    continue
                if kind == "one" and len(e["authors"]) > 1 and yy != y:
                    continue
                if (k, yy == y) not in hits:
                    hits.append((k, yy == y))
        exact = [k for k, ex in hits if ex]
        pick = exact or [k for k, _ in hits]
        out[key] = {"n": v["n"], "matched": pick}
    return out


# Archives, catalogs, tables and software the project reads that are not a single paper. Each is listed only when the harvest finds
# its pattern in the repository (or the extra folders), so the list never names something the project does not use.
ARCHIVES = (
    ("NIST Chemistry WebBook, NIST Standard Reference Database 69", "National Institute of Standards and Technology", "https://webbook.nist.gov/chemistry/", r"NIST (?:Chemistry )?WebBook|webbook\.nist"),
    ("NIST-JANAF Thermochemical Tables", "National Institute of Standards and Technology", "https://janaf.nist.gov/", r"JANAF"),
    ("Particle Data Group, atomic and nuclear properties", "Particle Data Group, Lawrence Berkeley National Laboratory", "https://pdg.lbl.gov/", r"Particle Data Group|\bPDG\b"),
    ("USGS Spectral Library Version 7 (splib07a)", "U.S. Geological Survey", "https://doi.org/10.3133/ds1035", r"splib07|USGS (?:spectral|splib)"),
    ("PubChem", "National Center for Biotechnology Information", "https://pubchem.ncbi.nlm.nih.gov/", r"PubChem"),
    ("Crystallography Open Database", "Crystallography Open Database project", "https://www.crystallography.net/cod/", r"Crystallography Open Database"),
    ("Webmineral Mineralogy Database", "David Barthelmy", "http://webmineral.com/", r"[Ww]ebmineral"),
    ("Handbook of Mineralogy", "Mineral Data Publishing and the Mineralogical Society of America", "https://www.handbookofmineralogy.org/", r"Handbook of Mineralogy"),
    ("BurnMan, a mineral physics toolkit", "Cottaar et al.; CIG Geodynamics", "https://geodynamics.github.io/burnman/", r"[Bb]urn[Mm]an"),
    ("ExoPlex, a Python toolkit for exoplanet interiors", "Unterborn et al.", "https://github.com/CaymanUnterborn/ExoPlex", r"ExoPlex"),
    ("SeaFreeze, thermodynamics of water and ice", "Journaux et al.", "https://github.com/Bjournaux/SeaFreeze", r"SeaFreeze"),
    ("MESA Isochrones and Stellar Tracks (MIST) v1.2", "Harvard-Smithsonian Center for Astrophysics", "https://waps.cfa.harvard.edu/MIST/", r"\bMIST\b"),
    ("ExoMol database and ExoMolOP opacities", "ExoMol project, University College London", "https://exomol.com/", r"ExoMol"),
    ("HITRAN molecular spectroscopic database", "Harvard-Smithsonian Center for Astrophysics", "https://hitran.org/", r"HITRAN"),
    ("Harris catalog of Milky Way globular cluster parameters (2010 edition)", "McMaster University", "https://physwww.mcmaster.ca/~harris/Databases.html", r"Harris (?:1996|catalog)"),
    ("Baraffe et al. (2015) BHAC15 evolutionary models", "Centre de Recherche Astrophysique de Lyon", "https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/", r"BHAC15"),
    ("Castelli and Kurucz ATLAS9 model-atmosphere grids", "Osservatorio Astronomico di Trieste", "https://wwwuser.oats.inaf.it/castelli/grids.html", r"ATLAS9|Castelli and Kurucz"),
    ("PHOENIX-ACES stellar spectra library", "Goettingen Spectral Library", "https://phoenix.astro.physik.uni-goettingen.de/", r"PHOENIX-ACES|PHOENIX"),
    ("MPI-Mainz UV/VIS Spectral Atlas of Gaseous Molecules", "Max Planck Institute for Chemistry, Mainz", "https://doi.org/10.5281/zenodo.6951", r"MPI-Mainz"),
    ("NOAA Space Weather Prediction Center solar proton events", "NOAA SWPC", "https://www.swpc.noaa.gov/products/solar-proton-events", r"NOAA SWPC|Solar Proton Events"),
    ("NASA Exoplanet Archive", "NASA Exoplanet Science Institute, Caltech/IPAC", "https://exoplanetarchive.ipac.caltech.edu/", r"Exoplanet Archive"),
    ("VizieR catalogue access tool", "Centre de Donnees astronomiques de Strasbourg", "https://vizier.cds.unistra.fr/", r"VizieR"),
    ("SIMBAD astronomical database", "Centre de Donnees astronomiques de Strasbourg", "https://simbad.cds.unistra.fr/", r"SIMBAD"),
    ("Gaia archive and data releases", "European Space Agency, Gaia Data Processing and Analysis Consortium", "https://gea.esac.esa.int/archive/", r"Gaia (?:DR|archive|Data Release)|Gaia Collaboration"),
    ("JPL Small-Body Database", "NASA Jet Propulsion Laboratory", "https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html", r"Small-Body Database|JPL SBDB|\bSBDB\b"),
    ("JPL Solar System Dynamics: planetary satellite physical and mean orbital parameters", "NASA Jet Propulsion Laboratory", "https://ssd.jpl.nasa.gov/sats/", r"ssd\.jpl\.nasa\.gov"),
    ("NASA Planetary Data System", "NASA", "https://pds.nasa.gov/", r"Planetary Data System|\bPDS\b"),
    ("NMSU planetary datasets (Kuiper-belt, giant-planet infrared spectra)", "New Mexico State University", "https://astronomy.nmsu.edu/", r"NMSU"),
    ("Local Volume Database", "Pace et al.", "https://github.com/apace7/local_volume_database", r"Local Volume Database"),
    ("Mean dwarf stellar color and effective temperature sequence (E. Mamajek, 2022.04.16)", "University of Rochester", "https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt", r"Mamajek"),
    ("UNSCEAR 2008 Report, Annex B", "United Nations Scientific Committee on the Effects of Atomic Radiation", "https://www.unscear.org/unscear/en/publications/2008_1.html", r"UNSCEAR"),
    ("ICRP Publication 103", "International Commission on Radiological Protection", "https://www.icrp.org/publication.asp?id=ICRP+Publication+103", r"ICRP"),
    ("Cantera, chemical kinetics and thermodynamics (NASA-7 polynomials)", "Cantera developers", "https://cantera.org/", r"[Cc]antera"),
    ("Edinburgh DataShare", "University of Edinburgh", "https://datashare.ed.ac.uk/", r"Edinburgh DataShare"),
    ("Zenodo open repository", "CERN and OpenAIRE", "https://zenodo.org/", r"[Zz]enodo"),
    ("Three.js JavaScript 3D library", "mrdoob and contributors", "https://threejs.org/", r"three\.js|[Tt]hree\.js"),
)


def find_archives(extra):
    counts = Counter()
    pats = [(i, re.compile(a[3])) for i, a in enumerate(ARCHIVES)]
    for p, label in iter_texts(extra):
        t = read(p)
        if not t:
            continue
        for i, rx in pats:
            if rx.search(t):
                counts[i] += 1
    return counts


def catalogue_datasets():
    sys.path.insert(0, os.path.join(ROOT, "src"))
    try:
        from planetgen.datasets import catalogue
    except Exception as exc:  # the package may need its dependencies
        log("catalogue import failed:", exc)
        return []
    out = []
    for d in catalogue.ALL:
        host = urllib.parse.urlsplit(d.pieces[0].url)
        out.append({"key": d.key, "title": d.title, "source": d.source, "url": "%s://%s/" % (host.scheme, host.netloc),
                    "file_url": d.pieces[0].url})
    return out


# ---------------------------------------------------------------- citation IDs (DOC.100)

ID_PREFIX = "R"
ID_WIDTH = 4


def format_id(number):
    """`R0123`: the ID of the 123rd source ever listed. Four digits at least; the number grows past 9999 without changing earlier IDs."""
    return "%s%0*d" % (ID_PREFIX, ID_WIDTH, number)


def identifiers(key, e):
    """Every name one source goes by in the cache, lower case: its key (a DOI, `arXiv:<id>` or a hand-made key), its DOI and its arXiv id.
    A source whose arXiv preprint later gets a DOI keeps the ID it had, because the arXiv name is still among these."""
    names = {str(key).lower()}
    if e.get("doi"):
        names.add(e["doi"].lower())
    if e.get("arxiv"):
        names.add("arxiv:" + e["arxiv"].lower())
    return sorted(names)


def assign_ids(c):
    """
    Gives every source in the list (papers, books, datasets, and the archives and setup datasets) one ID, once, and keeps it.

    The cache holds `ids` (each name a source has ever gone by -> its ID), `id_next` (the number the next new source takes; an ID is never
    reused, even when its source is removed) and `id_aliases` (an ID that was merged into another -> the one that stays, so an old ID in a text
    still resolves). A source is looked up by all of its names: the first ID found is its ID, and a second ID found means two records have
    become one paper, so the higher ID becomes an alias of the lower. A source with no ID takes the next number, in list order.
    Returns the number of IDs assigned.
    """
    ids = c.setdefault("ids", {})
    aliases = c.setdefault("id_aliases", {})
    nxt = c.get("id_next") or (1 + max([int(v[len(ID_PREFIX):]) for v in ids.values()] + [0]))
    given = 0
    sources = []
    for k, e in c["entries"].items():
        if e.get("authors") or e.get("title"):
            sources.append((k, identifiers(k, e)))
    for a in c.get("archives", []):
        sources.append(("archive:" + a["url"], ["archive:" + a["url"].lower()]))
    for d in c.get("setup_datasets", []):
        sources.append(("setup:" + d["url"], ["setup:" + d["url"].lower()]))
    # list order, so a first build numbers the sources in the order the page prints them
    ents = c["entries"]
    order = sorted((s for s in sources if s[0] in ents), key=lambda s: (auth_key(ents[s[0]]), int(ents[s[0]].get("year") or 0), (ents[s[0]].get("title") or "").lower(), s[0]))
    order += [s for s in sources if s[0] not in ents]
    for k, names in order:
        found = sorted({ids[n] for n in names if n in ids})
        if found:
            keep = found[0]
            for other in found[1:]:
                aliases[other] = keep
        else:
            keep = format_id(nxt)
            nxt += 1
            given += 1
        for n in names:
            ids[n] = keep
    c["id_next"] = nxt
    return given


def id_of(c, key):
    """The ID of the source whose key is `key` (a key of `entries`), after `assign_ids`."""
    return c["ids"][str(key).lower()]


def resolve_id(c, text):
    """The ID a text names, following an alias; `None` when it names no source."""
    text = c.get("id_aliases", {}).get(text, text)
    return text if text in c.get("ids", {}).values() else None


# ---------------------------------------------------------------- rendering

def auth_key(e):
    return [(fold(f), i) for f, i in e["authors"]] or [("", "")]


def fmt_authors(authors):
    """AAS/ADS author list: 'Family, I. I., Family, I., & Family, I.'; more than five authors shorten to 'First, I., et al.'"""
    names = ["%s, %s" % (f, i) if i else f for f, i in authors]
    if not names:
        return ""
    if len(names) == 1:
        return names[0]
    if len(names) > 5:
        return names[0] + ", et al."
    if len(names) == 2:
        return names[0] + " & " + names[1]
    return ", ".join(names[:-1]) + ", & " + names[-1]


def ads_title_url(first, year, title):
    """An ADS search link for a paper that has no DOI: first author, year and a few title words."""
    words_ = " ".join(re.findall(r"[A-Za-z]{4,}", title)[:6])
    return "https://ui.adsabs.harvard.edu/search/q=" + urllib.parse.quote('author:"^%s" year:%s title:(%s)' % (first, year, words_))


def entry_urls(e):
    """The primary URL (DOI when there is one, else arXiv, else the search link a hand-written record carries) and the ADS page for it."""
    doi, ax = e.get("doi"), e.get("arxiv")
    if not doi and not ax and e.get("url"):
        return e["url"], e["url"]
    if doi:
        return "https://doi.org/" + urllib.parse.quote(doi, safe="/()"), "https://ui.adsabs.harvard.edu/abs/" + urllib.parse.quote(doi, safe="/()") + "/abstract"
    if ax:
        return "https://arxiv.org/abs/" + ax, "https://ui.adsabs.harvard.edu/abs/arXiv:" + ax + "/abstract"
    return None, None


def entry_parts(e, suffix, sid=""):
    """The pieces of one reference, shared by the markdown and the HTML writers so the two pages carry the same text. `sid` is its ID."""
    authors = fmt_authors(e["authors"]) or (e.get("publisher") or "Anonymous")
    year = "%s%s" % (e.get("year") or "n.d.", suffix)
    journal = e.get("journal") or ""
    bits = []
    if e.get("book"):
        bits.append("in " + e["book"])
        if e.get("editors"):
            ed = fmt_authors(e["editors"]).replace(", et al.", " et al.")
            bits.append("ed. %s" % ed)
        if e.get("publisher"):
            bits.append(e["publisher"])
    elif journal:
        bits.append(journal)
        if e.get("volume"):
            bits.append(str(e["volume"]))
        if e.get("page"):
            bits.append(str(e["page"]))
    elif e.get("type") == "dataset" or e.get("publisher"):
        if e.get("publisher"):
            bits.append(e["publisher"])
    elif e.get("journal_ref"):
        bits.append("published as " + e["journal_ref"])
    ids = []
    if e.get("doi"):
        ids.append("doi:" + e["doi"])
    if e.get("arxiv"):
        ids.append("arXiv:" + e["arxiv"])
    url, ads = entry_urls(e)
    return {"authors": authors, "year": year, "title": e.get("title") or "", "bits": bits, "ids": ids, "url": url, "ads": ads,
            "kind": e.get("type") or "", "sid": sid}


def sorted_entries(ents):
    keyed = sorted(ents.items(), key=lambda kv: (auth_key(kv[1]), int(kv[1].get("year") or 0), (kv[1].get("title") or "").lower(), kv[0]))
    groups = defaultdict(list)
    for k, e in keyed:
        n = len(e["authors"])
        g = (fold(e["authors"][0][0]) if e["authors"] else "", int(e.get("year") or 0), min(n, 3), fold(e["authors"][1][0]) if n == 2 else "")
        groups[g].append(k)
    suffix = {}
    for g, ks in groups.items():
        if len(ks) > 1:
            for i, k in enumerate(ks):
                suffix[k] = chr(ord("a") + i) if i < 26 else "z%d" % i
    return keyed, suffix


def md_escape(t):
    return re.sub(r"([\\`*_{}\[\]<>#|])", r"\\\1", t).replace("&", "&amp;") if False else re.sub(r"([\\`*_\[\]<>|])", r"\\\1", t)


def paper_anchor_md(sid):
    """The ID of a source and the anchor `#paper-<id>` that the Papers page of the References section uses for the same source."""
    return '<a id="paper-%s"></a>**%s** ' % (sid.lower(), sid) if sid else ""


def md_entry(parts):
    line = "%s %s" % (md_escape(parts["authors"]), md_escape(parts["year"]))
    if parts["title"]:
        line += ", " + md_escape(parts["title"])
    for b in parts["bits"]:
        line += ", " + md_escape(b)
    for i in parts["ids"]:
        line += ", " + md_escape(i)
    line = line.rstrip(".") + "."
    links = []
    if parts["url"]:
        links.append("[%s](%s)" % ("paper" if parts["kind"] != "dataset" else "data", parts["url"]))
    if parts["ads"]:
        links.append("[ADS](%s)" % parts["ads"])
    if links:
        line += " " + " · ".join(links)
    return "- " + paper_anchor_md(parts["sid"]) + line


def html_entry(parts):
    h = html.escape
    line = "%s %s" % (h(parts["authors"]), h(parts["year"]))
    if parts["title"]:
        line += ", <i>" + h(parts["title"]) + "</i>"
    for b in parts["bits"]:
        line += ", " + h(b)
    for i in parts["ids"]:
        if i.startswith("doi:"):
            line += ', <a href="https://doi.org/%s">%s</a>' % (h(urllib.parse.quote(i[4:], safe="/()")), h(i))
        else:
            line += ', <a href="https://arxiv.org/abs/%s">%s</a>' % (h(i[6:]), h(i))
    line = line.rstrip(".") + "."
    links = []
    if parts["url"]:
        links.append('<a class="u" href="%s">%s</a>' % (h(parts["url"]), "paper" if parts["kind"] != "dataset" else "data"))
    if parts["ads"]:
        links.append('<a class="u" href="%s">ADS</a>' % h(parts["ads"]))
    return '<li id="paper-%s"><b class="rid">%s</b> %s <span class="links">%s</span></li>' % (parts["sid"].lower(), parts["sid"], line, " &middot; ".join(links))


INTRO = (
    "This project would not be possible without the dedicated work of the researchers, observers, instrument teams and data archives "
    "whose papers and datasets it builds on. Every number the generator draws from physics, every law of how stars, planets, moons, "
    "belts and atmospheres come to be, was first measured, modelled, checked and written down by people who spent years, often a whole "
    "career, on a single question and then shared the answer openly so that anyone could build on it. The telescopes and spacecraft "
    "teams, the laboratory spectroscopists, the survey collaborations, the modellers and the volunteers who keep archives alive are why a "
    "made-up galaxy can still obey the same physics as the real one. We are deeply grateful to every one of them, and this page is our "
    "attempt to name the work we stand on."
)

CONVENTIONS = (
    "References follow the AAS journal and NASA Astrophysics Data System (ADS) style: author-year, alphabetical by first author and then "
    "by the later authors and the year, with the journal in its ADS abbreviation, then volume and first page (or article number). Up to "
    "five authors are named; longer author lists shorten to the first author and \"et al.\" (the full list is on the linked page). Where "
    "one first author has more than one paper in a year, the years carry the letters a, b, c. The title is printed in italics, which the "
    "AAS journals now permit. Each entry gives the DOI, the arXiv identifier when the paper has one, and two links: the DOI or arXiv "
    "page, and the ADS abstract page. Metadata come from Crossref, DataCite and the arXiv API; nothing was filled in from memory."
)


def datasets_for_render(c):
    ids = c["ids"]
    rows = []
    for a in c.get("archives", []):
        rows.append((a["name"], a["provider"], a["url"], "", ids["archive:" + a["url"].lower()]))
    for d in c.get("setup_datasets", []):
        rows.append((d["title"], d["source"], d["url"], "downloaded at setup", ids["setup:" + d["url"].lower()]))
    seen, out = set(), []
    for r in sorted(rows, key=lambda r: r[0].lower()):
        if r[2] not in seen or r[3]:
            out.append(r)
            seen.add(r[2])
    return out


def _slug(heading):
    heading = re.sub(r"[*~]", "", heading).strip().lower()
    heading = re.sub(r"[^\w\- ]", "", heading)
    return heading.replace(" ", "-")


def render_all(c):
    ids = c["ids"]
    ents = {k: e for k, e in c["entries"].items() if e.get("authors") or e.get("title")}
    papers = {k: e for k, e in ents.items() if e.get("type") != "dataset"}
    datas = {k: e for k, e in ents.items() if e.get("type") == "dataset"}
    keyed, suffix = sorted_entries(papers)
    dkeyed, dsuffix = sorted_entries(datas)
    groups = citation_groups(c)
    npr, rvr = not_paper_rows(groups), review_rows(groups)
    places = {g: sum(n for _, n, _ in groups[g]) for g in groups}
    letters = defaultdict(list)
    for k, e in keyed:
        f = fold(e["authors"][0][0])[:1].upper() if e["authors"] else "0"
        letters[f or "0"].append((k, e))
    names = {"0": "No author listed (publisher or organisation)"}
    stats = {"papers": len(papers), "datasets": len(datas), "archives": len(datasets_for_render(c)), "tied": len(groups["tied"]), "not": len(groups["not"]), "review": len(groups["review"]),
             "bad_ids": len(c.get("unresolved_ids", {}))}
    # ---- markdown
    md = ["# References", "", INTRO, "", "## How to read this list", "", CONVENTIONS, "",
          "This page is written by `scripts/build_references.py` from `docs/reference/references-data.json`; the HTML twin "
          "[REFERENCES.html](REFERENCES.html) is built from the same data. Do not edit either by hand.", "",
          "**%d papers and books, %d datasets with their own DOI, and %d archives, catalogs and tools. Of the %s author-year citations in the "
          "notes, reports and code (%s places), %s are tied to a paper above, %s are labelled as not a paper, and %s are listed for review.**" %
          (stats["papers"], stats["datasets"], stats["archives"], format(sum(stats[k] for k in ("tied", "not", "review")), ","),
           format(sum(places[g] for g in ("tied", "not", "review")), ","), format(stats["tied"], ","), format(stats["not"], ","), format(stats["review"], ",")), "",
          "## Contents", "", "- [Papers and books](#papers-and-books) (" + ", ".join("[%s](#%s)" % (L if L != "0" else "other", _slug(names.get(L, L))) for L in sorted(letters)) + ")",
          "- [Datasets, catalogs, archives and tools](#datasets-catalogs-archives-and-tools)",
          "- [Author-year citations: where each one stands](#author-year-citations-where-each-one-stands)", "", "## Papers and books", ""]
    for L in sorted(letters):
        md += ["### %s" % names.get(L, L), ""]
        for k, e in letters[L]:
            md.append(md_entry(entry_parts(e, suffix.get(k, ""), ids[k.lower()])))
        md.append("")
    md += ["## Datasets, catalogs, archives and tools", "", "Archives and catalogs the project reads, with the page each is served from. A dataset the setup "
           "step downloads says so; the exact files and checksums are in `src/planetgen/datasets/catalogue.py`.", ""]
    for name, prov, url, note, sid in datasets_for_render(c):
        md.append("- %s%s. %s%s. <%s>" % (paper_anchor_md(sid), md_escape(name), md_escape(prov), (" (" + note + ")") if note else "", url))
    md.append("")
    if dkeyed:
        md += ["### Published datasets (with a DOI)", ""]
        for k, e in dkeyed:
            md.append(md_entry(entry_parts(e, dsuffix.get(k, ""), ids[k.lower()])))
        md.append("")
    md += ["## Author-year citations: where each one stands", "",
           "The design notes, reports and code cite many sources only as 'Author year'. Each of those citations is in exactly one of three "
           "groups. **Tied to a paper:** %s citations (%s places) are matched to a paper in the list above, by an identifier beside the mention, "
           "by the volume and page, or by the first author, year and words of the citing sentence; nothing is tied on a guess. **Not a paper:** "
           "%s are time stamps, organisations, standards, web pages, programs or sentence words, listed below with what each is and a link. "
           "**For review:** %s (%s places) could not be tied with confidence; each is listed with what was tried and a link that searches ADS "
           "for it." % (format(stats["tied"], ","), format(places["tied"], ","), format(stats["not"], ","), format(stats["review"], ","),
                        format(places["review"], ",")), "", "### Not a paper", ""]
    for a_, ys, n_, what, url, cls in npr:
        md.append("- %s %s (%s place%s): %s%s" % (md_escape(a_), years_text(ys), format(n_, ","), "" if n_ == 1 else "s", md_escape(what), (". <%s>" % url) if url else ""))
    md += ["", "### For review", ""]
    for text, n_, why, near, ads, tried in rvr:
        md.append("- %s (%d place%s). Tried: %s. %s%s ADS: <%s>" % (md_escape(text), n_, "" if n_ == 1 else "s", ", ".join(tried) or "nothing yet", md_escape(why) + ".",
                                                                 (" Nearest: " + md_escape(near_text(near)) + ".") if near else "", ads))
    md.append("")
    bad = sorted(c.get("unresolved_ids", {}).items())
    if bad:
        md += ["### Identifiers that did not resolve", ""]
        for k, why in bad:
            md.append("- %s: %s" % (md_escape(k), md_escape(why)))
        md.append("")
    # ---- html
    h = html.escape
    out = ['<!doctype html>', '<html lang="en">', '<head>', '<meta charset="utf-8">', '<meta name="viewport" content="width=device-width, initial-scale=1">',
           '<title>PlanetGen references</title>', '<meta name="description" content="Every paper, book, dataset and archive the PlanetGen design and generator cite, in AAS and ADS style, for anyone checking a source.">', '<style>', HTML_CSS, '</style>', '</head>', '<body><main>', '<p><a href="index.html">Docs index</a></p>', '<h1>References</h1>',
           '<p class="intro">%s</p>' % h(INTRO), '<h2>How to read this list</h2>', '<p>%s</p>' % h(CONVENTIONS),
           '<p>This page is written by <code>scripts/build_references.py</code> from <code>docs/reference/references-data.json</code>; its Markdown twin '
           '<a href="REFERENCES.md">REFERENCES.md</a> is built from the same data. Do not edit either by hand.</p>',
           '<p class="count"><b>%d papers and books, %d datasets with their own DOI, and %d archives, catalogs and tools. Of the %s author-year citations in the '
           'notes, reports and code (%s places), %s are tied to a paper above, %s are labelled as not a paper, and %s are listed for review.</b></p>' %
           (stats["papers"], stats["datasets"], stats["archives"], format(sum(stats[k] for k in ("tied", "not", "review")), ","),
            format(sum(places[g] for g in ("tied", "not", "review")), ","), format(stats["tied"], ","), format(stats["not"], ","), format(stats["review"], ",")),
           '<nav><b>Papers and books:</b> ' + " ".join('<a href="#L-%s">%s</a>' % (L, L) for L in sorted(letters)) +
           ' &middot; <a href="#datasets">Datasets and archives</a> &middot; <a href="#authoryear">Author-year citations</a></nav>', '<h2 id="papers">Papers and books</h2>']
    for L in sorted(letters):
        out += ['<h3 id="L-%s">%s</h3>' % (L, names.get(L, L)), "<ul>"]
        for k, e in letters[L]:
            out.append(html_entry(entry_parts(e, suffix.get(k, ""), ids[k.lower()])))
        out.append("</ul>")
    out += ['<h2 id="datasets">Datasets, catalogs, archives and tools</h2>', '<p>Archives and catalogs the project reads, with the page each is served from. '
            'A dataset the setup step downloads says so; the exact files and checksums are in <code>src/planetgen/datasets/catalogue.py</code>.</p>', "<ul>"]
    for name, prov, url, note, sid in datasets_for_render(c):
        out.append('<li id="paper-%s"><b class="rid">%s</b> %s. %s%s. <a href="%s">%s</a></li>' % (sid.lower(), sid, h(name), h(prov), h(" (" + note + ")") if note else "", h(url), h(url)))
    out.append("</ul>")
    if dkeyed:
        out += ["<h3>Published datasets (with a DOI)</h3>", "<ul>"]
        for k, e in dkeyed:
            out.append(html_entry(entry_parts(e, dsuffix.get(k, ""), ids[k.lower()])))
        out.append("</ul>")
    out += ['<h2 id="authoryear">Author-year citations: where each one stands</h2>',
            "<p>The design notes, reports and code cite many sources only as &lsquo;Author year&rsquo;. Each of those citations is in exactly one of three groups. "
            "<b>Tied to a paper:</b> %s citations (%s places) are matched to a paper in the list above, by an identifier beside the mention, by the volume and page, "
            "or by the first author, year and words of the citing sentence; nothing is tied on a guess. <b>Not a paper:</b> %s are time stamps, organisations, standards, "
            "web pages, programs or sentence words, listed below with what each is and a link. <b>For review:</b> %s (%s places) could not be tied with confidence; "
            "each is listed with what was tried and a link that searches ADS for it.</p>" % (format(stats["tied"], ","), format(places["tied"], ","),
            format(stats["not"], ","), format(stats["review"], ","), format(places["review"], ",")), "<h3>Not a paper</h3>", "<ul>"]
    for a_, ys, n_, what, url, cls in npr:
        out.append("<li>%s %s (%s place%s): %s%s</li>" % (h(a_), years_text(ys), format(n_, ","), "" if n_ == 1 else "s", h(what),
                                                         (' <a href="%s">%s</a>' % (h(url), h(url))) if url else ""))
    out += ["</ul>", "<h3>For review</h3>", "<ul>"]
    for text, n_, why, near, ads, tried in rvr:
        out.append('<li>%s (%d place%s). Tried: %s. %s%s <a href="%s">Search ADS</a></li>' % (
            h(text), n_, "" if n_ == 1 else "s", h(", ".join(tried) or "nothing yet"), h(why) + ".", (" Nearest: " + h(near_text(near)) + ".") if near else "", h(ads)))
    out.append("</ul>")
    if bad:
        out += ["<h3>Identifiers that did not resolve</h3>", "<ul>"] + ["<li>%s: %s</li>" % (h(k), h(w)) for k, w in bad] + ["</ul>"]
    out += ["</main></body>", "</html>", ""]
    return "\n".join(md), "\n".join(out)


HTML_CSS = """
:root{--bg:#fff;--fg:#1d2330;--mut:#5b6474;--link:#1a56a8;--rule:#d9dde5;--card:#f6f7fa}
@media (prefers-color-scheme:dark){:root{--bg:#12151c;--fg:#e4e7ee;--mut:#9aa3b5;--link:#7db2ff;--rule:#2b3140;--card:#1a1f2a}}
html{font-family:Georgia,'Times New Roman',serif;line-height:1.5}
body{background:var(--bg);color:var(--fg);margin:0}
main{max-width:62rem;margin:0 auto;padding:1.5rem 16px 4rem}
h1{font-size:2rem;margin:.2em 0 .5em}h2{margin-top:2.2em;border-bottom:1px solid var(--rule);padding-bottom:.2em}h3{margin-top:1.6em;color:var(--mut)}
.intro{font-size:1.1rem;background:var(--card);padding:1em 1.2em;border-left:4px solid var(--link);border-radius:4px}
nav{display:flex;flex-wrap:wrap;gap:.4em .7em;margin:1em 0;font-family:system-ui,sans-serif;font-size:.9rem}
ul{list-style:none;padding:0}li{margin:.55em 0;padding-left:2em;text-indent:-2em;font-size:.95rem;overflow-wrap:anywhere}
a{color:var(--link)}.rid{font-family:ui-monospace,monospace;font-size:.8rem;color:var(--mut)}.links{white-space:nowrap;font-family:system-ui,sans-serif;font-size:.8rem}
code{font-size:.9em}
"""


def ay_text(key):
    """'Hamilton and Burns 1991' for the key 'Hamilton|Burns|pair|1991'."""
    a, b, kind, y = key.split("|")
    return "%s%s %s" % (a, " and " + b if kind == "pair" else (" et al." if kind == "etal" else ""), y)


def citation_groups(c):
    """Every author-year mention the last harvest found, in three groups: tied to a paper in the list, labelled as not a paper (what it
    is, a link), and still for review (what was tried). A mention with none of the three is returned in 'open' and fails the test."""
    st = c.get("authoryear_status", {})
    out = {"tied": [], "not": [], "review": [], "open": []}
    for key, v in c.get("candidates", {}).get("authoryear", {}).items():
        s = st.get(key, {})
        n = v["n"]
        if s.get("matched"):
            out["tied"].append((key, n, s))
        elif s.get("not_a_paper"):
            out["not"].append((key, n, s))
        elif s.get("review"):
            out["review"].append((key, n, s))
        else:
            out["open"].append((key, n, s))
    return out


def ads_search_url(key):
    a, b, kind, y = key.split("|")
    q = 'author:"^%s" year:%s' % (a, y)
    return "https://ui.adsabs.harvard.edu/search/q=" + urllib.parse.quote(q)


def not_paper_rows(groups):
    """One row per name: (name, years, places, what, url, class), grouped so '5,874 time stamps' is one line."""
    by = defaultdict(lambda: {"years": set(), "n": 0, "what": "", "url": "", "class": ""})
    for key, n, s in groups["not"]:
        a, b, kind, y = key.split("|")
        r = by[a]
        r["years"].add(int(y))
        r["n"] += n
        r["what"], r["url"], r["class"] = s["not_a_paper"]["what"], s["not_a_paper"].get("url", ""), s["not_a_paper"].get("class", "")
    return sorted(((a, sorted(r["years"]), r["n"], r["what"], r["url"], r["class"]) for a, r in by.items()), key=lambda r: (r[5], fold(r[0])))


def review_rows(groups):
    """The mentions still for review, heaviest first: (text, places, why not tied, nearest records, ADS search link)."""
    rows = []
    for key, n, s in groups["review"]:
        r = s["review"]
        rows.append((ay_text(key), n, r.get("why", ""), r.get("near", []), ads_search_url(key), s.get("tried", [])))
    return sorted(rows, key=lambda x: (-x[1], fold(x[0])))


def years_text(ys):
    return ", ".join(str(y) for y in ys) if len(ys) <= 6 else "%d to %d" % (ys[0], ys[-1])


def near_text(near):
    """The nearest records Crossref returned, as 'Title (Journal vol page, year, doi:...)'."""
    out = []
    for x in near[:2]:
        where = " ".join(str(t) for t in (x.get("journal"), x.get("volume"), x.get("page")) if t)
        out.append("%s (%s%s; doi:%s)" % (x.get("title", "").rstrip("."), where + ", " if where else "", x.get("year", ""), x.get("doi", "")))
    return "; ".join(out)


def cmd_render(args):
    c = load_cache()
    given = assign_ids(c)
    if args.check and given:
        sys.exit("%d sources have no ID (run python scripts/build_references.py render)" % given)
    md, ht = render_all(c)
    if given:
        save_cache(c)
        log("gave IDs to", given, "sources")
    if args.check:
        bad = [p for p, t in ((OUT_MD, md), (OUT_HTML, ht)) if not os.path.exists(p) or open(p).read() != t]
        if bad:
            sys.exit("out of date: %s (run python scripts/build_references.py render)" % ", ".join(os.path.basename(b) for b in bad))
        print("references pages match the cache")
        return
    with open(OUT_MD, "w") as fh:
        fh.write(md)
    with open(OUT_HTML, "w") as fh:
        fh.write(ht)
    log("wrote", OUT_MD, OUT_HTML)


def cmd_match(args):
    """Tie author-year mentions to resolved entries by first author, year and (for pairs) second author. The decisions made later by
    scripts/settle_references.py (searches, labels, review notes) are kept; mentions the last harvest no longer finds are dropped."""
    c = load_cache()
    base = match_authoryear(c)
    old = c.get("authoryear_status", {})
    status = {}
    for k, v in base.items():
        s = dict(old.get(k, {}))
        s["n"] = v["n"]
        if v["matched"]:
            s["matched"] = v["matched"]
        s.setdefault("matched", [])
        status[k] = s
    # `resolve` merges an arXiv entry into the entry for its published DOI; follow a tie to the entry that now carries the paper
    ents = c["entries"]
    by_id = {}
    for k, e in ents.items():
        for i in (e.get("doi"), e.get("arxiv")):
            if i:
                by_id[i.lower()] = k
    for s in status.values():
        fixed = []
        for ref in s.get("matched", []):
            ref = ref if ref in ents else by_id.get(ref.lower().replace("arxiv:", ""), by_id.get(ref.lower()))
            if ref and ref not in fixed:
                fixed.append(ref)
        s["matched"] = fixed
    c["authoryear_status"] = status
    save_cache(c)
    log("author-year keys:", len(status), "tied:", sum(1 for v in status.values() if v["matched"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for n in ("harvest",):
        p = sub.add_parser(n)
        p.add_argument("--extra", nargs="*", default=[])
    sub.add_parser("match")
    sub.add_parser("render").add_argument("--check", action="store_true")
    sub.add_parser("resolve").add_argument("--offline", action="store_true")
    a = ap.parse_args()
    if a.cmd == "harvest":
        cmd_harvest(a)
    elif a.cmd == "resolve":
        cmd_resolve(a)
    elif a.cmd == "match":
        cmd_match(a)
    elif a.cmd == "render":
        cmd_render(a)
