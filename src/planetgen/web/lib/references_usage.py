# planetgen/web/lib/references_usage.py

"""
Where we use each source (DOC.89, References R3): which code and data files, stand-ins, and design notes cite it, and in which area.

A text cites a source as its ID in square brackets (`[R0123]`, DOC.98). This module reads those brackets from the files of the project and gives,
for every source ID, the places that cite it. The Papers page shows them (`references_site.papers_page`), counts them ("most cited by us"), and
groups the sources by area.

The result is computed when the page is first asked for and kept until the process ends, not committed. Every text that cites a source would
otherwise change a built file, and two lanes would collide on it for no reason. A test reads the same function, so the page and the check cannot
differ.

Three kinds of place:

* `values`: a code or data file (`src/planetgen/...`, `scripts/...`) that cites the source;
* `stand_ins`: a row of the stand-in registry (`docs/design/research/gap-filling-plan.md`, GEN.725) that cites it, by its registry key (`5.1:GR1`, the section number
  and the row ID, which is also the key of its page, DOC.90);
* `notes`: a design note (`docs/design/...`) that cites it.

The area of a place is read from its name: the first rule of `RULES` that matches the path of a file, or the section title of a registry row. A
path no rule matches is in the area `other`. A source cited in two areas is in both. Standard library only.
"""

import functools
import os
import re

from . import citations
from .reference_urls import AREAS

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", ".."))
"""str: The top folder of the project (the folder that holds `docs/` and `src/`)."""

OTHER = "other"
"""str: The area of a place that no rule of `RULES` names."""

AREA_TITLES = {
    "stars": "Stars", "galaxy": "The galaxy", "orbits": "Orbits and formation", "bodies": "Rocky and giant bodies", "moons-belts": "Moons and belts",
    "atmospheres": "Atmospheres", "climate": "Climate", "radiation": "Radiation", "rings-magnetism": "Rings and magnetism", "life": "Life",
    "appearance": "Color and maps", OTHER: "Other or general",
}
"""dict: The name of each area on the Papers page; the areas are those of the method pages (`reference_urls.AREAS`) and `other`."""

RULES = (
    ("life", r"life|habitab|biogenic|xenobio|survival|adaptation|protection|toxicity|thermal-respir|environmental-extremes|chemical-habit|hycean"),
    ("rings-magnetism", r"ring|magnet|dynamo"),
    ("radiation", r"radiation|radiogenic|radioisotope|solar_particles|ionizing|shielding|potassium|birth_radio"),
    ("appearance", r"colou?r|albedo_color|cloud-albedo|sky|photometry|mie|rayleigh|spectrum|spectra|stellar_library|render|3d-|map-ui|text-style|"
                   r"fly-through|appearance|darkening|describe|ozone_uv|chappuis"),
    ("climate", r"climate|saturation|seasonal|radiative|two_stream|greenhouse|infrared|line[-_]opacity|co2_far|hitran|ir_ktables|cia|thermosphere|"
                r"day_night|daylight|nightside|ocean[-_]tides|hydrosphere|aquifer|water_layer|surface[-_]liquid|surface_air|band_winds|ice_types|"
                r"clathrate|cloud|sulfuric|acid|aerosol|brewer|sulfur|waterworld|strange_solvents|substellar"),
    ("atmospheres", r"atmosph|gases|gas[-_]table|air_column|photolysis|escape|retention|mixture|ozone|surface[-_]pressure|supply|species_thermo|"
                    r"nasa_glenn|thermo|nh3|helium_fractionation|evolved|volatile|delivered|delivery|equilibrium-chem|condensate|envelope"),
    ("moons-belts", r"moon|belt|asteroid|comet|trojan|lagrange|irregular|captured|debris|small[-_]bod|meteorite|kuiper|taxon|interstellar-(object|aster)|"
                    r"field_|wanderer"),
    ("orbits", r"orbit|tidal|tide|resonan|stability|quadruple|hierarchy|collision|spin|rotation|course|closest|kepler|accretion|formation|"
               r"architecture|clearance|roche|scatter|impact|pair_history|mergers|sky-view|integrator|multi-?star|binary|mass_ratio|natal|rogue|"
               r"disintegrat|super_puff|planet[-_]mass|mass[-_]radius"),
    ("galaxy", r"galax|bar-|gen533|halo|globular|cluster|king|potential|star_formation|disk|nebula|superbubble|hii|black[-_]hole|explosion|"
               r"compact[-_]remnant|remnant|interstellar|core|intergalactic|imf|stellar-mass|stellar_mass|age[_-]metal|alpha|ms_lifetime|harris|thick"),
    ("stars", r"star|stellar|main_sequence|brown|white[-_]dwarf|neutron|isochrone|agb|carbon|pre_main|post_main|post-main|mass_loss|hosts|"
              r"helium_core|chemically|rprocess|mamajek|baraffe|phenomena|sequence|activity|mass_boundaries|wd[-_]"),
    ("bodies", r"bulk|interior|composition|planet|giant|rock|core_heat|crust|mantle|magma|materials|elements|isotopes|decay|compaction|roundness|"
               r"binding|strength|radius|body|surface|dry_carbon|oxygen_fug|redox|iron|mass|water|hydro|ice|salt|brine|lava|jwst|tuning|dwarf|size"),
)
"""tuple: `(area, pattern)` in the order tried, matched against the lower-case path of a file or the section title of a registry row."""

assert all(area in AREAS for area, _pattern in RULES)

SCAN_FOLDERS = ("src/planetgen", "docs/design", "scripts")
"""tuple: The folders whose files are read (relative to `ROOT`)."""

SKIP_FOLDERS = frozenset({"__pycache__", "static", "vendor", "migrations", "node_modules", "tests", "generated"})
"""frozenset: Folder names that are never read."""

SKIP_FILES = frozenset({
    # these files write source IDs as examples or build the pages that show them
    "citations.py", "references_site.py", "references_usage.py", "references_pages.py", "build_reference_pages.py", "build_references.py",
    "build_stand_in_values.py", "rewrite_citations.py", "rewrite_bibliographies.py", "references.md",
    # the registry is read as rows (`registry_rows`), and its generated twin repeats them
    "gap-filling-plan.md", "stand-in-values.md",
})
"""frozenset: File names that are never read as a place of use."""

TEXT_SUFFIXES = (".py", ".json", ".md", ".html", ".txt")
REGISTRY = os.path.join("docs", "design", "research", "gap-filling-plan.md")
_HEADING = re.compile(r"^#{2,4}\s+(.*\S)")
_RULE = re.compile(r"^\|\s*:?-{3,}")
_ROW_ID = re.compile(r"^[A-Z]{1,4}(\.[a-z])?\d+[a-z]?(, [a-z]?\d+)*$")


def area_of(text):
    """The area a path or a section title belongs to: the first matching rule of `RULES`, else `other`."""
    lower = text.lower()
    for area, pattern in RULES:
        if re.search(pattern, lower):
            return area
    return OTHER


def _kind(path):
    return "notes" if path.startswith("docs/") else "values"


def _files(root):
    for folder in SCAN_FOLDERS:
        for here, dirs, names in os.walk(os.path.join(root, folder)):
            dirs[:] = sorted(name for name in dirs if name not in SKIP_FOLDERS)
            for name in sorted(names):
                if name.endswith(TEXT_SUFFIXES) and name not in SKIP_FILES:
                    path = os.path.join(here, name)
                    yield os.path.relpath(path, root).replace(os.sep, "/"), path


def registry_rows(text):
    """
    `[(row id, section title, ids cited)]` for every row of the tables of the stand-in registry text (GEN.725). A row is a table line under a
    heading; the first cell is its ID (`GR1`, `BC12`, `AT.a1, a2`; a table whose first column is something else is not the registry); a header line (the one above the `---` line) is skipped.
    """
    rows, heading, lines = [], "", text.splitlines()
    for number, line in enumerate(lines):
        match = _HEADING.match(line)
        if match:
            heading = match.group(1)
            continue
        if not line.startswith("|") or _RULE.match(line):
            continue
        if number + 1 < len(lines) and _RULE.match(lines[number + 1]):
            continue
        cited = citations.ids_in(line)
        if cited:
            first = line.strip().strip("|").split("|")[0].replace("`", "").strip()
            if _ROW_ID.match(first):
                rows.append((first, heading, cited))
    return rows


def scan(root=None):
    """
    `{"uses": {id: {"values": [path], "notes": [path], "stand_ins": [section:row id], "areas": [area]}}, "files": number read}`: where every source ID
    that some text cites is cited. The lists are sorted, once each.
    """
    root = root or ROOT
    uses, files = {}, 0

    def place(sid, kind, name, area):
        entry = uses.setdefault(sid, {"values": [], "notes": [], "stand_ins": [], "areas": set()})
        if name not in entry[kind]:
            entry[kind].append(name)
        entry["areas"].add(area)

    for relative, path in _files(root):
        with open(path, encoding="utf-8", errors="replace") as handle:
            text = handle.read()
        files += 1
        area = area_of(relative)
        for sid in citations.ids_in(text):
            place(sid, _kind(relative), relative, area)
    registry = os.path.join(root, REGISTRY)
    if os.path.isfile(registry):
        with open(registry, encoding="utf-8") as handle:
            for row, section, cited in registry_rows(handle.read()):
                number = re.match(r"\d+\.\d+a?\b", section)
                key = f"{number.group(0)}:{row}" if number else row
                for sid in cited:
                    place(sid, "stand_ins", key, area_of(section))
    for entry in uses.values():
        for kind in ("values", "notes", "stand_ins"):
            entry[kind].sort()
        entry["areas"] = sorted(entry["areas"], key=lambda area: (area == OTHER, area))
    return {"uses": uses, "files": files}


@functools.lru_cache(maxsize=1)
def _cached():
    return scan()


def uses():
    """`scan()` for this project, read once and kept (a new process reads the files again)."""
    return _cached()


def cited_by(entry):
    """How many places cite a source: its files, registry rows, and notes together."""
    return len(entry["values"]) + len(entry["notes"]) + len(entry["stand_ins"])
