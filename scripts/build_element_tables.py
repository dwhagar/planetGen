#!/usr/bin/env python3
"""
Builds the element and isotope data files (GEN.438) from the published tables they come from.

    python scripts/build_element_tables.py                # fetch the sources and write the data files
    python scripts/build_element_tables.py --from DIR     # read the sources from DIR (names below), no network
    python scripts/build_element_tables.py --check        # build in memory and fail if the committed files differ

The data files are `src/planetgen/physics/data/elements.json`, `isotopes.json` and `decay_series.json`; each states its
sources and citations. The raw tables are not committed and are never redistributed; only the values read from
them are. Sources, all read on the day the files were built:

* `ciaaw-standard.html` and `ciaaw-abridged.html`: CIAAW, Standard Atomic Weights and Abridged Standard Atomic Weights
  (https://ciaaw.org/atomic-weights.htm, https://www.ciaaw.org/abridged-atomic-weights.htm), the "Atomic Weights 2021"
  report with the 2024 revisions of Gd, Lu and Zr (Prohaska et al., Pure Appl. Chem. 94, 573, 2022,
  [R3615]).
* `nist-compositions.txt`: NIST, Atomic Weights and Isotopic Compositions with Relative Atomic Masses (the ASCII listing
  of https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl, `isotype=some` for the natural composition and
  `isotype=all` for every isotope's mass); the isotopic compositions are IUPAC's, the masses the atomic mass evaluation's.
* `nubase2020.txt`: NUBASE2020 (Kondev et al., Chin. Phys. C 45, 030001, 2021, [R2443]), from
  https://www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt: half-lives, decay modes and branchings, and the mass
  excesses of the atomic mass evaluation AME2020 (Wang et al., Chin. Phys. C 45, 030003, 2021).
* The decay energies of the radionuclides come from Lugaro, Ott and Kereszturi, Prog. Part. Nucl. Phys. 102, 1 (2018),
  [R2768: Table 4] [R2768], and the heat constants of the four long-lived parents
  from Dye, Rev. Geophys. 50, RG3007 (2012), [R1203: Table 2] [R1203], typed in below.
"""

import argparse
import html
import json
import math
import os
import re
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))  # the standard User-Agent comes from the package (API.31)

from planetgen.util.user_agent import user_agent  # noqa: E402
DATA = os.path.join(HERE, "..", "src", "planetgen", "physics", "data")
RETRIEVED = "2026-10-17"

URLS = {
    "ciaaw-standard.html": "https://ciaaw.org/atomic-weights.htm",
    "ciaaw-abridged.html": "https://www.ciaaw.org/abridged-atomic-weights.htm",
    "nist-compositions.txt": "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=all",
    "nist-natural.txt": "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some",
    "nubase2020.txt": "https://www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt",
}

SECONDS = {"ys": 1e-24, "zs": 1e-21, "as": 1e-18, "fs": 1e-15, "ps": 1e-12, "ns": 1e-9, "us": 1e-6, "ms": 1e-3,
           "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0, "y": 365.25 * 86400.0}
for _prefix, _power in (("k", 3), ("M", 6), ("G", 9), ("T", 12), ("P", 15), ("E", 18), ("Z", 21), ("Y", 24)):
    SECONDS[_prefix + "y"] = SECONDS["y"] * 10.0 ** _power
"""dict: NUBASE2020's half-life units in seconds (the year is the Julian year of 365.25 days)."""

MEV_PER_KEV = 1e-3
ERG_J = 1e-7
KEV_J = 1.602176634e-16
"""float: One keV in joules (the exact value of the elementary charge)."""
HE4_MASS_EXCESS_KEV = None
"""float: Filled from NUBASE2020's 4He line."""

MIN_BRANCH_PERCENT = 0.005
"""float: A decay branch below this share is left out of the chains, as Lugaro et al. leave out branchings below one per
mil (their section 6.1); the threshold here is lower."""

NATURAL_SERIES = {"U-238": "Pb-206", "U-235": "Pb-207", "Th-232": "Pb-208"}
"""dict: The long-lived parents whose decay chain is followed to the stable lead isotope they end in."""

RADIONUCLIDES = (
    "Be-7", "Be-10", "Al-26", "Cl-36", "K-40", "Ca-41", "Mn-53", "Fe-60", "Co-60", "Rb-87", "Nb-92", "Pd-107", "I-129",
    "Sm-146", "Sm-147", "Hf-182", "Pb-205", "Lu-176", "Re-187", "Pu-244", "Cm-247", "Rn-220", "Rn-222")
"""tuple: The radioactive isotopes outside the three chains the table holds (the short-lived radionuclides of Lugaro et al.
Table 4, the long-lived ones it lists, the iron-60 daughter, and the two radon isotopes the ground's dose reads)."""

LUGARO_ERG = {
    "Al-26": 5.00e-6, "Be-10": 3.24e-7, "Mn-53": 8.23e-9, "Pd-107": 1.49e-8, "Hf-182": 2.96e-6, "I-129": 1.23e-7,
    "Cm-247": 2.67e-5, "Nb-92": 2.03e-6, "Sm-146": 4.05e-6, "Cl-36": 3.95e-7, "Fe-60": 4.16e-6, "Pu-244": 2.50e-5,
    "Be-7": 7.99e-8, "Ca-41": 4.31e-9, "Pb-205": 8.35e-9, "K-40": 1.05e-6, "Rb-87": 1.31e-7, "Sm-147": 3.70e-6,
    "Lu-176": 1.24e-6, "Re-187": 9.90e-10, "Th-232": 6.46e-5, "U-235": 7.10e-5, "U-238": 7.64e-5}
"""dict: Energy released per decay (erg), Table 4 of [R2768]: the full Q value for alpha decays, the electrons'
and radiation's energies without the neutrinos for beta decays; for a parent with a chain, the whole chain to the stable
end."""

DYE_2012 = {
    "U-238": {"abundance_percent": 99.2796, "decay_constant_per_s": 4.916e-18, "neutrinos": 6, "q_pj": 8.282,
              "q_neutrino_pj": 0.634, "q_heat_pj": 7.648, "heat_uw_per_kg": 95.13},
    "U-235": {"abundance_percent": 0.7204, "decay_constant_per_s": 31.210e-18, "neutrinos": 4, "q_pj": 7.434,
              "q_neutrino_pj": 0.325, "q_heat_pj": 7.108, "heat_uw_per_kg": 568.47},
    "Th-232": {"abundance_percent": 100.0, "decay_constant_per_s": 1.563e-18, "neutrinos": 4, "q_pj": 6.833,
               "q_neutrino_pj": 0.358, "q_heat_pj": 6.475, "heat_uw_per_kg": 26.28},
    "K-40": {"abundance_percent": 0.0117, "decay_constant_per_s": 17.200e-18, "neutrinos": 0.893, "q_pj": 0.213,
             "q_neutrino_pj": 0.103, "q_heat_pj": 0.110, "heat_uw_per_kg": 28.47},
}
"""dict: [R1203: Table 2]: the four parents of the Earth's radiogenic heat, per decay (the whole chain for U and Th)."""


def _clean(value):
    """The value with every float rounded to 12 significant digits, so a product like 1.05e-6 * 1e-7 stays readable."""
    if isinstance(value, float):
        return float(f"{value:.12g}")
    if isinstance(value, dict):
        return {key: _clean(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_clean(item) for item in value]
    return value


def _fetch(name, source_dir):
    path = os.path.join(source_dir, name) if source_dir else None
    if path and os.path.exists(path):
        with open(path, encoding="utf-8", errors="replace") as handle:
            return handle.read()
    request = urllib.request.Request(URLS[name], headers={"User-Agent": user_agent()})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8", errors="replace")


def _number(text):
    """`1.00782503223(9)` -> (1.00782503223, 9e-11); a plain number has no uncertainty; empty is `(None, None)`."""
    text = text.strip().replace(" ", "").replace("#", "")
    if not text:
        return None, None
    match = re.fullmatch(r"([0-9.]+(?:[eE][+-]?\d+)?)(?:\((\d+)\))?", text)
    if not match:
        return None, None
    value = float(match.group(1))
    if match.group(2) is None:
        return value, None
    mantissa = match.group(1).split("e")[0].split("E")[0]
    decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
    return value, float(match.group(2)) * 10.0 ** -decimals


def _cells(row):
    return [re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", cell)).replace("\xa0", " ")).strip()
            for cell in re.findall(r"<td[^>]*>(.*?)</td>", row, flags=re.S)]


def read_ciaaw(standard_html, abridged_html):
    """`{Z: dict}` of symbol, name, abridged weight and uncertainty, and the full standard weight or its interval; an
    element CIAAW gives no standard weight has only its symbol and name."""
    elements = {}
    for text, kind in ((standard_html, "standard"), (abridged_html, "abridged")):
        for row in re.findall(r"<tr>(.*?)</tr>", text, flags=re.S):
            cells = _cells(row)
            if len(cells) < 4 or not cells[0].isdigit():
                continue
            z, symbol, name, weight = int(cells[0]), cells[1], cells[2], cells[3].replace(" ", "")
            entry = elements.setdefault(z, {"z": z, "symbol": symbol, "name": name.capitalize()})
            interval = re.fullmatch(r"\[([0-9.]+),([0-9.]+)\]", weight)
            if weight == "—":
                continue        # an element with no standard weight (no isotope is stable enough, or common enough)
            elif kind == "abridged":
                value, uncertainty = re.fullmatch(r"([0-9.]+)(?:±([0-9.]+))?", weight).groups()
                entry["atomic_weight"] = float(value)
                entry["atomic_weight_uncertainty"] = None if uncertainty is None else float(uncertainty)
            elif interval:
                entry["standard_interval"] = [float(interval.group(1)), float(interval.group(2))]
            else:
                value, uncertainty = _number(weight)
                entry["standard_weight"] = value
                entry["standard_weight_uncertainty"] = uncertainty
    return elements


def read_nist(text):
    """`{(Z, A): dict}` of symbol, relative atomic mass and uncertainty, and the isotopic composition (atom fraction) and
    its uncertainty, from NIST's ASCII listing."""
    isotopes = {}
    block = {}
    for line in text.splitlines():
        if line.startswith("Atomic Number"):
            block = {"z": int(line.split("=")[1])}
        elif "=" in line and block:
            key, _sep, value = line.partition("=")
            block[key.strip()] = value.strip()
            if key.strip() == "Notes":
                a = int(block["Mass Number"])
                mass, mass_unc = _number(block.get("Relative Atomic Mass", ""))
                fraction, fraction_unc = _number(block.get("Isotopic Composition", ""))
                isotopes[(block["z"], a)] = {"z": block["z"], "a": a, "symbol": block["Atomic Symbol"], "mass_u": mass,
                                             "mass_uncertainty_u": mass_unc, "abundance": fraction,
                                             "abundance_uncertainty": fraction_unc,
                                             "mass_estimated": "#" in block.get("Relative Atomic Mass", "")}
                block = {}
    return isotopes


def read_nubase(text):
    """`{(Z, A): dict}` of the ground states: mass excess (keV), half-life (s) and uncertainty, `stable`, the decay modes
    `[[mode, percent]]` and the isotopic abundance (percent). Values NUBASE marks as extrapolated (#) are left out."""
    states = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        a = int(line[0:3])
        zzzi = line[4:8]
        if zzzi[3] != "0":
            continue
        z = int(zzzi[:3])
        excess = line[18:31].strip()
        half = line[69:78].strip()
        unit = line[78:80].strip()
        unc = line[81:88].strip()
        record = {"z": z, "a": a, "mass_excess_kev": None if not excess or "#" in excess else float(excess)}
        if half == "stbl":
            record["stable"] = True
            record["half_life_s"] = None
        elif half and re.fullmatch(r"[0-9.]+", half) and unit in SECONDS:
            record["stable"] = False
            record["half_life_s"] = float(half) * SECONDS[unit]
            try:
                record["half_life_uncertainty_s"] = float(unc) * SECONDS[unit]
            except ValueError:
                record["half_life_uncertainty_s"] = None
        else:
            record["stable"] = False
            record["half_life_s"] = None
        record["modes"], record["nubase_abundance_percent"] = _modes(line[119:].strip())
        states[(z, a)] = record
    return states


def _modes(field):
    """`([[mode, percent]], isotopic abundance)` from NUBASE2020's decay-mode field. A branch given as about 100 percent
    (`A~100`) takes what the other branches leave and gets a third item `"~"`."""
    modes, abundance = [], None
    for token in field.split(";"):
        token = token.strip()
        match = re.fullmatch(r"([A-Za-z0-9+\-]+?)\s*([=~])\s*([0-9.]+(?:[eE][+-]?\d+)?)(?:\s+\d+)?", token)
        if not match:
            continue
        name, sign, value = match.groups()
        if name == "IS":
            abundance = float(value)
        elif sign == "~" and float(value) == 100.0:
            modes.append([name, 100.0, "~"])
        else:
            modes.append([name, float(value)])
    others = math.fsum(mode[1] for mode in modes if len(mode) == 2)
    for mode in modes:
        if len(mode) == 3:
            mode[1] = round(100.0 - others, 9)
    return modes, abundance


def _name(z, a, elements):
    return f"{elements[z]['symbol']}-{a}"


def _daughters(state):
    """`[(mode, percent, (Z, A))]` for the alpha, beta and electron-capture branches of a NUBASE state."""
    z, a = state["z"], state["a"]
    found = {}
    for mode, percent, *_flag in state["modes"]:
        if mode == "A":
            target, label = (z - 2, a - 4), "alpha"
        elif mode == "B-":
            target, label = (z + 1, a), "beta-"
        elif mode in ("B+", "EC", "EC+B+"):
            target, label = (z - 1, a), "ec"
        else:
            continue
        if percent >= MIN_BRANCH_PERCENT:
            found[(label, target)] = found.get((label, target), 0.0) + percent
    return [(label, percent, target) for (label, target), percent in sorted(found.items())]


def _q_kev(parent, daughter, label, states):
    me_parent, me_daughter = states[parent]["mass_excess_kev"], states.get(daughter, {}).get("mass_excess_kev")
    if me_parent is None or me_daughter is None:
        return None
    if label == "alpha":
        return round(me_parent - me_daughter - HE4_MASS_EXCESS_KEV, 3)
    return round(me_parent - me_daughter, 3)


def build(source_dir=None):
    global HE4_MASS_EXCESS_KEV
    ciaaw = read_ciaaw(_fetch("ciaaw-standard.html", source_dir), _fetch("ciaaw-abridged.html", source_dir))
    nist_all = read_nist(_fetch("nist-compositions.txt", source_dir))
    nist_natural = read_nist(_fetch("nist-natural.txt", source_dir))
    nubase = read_nubase(_fetch("nubase2020.txt", source_dir))
    HE4_MASS_EXCESS_KEV = nubase[(2, 4)]["mass_excess_kev"]
    by_name = {f"{ciaaw[z]['symbol']}-{a}": (z, a) for z in ciaaw for a in range(1, 320)}

    # --- the chains ---
    series, wanted = {}, {}
    for parent, end in NATURAL_SERIES.items():
        parent_key = by_name[parent]
        steps, seen, queue = [], set(), [parent_key]
        while queue:
            key = queue.pop(0)
            if key in seen:
                continue
            seen.add(key)
            for label, percent, target in _daughters(nubase[key]):
                if target not in nubase:
                    continue
                steps.append({"from": _name(key[0], key[1], ciaaw), "to": _name(target[0], target[1], ciaaw),
                              "mode": label, "percent": percent, "q_kev": _q_kev(key, target, label, nubase)})
                queue.append(target)
        members = sorted({step["from"] for step in steps} | {step["to"] for step in steps},
                         key=lambda name: (-by_name[name][1], by_name[name][0]))
        series[parent] = {"parent": parent, "end": end, "members": members, "steps": steps}
        for name in members:
            wanted[by_name[name]] = True
    for name in RADIONUCLIDES:
        wanted[by_name[name]] = True

    # --- the isotopes: every natural isotope, and the radioactive ones the table is for ---
    keys = sorted(set(nist_natural) | set(wanted))
    isotopes = []
    for key in keys:
        z, a = key
        nist = nist_all.get(key) or nist_natural.get(key)
        state = nubase.get(key)
        if nist is None or state is None:
            raise SystemExit(f"{_name(z, a, ciaaw)} is missing from {'NIST' if nist is None else 'NUBASE2020'}")
        natural = nist_natural.get(key, {})
        abundance = natural.get("abundance")
        name = _name(z, a, ciaaw)
        row = {"z": z, "a": a, "name": name, "mass_u": nist["mass_u"], "mass_uncertainty_u": nist["mass_uncertainty_u"],
               "abundance": abundance, "abundance_uncertainty": natural.get("abundance_uncertainty"),
               "mass_estimated": nist["mass_estimated"],
               "stable": state["stable"], "half_life_s": state["half_life_s"],
               "half_life_uncertainty_s": state.get("half_life_uncertainty_s"),
               "modes": state["modes"]}
        if name in LUGARO_ERG:
            row["decay_energy_j"] = LUGARO_ERG[name] * ERG_J
        if name in DYE_2012:
            row["dye_2012"] = DYE_2012[name]
        isotopes.append(row)
    elements = []
    for z in sorted(ciaaw):
        entry = dict(ciaaw[z])
        entry["has_natural_isotopes"] = any(key[0] == z for key in nist_natural)
        if "atomic_weight" not in entry:
            lived = [(state["half_life_s"], key[1]) for key, state in nubase.items()
                     if key[0] == z and state["half_life_s"] is not None]
            entry["longest_lived_mass_number"] = max(lived)[1] if lived else None
        elements.append(entry)

    common = {"retrieved": RETRIEVED}
    element_file = {
        **common, "title": "The elements: atomic number, symbol, name and atomic weight (GEN.438)",
        "sources": [
            "CIAAW, Standard Atomic Weights (ciaaw.org/atomic-weights.htm): the 'Atomic Weights 2021' report (Prohaska et al., "
            "Pure Appl. Chem. 94, 573, 2022, [R3615]) with the 2024 revisions of Gd, Lu and Zr",
            "CIAAW, Abridged Standard Atomic Weights (ciaaw.org/abridged-atomic-weights.htm): the weight to five digits and "
            "its uncertainty, the one value an element with an interval standard weight has"],
        "notes": "atomic_weight is the abridged conventional value in g/mol and atomic_weight_uncertainty its stated "
                 "uncertainty; standard_weight (and standard_weight_uncertainty) is the full-precision value where the "
                 "element has one, standard_interval the interval where it has that instead; an element with no stable "
                 "isotope or no standard weight has neither and longest_lived_mass_number names its longest-lived isotope (NUBASE2020).",
        "elements": elements}
    isotope_file = {
        **common, "title": "The isotopes: masses, natural abundance, half-lives, decay modes and decay energies (GEN.438)",
        "sources": [
            "NIST, Atomic Weights and Isotopic Compositions with Relative Atomic Masses (physics.nist.gov/Comp): the "
            "relative atomic mass (u) of each isotope and its natural isotopic composition (IUPAC)",
            "NUBASE2020 (Kondev et al., Chin. Phys. C 45, 030001, 2021, [R2443]), "
            "www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt: ground-state half-lives (the year is the Julian year of "
            "365.25 days), decay modes and their branchings in percent",
            "Lugaro, Ott and Kereszturi, Prog. Part. Nucl. Phys. 102, 1 (2018), [R2768: Table 4]: "
            "decay_energy_j, the energy released per decay (the whole chain for a parent that has one), without the "
            "neutrinos",
            "Dye, Rev. Geophys. 50, RG3007 (2012), [R1203: Table 2]: dye_2012, the heat constants of "
            "238U, 235U, 232Th and 40K (kept as the cross-check: the newer evaluation wins where they differ)"],
        "notes": "abundance is an atom fraction of the element in normal material, null for an isotope that is not "
                 "primordial; modes are NUBASE2020's [mode, percent] pairs; masses are of the neutral atom.",
        "isotopes": isotopes}
    series_file = {
        **common, "title": "The decay chains of 238U, 235U and 232Th (GEN.438)",
        "sources": ["NUBASE2020 [R2443]: the decay modes and branchings, followed from the parent "
                    f"down every branch of at least {MIN_BRANCH_PERCENT}% to the stable lead isotope",
                    "q_kev is the energy of the step computed from the mass excesses of NUBASE2020 (the atomic mass "
                    "evaluation AME2020): for alpha decay ME(parent) - ME(daughter) - ME(4He), otherwise "
                    "ME(parent) - ME(daughter); for electron capture and beta-plus that is the electron-capture Q value"],
        "notes": "Steps are between ground states; an isomer fed by a step (234mPa) is not listed, the chain's totals are "
                 "the decay_energy_j of the parent in isotopes.json.",
        "series": list(series.values())}
    return {"elements.json": element_file, "isotopes.json": isotope_file, "decay_series.json": series_file}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--from", dest="source_dir", help="a directory holding the source tables")
    parser.add_argument("--check", action="store_true", help="fail if the committed files differ from a fresh build")
    args = parser.parse_args()
    files = build(args.source_dir)
    status = 0
    for name, content in files.items():
        text = json.dumps(_clean(content), indent=1, ensure_ascii=False) + "\n"
        path = os.path.join(DATA, name)
        if args.check:
            with open(path, encoding="utf-8") as handle:
                if handle.read() != text:
                    print(f"{name} differs from a fresh build", file=sys.stderr)
                    status = 1
        else:
            os.makedirs(DATA, exist_ok=True)
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
            print(f"wrote {path}")
    return status


if __name__ == "__main__":
    sys.exit(main())
