# planetgen/util/chem.py

"""
Chemical formulas (UX.149)
==========================

The one way a chemical formula is written out: `html` gives `H<sub>2</sub>O` for a page (and for Markdown and
wikitext, which take the same tags), `text` its plain-text twin `H₂O` for a terminal. Counts are subscripts, ion
charges superscripts after the formula and isotope mass numbers superscripts before their symbol:

- `H2O`, `C2H6`, `(CH3)2S`, `Ca5(PO4)3F`, `(Mg0.88,Fe0.12)SiO3`, `CxHy`: counts, brackets, decimal and variable counts.
- `CaSO4·2H2O` (or `CaSO4.2H2O`): a hydrate; the coefficient after the dot stays on the line.
- `3He`, `14C`, `13CO`: an isotope's mass number before the first symbol.
- `Na+`, `Cl-`, `NH4+`, `e-`: a charge of one; `Fe2+`, `O2+` (one element: the number is the charge, so `O2+` is
  doubly ionised oxygen) and `SO4^2-`, `O2^+` (a `^` names the charge of anything else).

These take one formula, never free text: a caller that prints a text naming chemicals (a nebula's contents, a limit's
meaning) splits it with `split`, which marks only what the registry (`known`) holds or a charged element, so a word that
happens to look like a formula elsewhere is never guessed at. The Jinja filters `formula` and `chem` (web/__init__.py)
and the Document's `chem` inline item (`planetgen.describe.model.chem`) are the ways a page or a description reaches
these.
"""

import functools
import html as _html
import json
import os
import re

_SUBSCRIPT = str.maketrans("0123456789+-=()x", "₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₓ")
_SUPERSCRIPT = str.maketrans("0123456789+-=()", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾")
MINUS = "−"
"""str: The minus sign a charge is written with on a page."""
HYDRATE_DOT = "\u00b7"
"""str: The middle dot between a salt and its water of crystallisation."""

_SYMBOL = re.compile(r"[A-Z][a-z]?")
_COUNT = re.compile(r"\d+(?:\.\d+)?|[xyn](?![a-z])")
_EXPRESSION_COUNT = re.compile(r"\((?:\d|[xyn])[0-9xyn.+-]*\)")
HEAVY_HYDROGEN = frozenset({"D", "T"})
"""frozenset: Symbols a formula may hold beside the elements: deuterium and tritium (`HDO`)."""
_HYDRATE = re.compile("[.\u00b7]" r"(\d*)(?=[A-Z(\[])")
_CHARGE_TAIL = re.compile(r"\^(\d*[+-])$|(?<=[A-Za-z\])])(\d*)([+-])$")
_OPEN, _CLOSE = "([", ")]"


class FormulaError(ValueError):
    """A string `parse` cannot read as a formula."""


def parse(formula):
    """
    `formula` as runs: `[(kind, text)]`, `kind` one of `"text"`, `"sub"` and `"sup"`, neighbouring runs of one kind joined.
    A charge's sign is `+` or `-` here; `html` writes the minus as `MINUS`.

    Raises:
        FormulaError: Not a formula (a space, an unknown character, an unbalanced bracket, an empty string).
    """
    if not isinstance(formula, str) or not formula or formula != formula.strip() or " " in formula:
        raise FormulaError(f"{formula!r} is not a formula")
    if formula == "e-":
        return [("text", "e"), ("sup", "-")]
    body, charge = _split_charge(formula)
    runs = []
    _read(body, runs, formula)
    if charge:
        runs.append(("sup", charge))
    return _joined(runs)


def _split_charge(formula):
    """`(body, charge)`: `Fe2+` is `("Fe", "2+")`, `NH4+` `("NH4", "+")`, `SO4^2-` `("SO4", "2-")`."""
    match = _CHARGE_TAIL.search(formula)
    if not match:
        return formula, ""
    if match.group(1) is not None:
        return formula[:match.start()], match.group(1)
    digits, sign = match.group(2), match.group(3)
    head = formula[:match.start()]
    one_element = _SYMBOL.fullmatch(head.lstrip("0123456789")) is not None
    if digits and one_element:
        return head, digits + sign           # Fe2+, O2+: one element, so the number is its charge
    return head + digits, sign               # NH4+: the number is a count and the charge is one


def _read(body, runs, whole):
    depth = 0
    i = 0
    start = True
    after_unit = False
    while i < len(body):
        char = body[i]
        if start and char.isdigit():
            match = re.match(r"\d+(?=[A-Z])", body[i:])
            if not match:
                raise FormulaError(f"{whole!r}: a number with no symbol after it")
            runs.append(("sup", match.group()))  # an isotope's mass number
            i += match.end()
            start = False
            continue
        hydrate = _HYDRATE.match(body, i) if after_unit and depth == 0 else None
        if hydrate and (char == HYDRATE_DOT or not _decimal_count_here(body, i)):
            runs.append(("text", HYDRATE_DOT + hydrate.group(1)))
            i = hydrate.end()
            start = True
            after_unit = False
            continue
        symbol = _symbol_at(body, i)
        if symbol:
            runs.append(("text", symbol))
            i += len(symbol)
        elif char in _OPEN:
            depth += 1
            runs.append(("text", char))
            i += 1
            start = True
            after_unit = False
            continue
        elif char in _CLOSE:
            depth -= 1
            if depth < 0:
                raise FormulaError(f"{whole!r}: a bracket closes that never opened")
            runs.append(("text", char))
            i += 1
        elif char == "," and depth > 0:
            runs.append(("text", char))
            i += 1
            start = True
            after_unit = False
            continue
        else:
            raise FormulaError(f"{whole!r}: {char!r} at {i} is not part of a formula")
        start = False
        after_unit = True
        count = _COUNT.match(body, i)
        expression = _EXPRESSION_COUNT.match(body, i)
        if count:
            number = count.group()
            if "." in number and not _decimal_count_here(body, i + number.index(".")):
                number = number[:number.index(".")]
            runs.append(("sub", number))
            i += len(number)
        elif expression:
            runs.append(("sub", expression.group()))   # Fe(1-x)O: a count that is an expression
            i = expression.end()
    if depth:
        raise FormulaError(f"{whole!r}: a bracket is never closed")


def _symbol_at(body, i):
    """The element symbol at `i` (two letters when those are an element, else one), or `None`."""
    symbols = _element_symbols() | HEAVY_HYDROGEN
    for size in (2, 1):
        candidate = body[i:i + size]
        if len(candidate) == size and _SYMBOL.fullmatch(candidate) and candidate in symbols:
            return candidate
    return None


def _decimal_count_here(body, dot):
    """Whether the `.` at `dot`, between digits, is a decimal point (`Mg0.88,`) rather than a hydrate's dot
    (`CaSO4.2H2O`): it is when the digits after it are not followed by a symbol or a bracket that opens."""
    if dot >= len(body) or body[dot] != "." or dot == 0 or not body[dot - 1].isdigit():
        return False
    match = re.match(r"\.(\d+)", body[dot:])
    if not match:
        return False
    after = dot + match.end()
    return after >= len(body) or not (body[after].isupper() or body[after] in "([")


def _joined(runs):
    out = []
    for kind, text in runs:
        if out and out[-1][0] == kind:
            out[-1] = (kind, out[-1][1] + text)
        else:
            out.append((kind, text))
    return out


def html(formula):
    """`formula` as escaped HTML with `<sub>` and `<sup>` (a str; the web's `formula` filter makes it Markup). The same
    tags are what Markdown and wikitext are given."""
    out = []
    for kind, text in parse(formula):
        text = _html.escape(text)
        if kind == "text":
            out.append(text)
        else:
            out.append(f"<{kind}>{text.replace('-', MINUS) if kind == 'sup' else text}</{kind}>")
    return "".join(out)


def text(formula):
    """`formula` as plain text with Unicode subscripts and superscripts (`H₂O`, `Fe²⁺`, `³He`): `html`'s twin for a
    terminal. A subscript Unicode has no form for (a decimal point, `y`, `n`) stays on the line."""
    out = []
    for kind, part in parse(formula):
        if kind == "sub":
            out.append(part.translate(_SUBSCRIPT))
        elif kind == "sup":
            out.append(part.translate(_SUPERSCRIPT))
        else:
            out.append(part)
    return "".join(out)


@functools.cache
def known():
    """
    frozenset: The formulas `split` marks: every gas (`physics/data/gases.json`) and every material's formula
    (`physics/data/materials.json`) that has a count, a charge, an isotope or two symbols to write out. A charged
    element (`Fe2+`, `N+`) and the electron (`e-`) are marked too, by `is_formula`.
    """
    data = os.path.join(os.path.dirname(os.path.dirname(__file__)), "physics", "data")
    with open(os.path.join(data, "gases.json"), encoding="utf-8") as handle:
        formulas = {row["formula"] for row in json.load(handle)["gases"]}
    with open(os.path.join(data, "materials.json"), encoding="utf-8") as handle:
        formulas |= {row["formula"] for row in json.load(handle)["materials"] if row.get("formula")}
    return frozenset(f for f in formulas if parse(f) != [("text", f)])


@functools.cache
def _element_symbols():
    from planetgen.physics.elements import ELEMENTS
    return frozenset(element.symbol for element in ELEMENTS)


def is_formula(token):
    """Whether `split` marks `token`: a `known` formula with or without a charge, a charged element, or `e-`."""
    if token == "e-":
        return True
    try:
        body, charge = _split_charge(token)
    except FormulaError:
        return False
    if body in known():
        return True
    return bool(charge) and body in _element_symbols()


def split(words):
    """
    A text that names chemicals, as a list of plain strings and `Formula` items: the tokens `is_formula` accepts (see
    the module docstring: the registry, never a guess). `"H+, e-, O2+ and S+"` is
    `[Formula("H+"), ", ", Formula("e-"), ", ", Formula("O2+"), " and ", Formula("S+")]`.
    """
    out = []
    last = 0
    for match in re.finditer(r"[^\s,;:]+", words):
        token = match.group().rstrip(".)")
        token = token[1:] if token.startswith("(") and token.count("(") > token.count(")") else token
        begin = match.start() + match.group().index(token) if token else match.start()
        if token and is_formula(token):
            if begin > last:
                out.append(words[last:begin])
            out.append(Formula(token))
            last = begin + len(token)
    if last < len(words):
        out.append(words[last:])
    return out


class Formula(str):
    """A formula `split` found in a text: a str (so it reads as itself) that the writers below mark up."""

    __slots__ = ()


def mark_html(words):
    """A text naming chemicals as escaped HTML, each formula `split` finds written by `html`."""
    return "".join(html(part) if isinstance(part, Formula) else _html.escape(part) for part in split(words))


def mark_text(words):
    """`mark_html`'s plain-text twin: each formula `split` finds written by `text`."""
    return "".join(text(part) if isinstance(part, Formula) else part for part in split(words))
