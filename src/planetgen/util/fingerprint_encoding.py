# planetgen/util/fingerprint_encoding.py

"""
The canonical encoding every fingerprint hashes (GEN.304)
=========================================================

One encoder for the digests of generated values, so a galaxy's fingerprint (`db/fingerprint.py`, GEN.58), the
generator battery's digest (OPS.28) and the golden tests (TEST.77) hash a value one way, and a 1-ulp difference in
`exp` or `sin` between two machines does not show (docs/design/reproducible-galaxies.md, generation-determinism.md
section 4).

A value is a tree of `dict` (string keys), `list` or `tuple`, `str`, `bool`, `None`, `int`, `float`, `Decimal` and
`bytes`; `encode` writes it as text:

* a dict's keys are sorted and a list keeps its order (the generation order), so the caller sorts what has no order;
* an integer is exact: decimal up to 2**53, `0x` hexadecimal above it (and `-0x` below its negative);
* a float is `format(x, ".8e")`, nine significant digits, with -0.0 as 0.0 and NaN and the infinities as `nan`,
  `inf` and `-inf`; a `Decimal` is a float to the encoder;
* a string is JSON-quoted with every non-ASCII character escaped, so no Unicode normalisation happens and the text
  holds on any platform;
* bytes are `b"<hex>"`; there are no ids, timestamps or other fields of the encoder's own.

A digest is the SHA-256 of `FORMAT_TAG`, the kind and the text, so a format change is a new tag (`pgfp2`) and never a
silent drift. Three kinds are kept: `full` (the text above), `structural` (every float written as `~`, so only the
tree, the integers and the strings count) and `floats` (the float tokens alone, in tree order). A failing comparison
reads the structural digest first: if it differs the shape of the data changed, if only the float digest differs a
number did.

The rounding is of the text, not of the value: a float within one ulp of a rounding boundary of the ninth digit can
still change its token (a few in a hundred million values), so a golden test that sees a lone float flip across
machines reads the value, not the digest.
"""

import dataclasses
import decimal
import hashlib

FORMAT_TAG = "pgfp1"
"""str: The version of the encoding, hashed into every digest."""

FLOAT_FORMAT = ".8e"
"""str: A float's text, one digit and eight decimals of the mantissa: nine significant digits."""

EXACT_INT_LIMIT = 2 ** 53
"""int: The largest magnitude an integer is written in decimal; above it, in hexadecimal."""

KINDS = ("full", "structural", "floats")

_FLOAT_TYPES = (float, decimal.Decimal)

_ESCAPES = {'"': '\\"', "\\": "\\\\", "\n": "\\n", "\r": "\\r", "\t": "\\t"}


def _string(text):
    out = ['"']
    for char in text:
        code = ord(char)
        if char in _ESCAPES:
            out.append(_ESCAPES[char])
        elif code < 0x20 or code > 0x7e:
            out.append(f"\\u{code:04x}" if code < 0x10000 else f"\\U{code:08x}")
        else:
            out.append(char)
    out.append('"')
    return "".join(out)


def float_token(value):
    """A float (or `Decimal`) as its canonical text: nine significant digits, one zero, NaN and infinities named."""
    if isinstance(value, decimal.Decimal):
        if value.is_nan():
            return "nan"
        if value.is_infinite():
            return "inf" if value > 0 else "-inf"
        value = float(value)
    if value != value:
        return "nan"
    if value in (float("inf"), float("-inf")):
        return "inf" if value > 0 else "-inf"
    if value == 0:
        return format(0.0, FLOAT_FORMAT)
    return format(value, FLOAT_FORMAT)


def int_token(value):
    """An integer as its canonical text: exact decimal up to 2**53, hexadecimal above."""
    if abs(value) <= EXACT_INT_LIMIT:
        return str(value)
    return f"-0x{-value:x}" if value < 0 else f"0x{value:x}"


def _walk(value, floats, out):
    if value is None:
        out.append("null")
    elif value is True:
        out.append("true")
    elif value is False:
        out.append("false")
    elif isinstance(value, int):
        out.append(int_token(value))
    elif isinstance(value, _FLOAT_TYPES):
        out.append(float_token(value) if floats else "~")
    elif isinstance(value, str):
        out.append(_string(value))
    elif isinstance(value, (bytes, bytearray, memoryview)):
        out.append(f'b"{bytes(value).hex()}"')
    elif isinstance(value, dict):
        out.append("{")
        for number, key in enumerate(sorted(value)):
            if not isinstance(key, str):
                raise TypeError(f"a fingerprint key must be a string, not {type(key).__name__}")
            if number:
                out.append(",")
            out.append(_string(key))
            out.append(":")
            _walk(value[key], floats, out)
        out.append("}")
    elif isinstance(value, (list, tuple)):
        out.append("[")
        for number, item in enumerate(value):
            if number:
                out.append(",")
            _walk(item, floats, out)
        out.append("]")
    else:
        raise TypeError(f"cannot fingerprint a {type(value).__name__}")


def encode(value, floats=True):
    """
    The canonical text of `value`.

    Args:
        value: A tree of the types in the module docstring.
        floats (bool): `False` writes every float as `~` (the structural text).

    Returns:
        str: The text.
    """
    out = []
    _walk(value, floats, out)
    return "".join(out)


def float_tokens(value):
    """The canonical float tokens of `value`, in tree order (dict keys sorted)."""
    tokens = []

    def visit(item):
        if isinstance(item, _FLOAT_TYPES) and not isinstance(item, bool):
            tokens.append(float_token(item))
        elif isinstance(item, dict):
            for key in sorted(item):
                visit(item[key])
        elif isinstance(item, (list, tuple)):
            for entry in item:
                visit(entry)

    visit(value)
    return tokens


def text_digest(kind, text):
    """The hex SHA-256 of `FORMAT_TAG`, `kind` and `text`."""
    if kind not in KINDS:
        raise ValueError(f"unknown digest kind {kind!r}")
    return hashlib.sha256(f"{FORMAT_TAG}:{kind}\n{text}".encode("ascii")).hexdigest()


@dataclasses.dataclass(frozen=True)
class Digests:
    """
    The three digests of one value.

    Attributes:
        full (str): Over the whole text, floats at nine significant digits.
        structural (str): Over the text with every float written as `~`.
        floats (str): Over the float tokens alone, in tree order.
    """

    full: str
    structural: str
    floats: str


def digests(value):
    """The `Digests` of `value`."""
    return Digests(text_digest("full", encode(value)), text_digest("structural", encode(value, floats=False)),
                   text_digest("floats", "\n".join(float_tokens(value))))


def digest(value):
    """The full digest of `value`."""
    return text_digest("full", encode(value))


def sorted_rows(rows):
    """
    `rows` (values with no order of their own) in the order every digest of them uses: by their structural text, then
    by their full text, so a float's last digits never decide an order that the structure settles.
    """
    return sorted(rows, key=lambda row: (encode(row, floats=False), encode(row)))


def row_set_digests(rows):
    """The `Digests` of an unordered set of rows, taken in the order of `sorted_rows`."""
    return digests(sorted_rows(rows))
