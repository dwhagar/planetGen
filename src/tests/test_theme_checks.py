# tests/test_theme_checks.py

"""UX.82: the theme checks carried over from PR #800: the icon files draw in `currentColor`, the Shoelace hover colour clears 4.5:1 in
both themes, the placeholder does too, and `--bg-subtle` is an opaque colour in every theme."""

import pathlib
import re

STATIC = pathlib.Path(__file__).resolve().parents[1] / "html" / "static"
CSS = (STATIC / "style.css").read_text(encoding="utf-8")
SHOELACE_THEME = (STATIC / "shoelace-theme.css").read_text(encoding="utf-8")
LIGHT_ROOT = re.search(r"^:root \{\n(.*?)^\}", CSS, re.S | re.M).group(1)
DARK_ROOT = re.search(r':root\[data-theme="dark"\] \{\n(.*?)^\}', CSS, re.S | re.M).group(1)
BG_SUBTLE_VALUES = re.findall(r"^\s*--bg-subtle: (#[0-9a-fA-F]+);", CSS, re.M)


def _token(block, name):
    return re.search(rf"^\s*--{name}: (#[0-9a-fA-F]{{6}});", block, re.M).group(1)


def _rgb(hex_colour):
    return tuple(int(hex_colour[i : i + 2], 16) / 255.0 for i in (1, 3, 5))


def _mix(first, share, second):
    """`color-mix(in srgb, first share%, second)` on 0..1 channels."""
    return tuple(a * share + b * (1.0 - share) for a, b in zip(first, second))


def _luminance(rgb):
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def _contrast(first, second):
    high, low = sorted((_luminance(first), _luminance(second)), reverse=True)
    return (high + 0.05) / (low + 0.05)


def _primary_500(block):
    """The hover fill of a primary button, from the formula in shoelace-theme.css, in the tokens of `block`."""
    match = re.search(r"--sl-color-primary-500: color-mix\(in srgb, var\(--accent\) (\d+)%, var\(--(\w[\w-]*)\)\);", SHOELACE_THEME)
    share, other = int(match.group(1)) / 100.0, match.group(2)
    return _mix(_rgb(_token(block, "accent")), share, _rgb(_token(block, other)))


def test_every_icon_file_draws_in_current_color():
    files = [STATIC / "icons.svg"] + sorted((STATIC / "vendor" / "shoelace" / "assets").rglob("*.svg"))
    assert len(files) > 10
    for path in files:
        text = path.read_text(encoding="utf-8")
        assert "currentColor" in text, f"{path.name} has no currentColor"
        assert not re.search(r'(?:fill|stroke)="#(?:[0-9a-fA-F]{3}){1,2}"', text), f"{path.name} hard-codes a colour"


def test_the_primary_button_hover_clears_four_and_a_half_to_one_in_both_themes():
    for name, block in (("light", LIGHT_ROOT), ("dark", DARK_ROOT)):
        fill = _primary_500(block)
        text = _rgb(_token(block, "bg-elevated"))
        assert _contrast(fill, text) >= 4.5, f"{name} theme primary hover {_contrast(fill, text):.2f}:1"


def test_the_placeholder_clears_four_and_a_half_to_one_in_both_themes():
    assert "--sl-input-placeholder-color: var(--text-muted);" in SHOELACE_THEME
    for name, block in (("light", LIGHT_ROOT), ("dark", DARK_ROOT)):
        assert _contrast(_rgb(_token(block, "text-muted")), _rgb(_token(block, "bg-elevated"))) >= 4.5, name


def test_bg_subtle_is_an_opaque_colour_in_every_theme():
    # The light root, the dark theme by preference, the dark theme by choice and the map viewports.
    assert len(BG_SUBTLE_VALUES) == 4
    assert all(len(value) == 7 for value in BG_SUBTLE_VALUES), BG_SUBTLE_VALUES
