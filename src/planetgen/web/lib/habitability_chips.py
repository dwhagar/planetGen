# planetgen/web/lib/habitability_chips.py

"""
The habitability chips of a planet's or moon's page (UX.181, UX.154): the "Habitable", "Habitable class", "Habitable moon",
"Habitable-class moon" and "Inhabited" flags and the equipment chips a human needs on the body. They sat in the system page's rows until the
object pages took the rows' place (UX.173); the symbols beside them are `habitability_symbols.symbol_row`.

Each chip links to its entry on the habitability page when `help_url` is given. Every value is escaped here.
"""

from planetgen.physics.habitability_world import DOMAINS, TIER_NAMES
from planetgen.web.lib.fmt import esc
from planetgen.web.lib.habitability_symbols import CHIP_ANCHORS, PHI6_CHIP_ANCHOR, equipment_anchor, help_link


def flag_chip(text, help_url):
    """A "flag-yes" chip, linked to its entry on the habitability page when `help_url` is given."""
    return help_link(f'<span class="flag flag-yes">{text}</span>', help_url, CHIP_ANCHORS[text.lower()])


def equipment_chips(body, help_url=None):
    """GEN.89, GEN.655, UX.154: the equipment a human needs on this body as one chip per item, solid for a required item
    and outlined for a recommended one, each in the colour of the PHI-6 domain it comes from; a body that needs nothing
    shows one chip of its worst PHI-6 domain's colour; empty for a gas giant or a body not scored. With `help_url` each
    chip links to the habitability page's entry for its item."""
    levels = body.get("equipment")
    if levels is None:
        return ""
    chips = []
    for item in levels["items"]:
        tier = max(body.get(f"tier_{domain.lower()}") or 0 for domain in item.get("domains") or [item["domain"]])
        text = item["item"][:1].upper() + item["item"][1:]
        hidden = ' <span class="sr-only">(recommended)</span>' if item["severity"] == "recommended" else ""
        chip = f'<span class="flag flag-eq flag-eq-{item["severity"]} flag-tier flag-tier-{tier}">{esc(text)}{hidden}</span>'
        chips.append(help_link(chip, help_url, equipment_anchor(item["item"])))
    if not chips:
        worst = max(body.get(f"tier_{domain}") or 0 for domain in DOMAINS)
        return help_link(f'<span class="flag flag-tier flag-tier-{worst}">PHI-6 {TIER_NAMES[worst]}</span>', help_url,
                         PHI6_CHIP_ANCHOR)
    return f'<span class="eq-chips">{"".join(chips)}</span>'


def own_flag(body):
    """UX.224: the one flag a body's own state earns, the most specific of three (an inhabited world is habitable and of a habitable class, a habitable one is, by its rating, a world a
    class could carry): "Inhabited", else "Habitable", else "Habitable class", else `None`."""
    if body.get("inhabited"):
        return "Inhabited"
    if body.get("habitable"):
        return "Habitable"
    if body.get("habitable_class"):
        return "Habitable class"
    return None


def moon_flag(moons):
    """UX.224: the one flag a planet's moons earn, "Habitable moon" over "Habitable-class moon", else `None`."""
    if any(moon.get("habitable") for moon in moons):
        return "Habitable moon"
    if any(moon.get("habitable_class") for moon in moons):
        return "Habitable-class moon"
    return None


def chips_html(body, help_url=None):
    """All the chips of a planet or moon row, one `<span class="stat">` each, or `""` for none: the one flag the body itself earns (`own_flag`: "Inhabited" over "Habitable" over "Habitable class"; "Habitable"
    is a Blue or Green overall rating, GEN.643, and "Habitable class" a class that can carry life), the one its moons earn (`moon_flag`), then the equipment chips (UX.224)."""
    chips = []
    for flag in (own_flag(body), moon_flag(body.get("moons") or ())):
        if flag:
            chips.append(flag_chip(flag, help_url))
    chips.append(equipment_chips(body, help_url))
    return "".join(f'<span class="stat">{chip}</span>' for chip in chips if chip)
