# planetgen/web/lib/habitability_ref.py

"""
The text of the habitability levels page (`/classes/habitability`, UX.90):
the five equipment levels, the four PHI-4 domains with their Blue, Green,
Yellow and Red thresholds, the compact-host rule and the chips a planet row
shows. The equipment labels come from `planetgen.physics.habitability_world`;
the thresholds repeat `docs/design/habitability-index.md` section 2 and
`planetgen.physics.habitability`, where they are named with their sources.
"""

from planetgen.physics.habitability_world import EQUIPMENT_LABELS, LETHAL_DOSE_SV_YR

EQUIPMENT_MEANINGS = (
    "Ordinary weather clothing is enough: the air has the right pressure and oxygen, no harmful gas, "
    "a comfortable temperature and a low radiation dose.",
    "A breathing mask: oxygen is low or the pressure is a little off, but a mask makes the air safe.",
    "A mask with a scrubber: a gas such as carbon dioxide, carbon monoxide, hydrogen sulfide or sulfur dioxide "
    "is above its long-term limit and has to be filtered out.",
    "A sealed suit: the pressure is too low for a mask to work, or the temperature or the air's chemistry "
    "is outside what a person can stand.",
    "Full life support with radiation hardening: the radiation dose is beyond what any suit can shield, "
    "so a visitor needs a hardened habitat.",
)
"""tuple: What each equipment level means, in the order of `EQUIPMENT_LABELS`."""

DOMAIN_TIERS = (
    ("Blue", "The best band: a visitor notices nothing wrong."),
    ("Green", "Workable with ordinary care."),
    ("Yellow", "Needs equipment or limits the time spent outside."),
    ("Red", "Beyond what a person can survive without heavy protection."),
)
"""tuple: `(colour, meaning)` of a PHI-4 domain's score, best first."""

DOMAINS = (
    ("Pressure", [
        ("Air pressure", "50 to 250 kPa", "14.3 to 400 kPa", "6.3 to 10,000 kPa"),
        ("Oxygen the lungs get (the pressure minus 6.3 kPa of water vapour)", "8 to 50 kPa",
         "under 8 kPa (a mask works)", "50 to 160 kPa"),
    ]),
    ("Temperature", [
        ("Air temperature", "0 to 31 C", "-20 to 45 C", "-50 to 122 C"),
        ("Wet-bulb temperature (heat and humidity together)", "under 31 C",
         "under 31 C", "31 C or more (a body can no longer shed heat)"),
    ]),
    ("Chemistry", [
        ("Carbon dioxide", "under 0.93 kPa", "up to 2.0 kPa", "up to 5.0 kPa"),
        ("Carbon monoxide, hydrogen sulfide, sulfur dioxide", "under the long-term limit", "",
         "up to ten times the short-term limit"),
        ("Water acidity (pH)", "6 to 8.5", "5 to 9.5", "1 to 11.5"),
        ("Water activity (how freely water is available)", "0.90 or more", "0.75 or more", "0.605 or more"),
    ]),
    ("Radiation", [
        ("Dose per year", "under 0.05 Sv", "up to 0.1 Sv", "up to 10 Sv"),
    ]),
)
"""tuple: `(domain, [(variable, blue, green, yellow)])`: beyond the Yellow
column a variable is Red."""

SCORE_RULE = (
    "Each variable scores 1 in Blue and falls to 0.75 at the edge of Green, 0.4 at the edge of Yellow and 0 "
    "at the outer edge of Red. A domain takes the score of its worst variable: Blue at 1, Green from 0.75, "
    "Yellow from 0.4, Red below. The PHI-4 score on a "
    "planet's tooltip is the geometric mean of the four domains, so one domain at 0 makes the whole score 0."
)

EQUIPMENT_RULE = (
    "The equipment level comes from the worst domain: radiation in the Red band needs full life support; "
    "pressure below Green, or any domain in the Red band, needs a sealed suit; a gas above its long-term limit "
    "needs a mask with a scrubber; low oxygen or a pressure below Blue needs a mask; anything else is ideal."
)

COMPACT_HOST_RULE = (
    f"A planet of a pulsar, a neutron star or a black hole is rated as if its surface took "
    f"{LETHAL_DOSE_SV_YR:,.0f} Sv a year. That puts radiation in Red, so the planet needs full life support "
    "and its life scores are 0. Its note names the host."
)

CHIPS = (
    ("Habitable", "The planet's class can carry life (an Earth-like or similar world). It is not a promise "
                  "that anything lives there."),
    ("Habitable moon", "A moon of this planet is habitable in that sense."),
    ("Inhabited", "A species lives here."),
)
"""tuple: The chips a planet or moon row can show beside its equipment level."""


def page():
    """Everything the page shows."""
    return {
        "equipment": list(zip(EQUIPMENT_LABELS, EQUIPMENT_MEANINGS)),
        "domain_tiers": DOMAIN_TIERS,
        "domains": DOMAINS,
        "score_rule": SCORE_RULE,
        "equipment_rule": EQUIPMENT_RULE,
        "compact_host_rule": COMPACT_HOST_RULE,
        "chips": CHIPS,
    }
