# planetgen/generation/star_labels.py

"""
How a run's logs name kinds of star: the sector fill's summary line ("3 G-type,
1 K-type giant") and the bright-star scatter's per-shell lines (GEN.131) use
the same words, so a log reads the same whichever step wrote it.
"""

GIANT_YERKES = ("III", "II")
SUPERGIANT_YERKES = ("IB", "IAB", "IA", "IA+", "0")


MEMBER_KIND_LABELS = {
    "stripped-helium": ("stripped helium star", "stripped helium stars"),
    "thorne-zytkow": ("Thorne-Zytkow object", "Thorne-Zytkow objects"),
    "wolf-rayet": ("Wolf-Rayet star", "Wolf-Rayet stars"),
}
"""dict: The member kinds of the multi-star catalogue that have no Yerkes class of their own, named the same way."""

NAMED_CLASS_LABELS = {
    **MEMBER_KIND_LABELS,
    "D": ("white dwarf", "white dwarfs"), "VII": ("white dwarf", "white dwarfs"),
    "protostar": ("protostar", "protostars"),
    "W": MEMBER_KIND_LABELS["wolf-rayet"],
    "PMS": ("pre-main-sequence star", "pre-main-sequence stars"),
    "BH": ("black hole", "black holes"),
    "NS": ("neutron star", "neutron stars"),
}
"""dict: The classes and kinds named by what they are, not by their spectral letter: `(label, plural)` by Yerkes class or member kind."""


def star_label(star_type, yerkes_class):
    """
    `(label, plural)` for one star in a summary (UX.34): white dwarfs, black
    holes and neutron stars by name, giants and supergiants under their
    letter, others by spectral letter alone.

    Args:
        star_type (str): The spectral type ("G2", "K5"), or `None`.
        yerkes_class (str): The luminosity class ("V", "III", "D"), or `None`; a key of
            `MEMBER_KIND_LABELS` names that kind.
    """
    named = NAMED_CLASS_LABELS.get(yerkes_class)
    if named is not None:
        return named
    letter = (star_type or "?")[0]
    if yerkes_class in GIANT_YERKES:
        return f"{letter}-type giant", f"{letter}-type giants"
    if yerkes_class in SUPERGIANT_YERKES:
        return f"{letter}-type supergiant", f"{letter}-type supergiants"
    return f"{letter}-type", f"{letter}-type"


def describe_types(counts):
    """`"3 G-type, 1 K-type giant"` from a `Counter` of `star_label` pairs, in label order."""
    return ", ".join(f"{count} {label if count == 1 else plural}"
                     for (label, plural), count in sorted(counts.items()))
