# planetgen/generation/star_labels.py

"""
How a run's logs name kinds of star: the sector fill's summary line ("3 G-type,
1 K-type giant") and the bright-star scatter's per-layer lines (GEN.131) use
the same words, so a log reads the same whichever step wrote it.
"""

GIANT_YERKES = ("III", "II")
SUPERGIANT_YERKES = ("IB", "IAB", "IA", "IA+", "0")


def star_label(star_type, yerkes_class):
    """
    `(label, plural)` for one star in a summary (UX.34): white dwarfs, black
    holes and neutron stars by name, giants and supergiants under their
    letter, others by spectral letter alone.

    Args:
        star_type (str): The spectral type ("G2", "K5"), or `None`.
        yerkes_class (str): The luminosity class ("V", "III", "D"), or `None`.
    """
    if yerkes_class in ("D", "VII"):
        return "white dwarf", "white dwarfs"
    if yerkes_class == "BH":
        return "black hole", "black holes"
    if yerkes_class == "NS":
        return "neutron star", "neutron stars"
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
