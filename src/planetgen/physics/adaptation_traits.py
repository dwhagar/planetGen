# planetgen/physics/adaptation_traits.py

"""
Adaptation traits (GEN.876, docs/design/habitability/adaptation-traits.md): how the life of a hard world is tougher than the life of an easy one.

A trait widens one named limit of the habitability index that the body's own stress presses on, and nothing else. Three stresses have a limit in
the index to widen, so three traits are built:

* `radiation`: the absorbed dose where the life is (`habitability.l_rad`); a microbe's threshold and an animal's are separate numbers
  (`habitability.MICROBE_DOSE_THRESH_GY`, `habitability.DOSE_THRESH_GY`), and a trait multiplies each side's own;
* `sulfide`: the hydrogen sulfide scale of the toxins of the air that animals meet (`habitability.H2S_TOX_KPA` in `habitability.toxin_factor`);
* `carbon_dioxide`: its carbon dioxide scale (`habitability.CO2_TOX_KPA`).

The other stresses of the design have no limit in the index to move, so no trait is drawn for them (registry rows LF70 and LF71): brine (the
water-activity floor is one number for all life and the class table reads it), heat and cold (the class windows are the life's own bands), ultraviolet,
acid and alkali (no class term reads them), low oxygen (the life builds the oxygen, and each class has its own floors, GEN.877), and the carbon monoxide
and sulfur dioxide toxins (only animals are scored against them, and no animal limit exists for a microbe to widen).

A trait has a grade: 1 ordinary (none), 2 tolerant, 3 extreme. It is stored as text (`life_traits`, `radiation:2,sulfide:3`, `None` for none) and shown as
"adapted to: radiation (tolerant), sulfide (extreme)". A trait never makes a world habitable to people: the human scores never read one.

This module holds the grades, their factors and the text; the draw, which reads the body's stress and the limits it moves, is `physics/adaptation_draw.py`.
Every number here is a stand-in with a registry row (`docs/design/research/gap-filling-plan.md` section 6.18, LF66 to LF69).
"""

STRESSES = ("radiation", "sulfide", "carbon_dioxide")
"""tuple: The stresses that have a trait, in the order they are drawn and listed."""

GRADES = (1, 2, 3)
GRADE_WORDS = {1: "ordinary", 2: "tolerant", 3: "extreme"}
"""dict: The words for a grade, as the object page shows it (grade 1 is never shown)."""

STRESS_WORDS = {"radiation": "radiation", "sulfide": "sulfide", "carbon_dioxide": "carbon dioxide"}
"""dict: Each stress as a person reads it."""

MICROBE_FACTOR = {"radiation": {1: 1.0, 2: 10.0, 3: 100.0}}
"""dict: By stress, by grade, how far a trait moves the limit a microbe meets. Only radiation has one: no microbe limit on sulfide or carbon dioxide exists in
the index (the toxin terms are animals'). Tier T9, tweak level L3 (registry LF66)."""

ANIMAL_FACTOR = {
    "radiation": {1: 1.0, 2: 10.0, 3: 10.0},
    "sulfide": {1: 1.0, 2: 10.0, 3: 100.0},
    "carbon_dioxide": {1: 1.0, 2: 2.0, 3: 5.0},
}
"""dict: By stress, by grade, how far a trait moves the limit an animal meets: capped below a microbe's, because the best animals are far less tolerant
than the best microbes (radiation tops out at grade 2, carbon dioxide at 5 times). Tier T7 for the shape, T9 for the numbers, tweak level L3 (registry LF67)."""

P_CAP = 0.8
"""float: The most likely a trait is on the harshest world (it is `P_CAP` times the stress score). Tier T9, tweak level L3 (registry LF68)."""
GRADE3_SHARE = 0.5
"""float: A trait is extreme (grade 3) if its second draw is under `GRADE3_SHARE` times the square of the stress score, tolerant (grade 2) otherwise.
Tier T9, tweak level L3 (registry LF68)."""


def factor(traits, stress, animal):
    """How far the traits (`{stress: grade}`, or `None`) move the limit of `stress` for an animal (`animal=True`) or a microbe; 1 for no trait."""
    grade = (traits or {}).get(stress, 1)
    table = ANIMAL_FACTOR if animal else MICROBE_FACTOR
    return table.get(stress, {}).get(grade, 1.0)


def parse(text):
    """`{stress: grade}` of a stored `life_traits` text (`None` or empty text is none). An unknown stress or a grade that is not 2 or 3 is an error."""
    traits = {}
    for part in (text or "").split(","):
        if not part.strip():
            continue
        stress, _, grade = part.strip().partition(":")
        if stress not in STRESSES or grade not in ("2", "3") or stress in traits:
            raise ValueError(f"not an adaptation trait: {part!r}")
        traits[stress] = int(grade)
    return traits


def text_of(traits):
    """The `life_traits` text of `{stress: grade}`, in the order of `STRESSES`, grades 2 and 3 only; `None` for none."""
    parts = [f"{stress}:{traits[stress]}" for stress in STRESSES if traits.get(stress, 1) > 1]
    return ",".join(parts) or None


def sentence(text):
    """The stored `life_traits` as a person reads it ("radiation (tolerant), sulfide (extreme)"), or `None` for none."""
    traits = parse(text)
    return ", ".join(f"{STRESS_WORDS[stress]} ({GRADE_WORDS[traits[stress]]})" for stress in STRESSES if stress in traits) or None
