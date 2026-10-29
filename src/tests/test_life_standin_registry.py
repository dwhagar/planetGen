"""
The stand-in registry of the life constants (GEN.870, docs/design/research/gap-filling-plan.md section 6.18): every row names a
constant that exists, and every constant of the life code either has a row or is on the list of those that need none.
"""

import importlib
import pathlib
import re

import pytest

from planetgen import tuning
from planetgen.physics import habitability, habitability_world

PLAN = pathlib.Path(__file__).resolve().parents[2] / "docs" / "design" / "research" / "gap-filling-plan.md"
MODULES = {"tuning": tuning, "habitability": habitability, "habitability_world": habitability_world}
NAME = re.compile(r"\b(tuning|habitability_world|habitability)\.([A-Z][A-Z0-9_]*)\b")

TUNING_PREFIXES = ("LIFE_", "STAR_EVOLUTION", "EVOLUTIONARY_", "HABITABILITY_LIFE_")

NEEDS_NO_ROW = {
    "tuning.LIFE_METABOLISMS": "words for the photosynthesis milestone by metabolism: text, no number",
    "tuning.EVOLUTIONARY_TEXT": "the sentence of each milestone: text",
    "tuning.EVOLUTIONARY_PACE_LABELS": "how each pace is worded: text",
    "tuning.LIFE_MILESTONES": "the names of the milestones in order: structure",
    "tuning.LIFE_SCHEDULE_FIELDS": "the column names of the milestone ages: structure",
    "tuning.LIFE_WORLD_STAGES": "which milestones make a world's life stage: structure",
    "habitability.TIERS": "the four colour names: structure",
    "habitability.POTENTIALLY_HABITABLE_TIERS": "which colours count as habitable: structure",
    "habitability.FALLBACK_TIERS": "which colours a fallback request accepts: structure",
    "habitability.FALLBACK_SUITABILITY_TIERS": "which Suitability colours a fallback accepts: structure",
    "habitability.LIFE_STAGES": "the names of the four life stages: structure",
    "habitability.CORE_DOMAINS": "the names of the five core domains: structure",
    "habitability.DOMAINS": "the names of the six domains: structure",
    "habitability.HUMAN": "the human visitor's bands, each edge read from the human limits file with its own source",
    "habitability.ARMSTRONG_KPA": "read from the human limits file (the Armstrong line)",
    "habitability.INSPIRED_O2_MIN_KPA": "read from the human limits file",
    "habitability.INSPIRED_O2_CHRONIC_KPA": "read from the human limits file",
    "habitability.INSPIRED_O2_ACUTE_KPA": "read from the human limits file",
    "habitability.MASK_MIN_KPA": "the sum of two limits read from the human limits file",
    "habitability.WET_BULB_CRIT_C": "read from the human limits file",
    "habitability.PRESSURE_BAND": "read from the human limits file",
    "habitability.TEMPERATURE_BAND": "read from the human limits file",
    "habitability.DOSE_BAND": "read from the human limits file",
    "habitability.UV_BAND": "read from the human limits file",
    "habitability.SINGLE_YEAR_DOSE_SV": "read from the human limits file",
    "habitability.O2_HIGH_BAND": "built from the human limits read from the limits file",
    "habitability_world.TIER_INDEX": "the colour names numbered: structure",
    "habitability_world.TIER_NAMES": "the colour names: structure",
    "habitability_world.DOMAINS": "the domain names: structure",
    "habitability_world.CORE_DOMAINS": "the domain names: structure",
    "habitability_world.BODY_FIELDS": "the stored column names: structure",
    "habitability_world.SCORE_FIELDS": "the stored column names: structure",
    "habitability_world.LIKELIHOOD_FIELDS": "the names of the two likelihoods: structure",
    "habitability_world.COMPACT_HOSTS": "what a compact host is called: text",
    "habitability_world.MOLAR_KG": "standard molar masses of the stored gases",
    "habitability_world.BALANCE_MOLAR_KG": "the molar mass of the balance gas: a standard molar mass",
    "habitability_world.CARBON_KG_MOL": "the molar mass of carbon",
    "habitability_world.NITROGEN_KG_MOL": "the molar mass of nitrogen",
    "tuning.HABITABLE_PLANET_CLASSES": "a legacy list of class codes, retired by the class cut-over (GEN.646)",
}
"""dict: The names of the life code that carry no stand-in value, and why."""


def _rows():
    text = PLAN.read_text()
    block = text.split("<!-- life-standins:begin -->")[1].split("<!-- life-standins:end -->")[0]
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in block.strip().splitlines()]
    return [row for row in rows[2:] if len(row) == 10]


def _named():
    return {f"{module}.{name}" for row in _rows() for module, name in NAME.findall(row[4])}


def _life_names():
    names = {f"tuning.{name}" for name in vars(tuning) if name.startswith(TUNING_PREFIXES)}
    for label, module in (("habitability", habitability), ("habitability_world", habitability_world)):
        names |= {f"{label}.{name}" for name in vars(module) if re.fullmatch(r"[A-Z][A-Z0-9_]*", name)}
    return names


def test_the_registry_has_rows():
    assert len(_rows()) >= 30


@pytest.mark.parametrize("row", _rows(), ids=lambda row: row[0])
def test_a_row_names_a_constant_that_exists(row):
    found = NAME.findall(row[4])
    assert found, f"{row[0]} names no module.CONSTANT"
    for module, name in found:
        assert hasattr(MODULES[module], name), f"{row[0]}: {module}.{name} does not exist"


def test_every_life_constant_has_a_row_or_a_reason_to_have_none():
    missing = sorted(_life_names() - _named() - set(NEEDS_NO_ROW))
    assert not missing, f"a life constant with no stand-in row (GEN.725): add it to section 6.18 of gap-filling-plan.md: {missing}"


def test_no_constant_is_both_in_the_registry_and_on_the_list_of_those_that_need_none():
    assert not sorted(_named() & set(NEEDS_NO_ROW))


def test_every_name_on_the_list_of_those_that_need_none_exists():
    for dotted in NEEDS_NO_ROW:
        module, name = dotted.split(".")
        assert hasattr(MODULES[module], name), dotted
