"""
The page of life class rules is the rule file as built (GEN.886, docs/design/habitability/life-class-rules.md,
scripts/build_life_class_rules.py).
"""

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_life_class_rules as builder  # noqa: E402


def test_the_page_is_the_rule_file_as_built():
    assert pathlib.Path(builder.OUT).read_text(encoding="utf-8") == builder.render(builder.load()), \
        "run python scripts/build_life_class_rules.py to rebuild docs/design/habitability/life-class-rules.md"


def test_every_class_of_the_rule_file_has_a_section_and_a_row():
    document = builder.load()
    text = builder.render(document)
    for row in document["classes"]:
        assert f"### Class {row['id']}: " in text and f"| {row['id']} |" in text
