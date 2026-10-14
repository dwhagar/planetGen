"""The docs site (GitHub Pages, Jekyll) reads every Markdown page under docs/ as a Liquid template, so a literal
`{%` or `{{` outside a raw block (a Jinja example, say) stops the whole site from building (OPS.44)."""

import glob
import os
import re

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RAW_BLOCK = re.compile(r"\{%-?\s*raw\s*-?%\}.*?\{%-?\s*endraw\s*-?%\}", re.S)


def _liquid_text(text):
    """`text` with each raw block cut out, the line breaks kept so a hit keeps its line number."""
    return RAW_BLOCK.sub(lambda match: "\n" * match.group(0).count("\n"), text)


def test_no_docs_page_holds_liquid_delimiters_outside_a_raw_block():
    pages = glob.glob(os.path.join(ROOT, "docs", "**", "*.md"), recursive=True)
    found = []
    for path in sorted(pages):
        with open(path, encoding="utf-8") as f:
            for number, line in enumerate(_liquid_text(f.read()).split("\n"), 1):
                if "{%" in line or "{{" in line:
                    found.append(f"{os.path.relpath(path, ROOT)}:{number}")
    assert not found, ("Jekyll reads these as Liquid and fails the Pages build; wrap the example in "
                       f"raw/endraw tags or reword it: {found}")


def test_a_raw_block_hides_its_delimiters_and_nothing_else():
    text = "a {{ x }}\n{% raw %}b {{ y }}\n{% endraw %}\nc {% tag %}"
    left = _liquid_text(text).split("\n")
    assert "{{" in left[0] and "{{" not in left[1] and "{%" in left[3]
