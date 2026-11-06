"""static/starthin.js holds the server's per-tile star caps (MAP.366): a tile counts as capped when it lists as many stars as the server lets it."""

import os
import re

from planetgen.db import nodes

SOURCE = os.path.join(os.path.dirname(__file__), "..", "html", "static", "starthin.js")


def _constant(name):
    with open(SOURCE, encoding="utf-8") as handle:
        return int(re.search(rf"export const {name} = (\d+);", handle.read()).group(1))


def test_the_page_caps_are_the_servers():
    assert nodes.NODE_STARS == nodes.NODE_GENERATED_STARS == _constant("NODE_LIST_CAP")
