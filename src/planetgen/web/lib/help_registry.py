# planetgen/web/lib/help_registry.py

"""
The one list of help pages (`/help`): every page under the Help menu is declared here once, with its category, its label and a
one-line summary. The Help menu (`web/helpers.py`), the help index page (`/help`), the sitemap and the coverage test
(`tests/test_help_registry.py`) all read this list, so a page cannot be in one place and missing from another.

The children of a page (the class types under Classes, the life topics and life classes under Life) are not written here: each topic names a
function that reads them from the tables the program uses (`web/lib/classref.py`, `web/lib/life_help.py`), so a new class appears in the index
the moment the program knows it. Imports no Flask, and nothing heavy at import time.

    topics()            every `Topic`, in the order of the menu
    categories()        `[{"slug", "title", "summary", "topics": [Topic]}]`, in the order of `CATEGORIES`
    menu()              `((endpoint, label), ...)`, what the Help menu opens
    children(topic)     `[Child]` of a topic, read from the program's tables
"""

from collections import namedtuple

Topic = namedtuple("Topic", "endpoint label category summary children")
"""namedtuple: A help page. `children` is the name of a function of this module that lists the pages under it, or `None`."""

Child = namedtuple("Child", "endpoint values label summary")
"""namedtuple: A page under a topic: the endpoint, its URL values (a dict), its label and a one-line summary ('' for none)."""

CATEGORIES = (
    ("objects", "Objects and classes", "What the objects of the galaxy are, how they are sorted, and where their values come from."),
    ("worlds", "Worlds and life", "How a visitor's gear is judged on a world and how life on it is described."),
)
"""tuple: `(slug, title, summary)` of each help category, in the order the index lists them."""

TOPICS = (
    Topic("classes", "Classes", "objects",
          "Every class of star, planet, nebula and other object the survey assigns.", "class_children"),
    Topic("habitability_levels", "Habitability", "worlds",
          "What the equipment levels, the Blue to Red domain colors and the Habitable chips mean.", None),
    Topic("property_help_page", "Properties", "objects",
          "Where the orbit, heat, spin, field, ring and life values on a star, planet or moon page come from.", None),
    Topic("life_help_page", "Life", "worlds",
          "The kinds of life, the scales a planet page uses, the stages and the life classes.", "life_children"),
)
"""tuple: Every help page, in the order of the Help menu. The endpoint is a function of `web/class_pages.py`."""

INDEX_ENDPOINT = "help_index"
"""str: The endpoint of the help index page."""


def topics():
    """Every `Topic`, in the order of the menu."""
    return TOPICS


def menu():
    """What the Help menu opens: `((endpoint, label), ...)`, the topics in order."""
    return tuple((topic.endpoint, topic.label) for topic in TOPICS)


def categories():
    """The categories with their topics: `[{"slug", "title", "summary", "topics": [Topic]}]`. A category with no topic is left out."""
    found = [{"slug": slug, "title": title, "summary": summary, "topics": [t for t in TOPICS if t.category == slug]}
             for slug, title, summary in CATEGORIES]
    return [entry for entry in found if entry["topics"]]


def class_children():
    """The class types, alphabetically, each with how many classes it has (`web/lib/classref.py`)."""
    from planetgen.web.lib import classref
    rows = classref.alphabetical([{"name": entry["name"], "slug": slug, "count": len(entry["classes"]), "summary": entry["summary"]}
                                  for slug, entry in classref.catalog().items()])
    return [Child("class_type_page", {"type_slug": row["slug"]}, row["name"], row["summary"]) for row in rows]


def life_children():
    """The life help pages (`life_help.TOPICS`, the first being the Life page itself) and one page for each life class number."""
    from planetgen.web.lib import life_help
    found = [Child("life_help_topic", {"topic": slug}, title, summary) for slug, title, summary in life_help.TOPICS if slug != "kinds"]
    found += [Child("life_help_class", {"number": number}, label, summary) for number, label, summary in life_help.class_numbers()]
    return found


def children(topic):
    """The pages under `topic`, read from the program's tables; `[]` for a topic with none."""
    return globals()[topic.children]() if topic.children else []
