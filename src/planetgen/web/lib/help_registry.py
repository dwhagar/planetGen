# planetgen/web/lib/help_registry.py

"""
The help registry: each part of the program registers its own help pages, categories and link terms here, and nothing else keeps a list.
The Help menu (`web/helpers.py`), the help index (`/help`), the sitemap and the coverage test (`tests/test_help_registry.py`) read what has
been registered, so a page cannot be in one place and missing from another, and a new part of the program joins the help by registering,
with no edit to a central list.

A module registers when it is imported (the web app imports them all at start):

    register_category("objects", "Objects and classes", "What the objects of the galaxy are.", order=10)
    register_topic("classes", "Classes", "objects", "Every class the survey assigns.", order=10, children=class_children)
    register_terms("life", life_terms)      # `life_terms()` returns [Term(text, endpoint, values)] for the auto-linker

    topics()            every `Topic`, by `order`
    categories()        `[{"slug", "title", "summary", "topics": [Topic]}]`, by `order`, a category with no topic left out
    menu()              `((endpoint, label), ...)` for the Help menu
    children(topic)     `[Child]` of a topic, from the function its module gave
    terms()             every `Term` the registered providers give, for linking descriptions to their help pages

Registering the same slug or endpoint again replaces the earlier entry (a module reloaded in a test), so nothing doubles.
Imports no Flask, and nothing heavy.
"""

from collections import namedtuple

Topic = namedtuple("Topic", "endpoint label category summary order children")
"""namedtuple: A help page. `children` is a function returning its `Child` list, or `None`."""

Child = namedtuple("Child", "endpoint values label summary")
"""namedtuple: A page under a topic: the endpoint, its URL values (a dict), its label and a one-line summary ('' for none)."""

Term = namedtuple("Term", "text endpoint values after", defaults=("",))
"""namedtuple: A word or phrase a description may mention, and the help page (endpoint and URL values) it links to. `after` is text that must come right before it (a bare number is a term only after "Class "); the linker leaves it outside the link."""

Category = namedtuple("Category", "slug title summary order")

INDEX_ENDPOINT = "help_index"
"""str: The endpoint of the help index page."""

_categories = {}
_topics = {}
_term_sources = {}


def register_category(slug, title, summary, order=100):
    """Registers (or replaces) a help category."""
    _categories[slug] = Category(slug, title, summary, order)


def register_topic(endpoint, label, category, summary, order=100, children=None):
    """Registers (or replaces) a help page: its Flask endpoint, label, category slug, one-line summary, place in the menu and, for a page with
    pages beneath it, the function that lists them (`[Child]`, read from the program's tables when called)."""
    _topics[endpoint] = Topic(endpoint, label, category, summary, order, children)


def register_terms(source, provider):
    """Registers (or replaces) a provider of link terms: `provider()` returns `[Term]`, read from the tables the module owns."""
    _term_sources[source] = provider


def topics():
    """Every registered `Topic`, by `order` and then label."""
    return tuple(sorted(_topics.values(), key=lambda t: (t.order, t.label)))


def menu():
    """What the Help menu opens: `((endpoint, label), ...)`, the topics in order."""
    return tuple((topic.endpoint, topic.label) for topic in topics())


def categories():
    """The categories with their topics: `[{"slug", "title", "summary", "topics": [Topic]}]`, in order. A category with no topic is left out;
    a topic naming a category nobody registered is a mistake and raises `KeyError`."""
    for topic in _topics.values():
        if topic.category not in _categories:
            raise KeyError(f"help topic {topic.endpoint!r} names the unregistered category {topic.category!r}")
    found = [{"slug": c.slug, "title": c.title, "summary": c.summary, "topics": [t for t in topics() if t.category == c.slug]}
             for c in sorted(_categories.values(), key=lambda c: (c.order, c.title))]
    return [entry for entry in found if entry["topics"]]


def children(topic):
    """The pages under `topic`, read from the program's tables; `[]` for a topic with none."""
    return list(topic.children()) if topic.children else []


def terms():
    """Every `Term` the registered providers give, in registration order."""
    return [term for provider in _term_sources.values() for term in provider()]
