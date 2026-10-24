# tests/render_support.py

"""The text a system is written as, for the tests that compare one system's text with another's or look for a line in it."""

from planetgen import describe


def text(system):
    """The system's plain text (`describe.render` of its Document)."""
    return describe.render(describe.describe(system), "text")


def prose_text(section):
    """All the sentences of a Document section's prose blocks, in one string."""
    return " ".join(describe.quantity.inline_text(sentence) for block in section["blocks"] if block["type"] == "prose"
                    for sentence in block["text"])


def section_prose(system, name):
    """The prose of the section titled `name` in the system's Document (a star, planet, moon, belt or comet)."""
    document = describe.describe(system)
    found = next(section for section in describe.model.sections(document) if section["title"] == name)
    return prose_text(found)


def facts_text(section):
    """The labels and values of a Document section's fact rows, in one string."""
    return " ".join(f"{row['label']} {describe.quantity.inline_text(row['value'])}" for block in section["blocks"]
                    if block["type"] == "facts" for row in block["rows"])


def phenomenon_section(obj):
    """The Document section of a phenomenon (`describe.describe_phenomenon`)."""
    return describe.describe_phenomenon(obj)["sections"][0]
