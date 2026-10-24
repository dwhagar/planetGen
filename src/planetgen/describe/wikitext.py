# planetgen/describe/wikitext.py

"""A Document as MediaWiki text, on request (docs/design/architecture/system-description.md, section 4)."""

from planetgen.describe.render import rating_lines, render_document, share_text

TEMPLATES = {"star": "Star Data", "black-hole": "Black Hole Data", "neutron-star": "Neutron Star Data",
             "planet": "Planet Data", "moon": "Class Data"}
"""dict: The wiki's data template for a section kind that has one; the pair kinds are named below."""

PAIR_TEMPLATES = {"Binary system": "Binary System Data", "Wide binary system": "Wide Binary System Data"}


def _wikitable(header, rows):
    lines = ['{| class="wikitable"', "! " + " !! ".join(header)]
    for row in rows:
        lines += ["|-", "| " + " || ".join(row)]
    lines.append("|}")
    return "\n".join(lines)


class Style:
    @staticmethod
    def heading(depth, section):
        marks = "=" * depth
        return f"{marks} {section['title']} {marks}"

    @staticmethod
    def facts(section, block, rows):
        template = TEMPLATES.get(section["kind"]) or PAIR_TEMPLATES.get(block["caption"])
        if template is None:
            return _wikitable(["Property", "Value"], [(label, value) for _key, label, value in rows])
        return "\n".join(["{{" + template] + [f"|{key}={value}" for key, _label, value in rows] + ["}}"])

    @staticmethod
    def composition(block):
        lines = []
        if block["summary"]:
            lines.append(f"'''{block['label']}''': {block['summary']}.")
        if block["parts"]:
            lines.append(_wikitable(["Component", "Share"], [(part["name"].capitalize(), share_text(part)) for part in block["parts"]]))
        return "\n\n".join(lines)

    @staticmethod
    def ratings(block):
        headline, rows = rating_lines(block, "wikitext")
        return headline + "\n\n" + _wikitable(["Category", "Rating", "Score", "Notes"], rows)

    @staticmethod
    def join(chunks):
        return "\n\n".join(chunks)


CATEGORIES = {"system": "Star Systems", "phenomenon": "Phenomena"}
"""dict: The wiki category each kind of Document is filed under."""


def render(document):
    """The Document as wikitext."""
    return render_document(document, Style, "wikitext") + f"\n\n[[Category:{CATEGORIES[document['kind']]}]]\n"
