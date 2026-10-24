# planetgen/describe/markdown.py

"""A Document as Markdown, on request (docs/design/architecture/system-description.md, section 4)."""

from planetgen.describe.quantity import inline_text
from planetgen.describe.render import rating_lines, render_document, share_text


def _cell(text):
    return text.replace("|", "\\|")


def _table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(_cell(cell) for cell in row) + " |" for row in rows]
    return "\n".join(lines)


class Style:
    @staticmethod
    def heading(depth, section):
        return f"{'#' * depth} {section['title']}"

    @staticmethod
    def facts(section, block, rows):
        return _table(["Property", "Value"], [(label, value) for _key, label, value in rows])

    @staticmethod
    def composition(block):
        lines = []
        if block["summary"]:
            lines.append(f"**{block['label']}**: {block['summary']}.")
        if block["parts"]:
            lines.append(_table(["Component", "Share"], [(part["name"].capitalize(), share_text(part)) for part in block["parts"]]))
        return "\n\n".join(lines)

    @staticmethod
    def ratings(block):
        headline, rows = rating_lines(block, "markdown")
        return headline + "\n\n" + _table(["Category", "Rating", "Score", "Notes"], rows)

    @staticmethod
    def join(chunks):
        return "\n\n".join(chunks) + "\n"


def render(document):
    """The Document as Markdown."""
    return render_document(document, Style, "markdown")
