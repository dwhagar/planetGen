# planetgen/describe/text.py

"""A Document as plain text for a terminal, with no markup (docs/design/system-description.md, section 4)."""

from planetgen.describe.render import rating_lines, render_document, share_text


def _table(header, rows):
    rows = [tuple(row) for row in rows]
    widths = [max(len(str(line[column])) for line in [header] + rows) for column in range(len(header))]
    lines = ["  ".join(str(cell).ljust(width) for cell, width in zip(line, widths)).rstrip() for line in [header] + rows]
    lines.insert(1, "  ".join("-" * width for width in widths))
    return "\n".join(lines)


class Style:
    @staticmethod
    def heading(depth, section):
        title = section["title"]
        if depth <= 2:
            return f"{title}\n{('=' if depth == 1 else '-') * len(title)}"
        return title

    @staticmethod
    def facts(section, block, rows):
        return _table(("Property", "Value"), [(label, value) for _key, label, value in rows])

    @staticmethod
    def composition(block):
        lines = []
        if block["summary"]:
            lines.append(f"{block['label']}: {block['summary']}.")
        if block["parts"]:
            lines.append(_table(("Component", "Share"), [(part["name"].capitalize(), share_text(part)) for part in block["parts"]]))
        return "\n\n".join(lines)

    @staticmethod
    def ratings(block):
        headline, rows = rating_lines(block, "text")
        return headline + "\n\n" + _table(("Category", "Rating", "Score", "Notes"), rows)

    @staticmethod
    def join(chunks):
        return "\n\n".join(chunks) + "\n"


def render(document):
    """The Document as plain text."""
    return render_document(document, Style, "text")
