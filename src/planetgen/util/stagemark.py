# planetgen/util/stagemark.py

"""
Stage markers (UX.98)
=====================

The one wording for "where a job is": `"Stage 4 of 12"` for one stage and
`"Stages 4 to 12 of 12"` for a run of them. The web log's section headers,
the stage lines the commands print, the CLI progress bar, the status panel
and the downloaded log all take their text from here (the page's script shows
the `step_marker` the server sends), so they agree.

Standard library only, so the job runner can import it before anything else.
"""


def run(first, count=1):
    """`"Stage 4"`, or `"Stages 4 to 12"` for `count` stages from `first`."""
    return f"Stage {first}" if count == 1 else f"Stages {first} to {first + count - 1}"


def marker(first, total, count=1):
    """`"Stage 4 of 12"`, or `"Stages 4 to 12 of 12"`."""
    return f"{run(first, count)} of {total}"


def section(first, total, count, label):
    """A section separator line: `"=== Stages 4 to 12 of 12: label ==="`."""
    return f"=== {marker(first, total, count)}: {label} ==="


def finished(first, count, code, seconds):
    """The line closing a section: `"=== Stages 4 to 12 exited with status 0 after 31 s ==="`."""
    return f"=== {run(first, count)} exited with status {code} after {seconds:.0f} s ==="
