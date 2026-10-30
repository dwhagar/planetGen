# planetgen/web/job_steps.py

"""
The steps of a Generate job as they run (PERF.203 step 4, docs/design/operations/one-run-per-job.md).

The Generate page builds a job as a list of steps, one command line each (`generate_page.build_job`): that list is what the Review
page shows and what the size estimate reads. When the job starts (`jobs.start_job`) the steps that are generation commands a run can
hold (`generation.job.JOB_COMMANDS`) become the one step `planetgen run-job SPEC`, so the job's generating is ONE run in the work
queue with one stage list and one estimate: the data check and the reset stay steps of their own, and a step that is not a stage of
a run (the band, `finish-run`, a database check) ends the group it stood in.
"""

import json

from planetgen.web import jobs

SHARED_FLAGS = ("--disregard-missing-data", "--debug")
"""tuple: The flags the page adds to every generating step; the one `run-job` takes them once for all its commands."""

REPEATED_STAGE = "Make the block light totals"
"""str: The label of the stage that the scatters each end with and a run does once, after the last of them (`followed.dedupe`)."""


def entry_of(step):
    """
    `(command, options, flags)` of a step that is a generation command a run can hold, else `None`: the command's name, the
    options after it and the `SHARED_FLAGS` the page added.
    """
    from planetgen.generation import job

    argv = step["argv"]
    start = len(jobs.GENERATE_COMMAND) + 1
    if len(argv) <= start or argv[1:start] != jobs.GENERATE_COMMAND or argv[start] not in job.JOB_COMMANDS:
        return None
    command = argv[start]
    options = [word for word in argv[start + 1:] if word not in SHARED_FLAGS]
    flags = [word for word in argv[start + 1:] if word in SHARED_FLAGS]
    if command == "scatter" and ({"--band", "--nebula-field"} & set(options)):
        return None
    if command == "settle" and "--sector-id" in options:
        return None
    if command == "sectors" and "--region" in options:
        return None
    return command, options, flags


def merged_stages(steps):
    """The stage list of a group of steps: each step's own, the stage every scatter ends with only at its last place."""
    listed = []
    for step in steps:
        listed += [dict(stage) for stage in step.get("stages") or [{"label": step["label"], "skipped": None}]]
    last = max((place for place, stage in enumerate(listed) if stage["label"] == REPEATED_STAGE), default=None)
    return [stage for place, stage in enumerate(listed) if stage["label"] != REPEATED_STAGE or place == last]


def collapse(steps, title):
    """
    `steps` as they run: every group of two or more neighbouring steps that are generation commands a run can hold is one step
    `planetgen run-job SPEC` (labelled `title`, or "Rest of the job" when the job has other generating steps outside the group),
    its stages the group's. Any other step, and a group of one (already the one run of its own command), is left as it is.
    """
    entries = [entry_of(step) for step in steps]
    out, place = [], 0
    groups = []
    while place < len(steps):
        end = place
        while end < len(steps) and entries[end] is not None:
            end += 1
        if end - place >= 2:
            groups.append((place, end))
            place = end
        else:
            place = max(end, place + 1)
    if not groups:
        return list(steps)
    whole = len(groups) == 1 and sum(1 for entry in entries if entry is not None) == groups[0][1] - groups[0][0]
    place = 0
    for first, end in groups:
        out += steps[place:first]
        group = steps[first:end]
        spec = [{"command": entries[at][0], "argv": entries[at][1]} for at in range(first, end)]
        flags = list(dict.fromkeys(flag for at in range(first, end) for flag in entries[at][2]))
        out.append({"label": title if whole else "Rest of the job",
                    "argv": [group[0]["argv"][0], *jobs.GENERATE_COMMAND, "run-job", json.dumps(spec), *flags],
                    "stages": merged_stages(group)})
        place = end
    return out + steps[place:]
