"""OPS.226: the memory the template saves, measured on the real sampling tables (stage 0 measured it once, `reports/job-engine-design-2026-11-13.md` section 10.4): jobs forked from one template share the
tables it built, where as many workers would each build their own. Slow: the tables take about a minute to build, and the test builds them with no kept copy, which is the size a worker with none holds."""

import os
import time

import pytest

from planetgen.queue import allowlist, jobhost

pytestmark = [pytest.mark.slow, pytest.mark.skipif(not os.path.exists("/proc/self/smaps_rollup"), reason="needs /proc/<pid>/smaps_rollup (Linux)")]

R = "tests.template_runners"
JOBS = 3


def _rollup_mb(pid, field):
    with open(f"/proc/{pid}/smaps_rollup", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(field + ":"):
                return int(line.split()[1]) / 1024.0


def test_jobs_forked_from_one_template_hold_the_star_tables_once(monkeypatch):
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    monkeypatch.setenv("PLANETGEN_BAND_TABLE_CACHE", "")        # the tables are built, as a worker with no kept copy builds them (the suite keeps them, and the kept form is a fraction of the size)
    templates = jobhost.TemplateSet(ready_seconds=900)
    try:
        held = [templates.fork("scatter-mass", f"{R}.draws_stars", {"seed": 100 + number, "draws": 1000})[1] for number in range(JOBS)]
        host = templates.host("tables")
        said, deadline = {}, time.monotonic() + 600
        while len(said) < JOBS and time.monotonic() < deadline:
            for job in held:
                for message in job.messages():
                    said[job.id] = message
            time.sleep(0.1)
        assert len(said) == JOBS and all("pss_mb" in message for message in said.values()), f"every job drew its stars: {said}"
        template_mb = _rollup_mb(host.template_pid, "Rss")          # what the template holds: the code and the tables, whoever shares them
        family_mb = _rollup_mb(host.template_pid, "Pss") + sum(message["pss_mb"] for message in said.values())       # what the template and its jobs hold together, shared pages counted once
        separate_mb = JOBS * template_mb          # as many workers, each holding the tables it built for itself (the tables are most of what the template holds)
        assert template_mb > 300, f"the tables are in the template ({template_mb:.0f} MB)"
        assert family_mb < 0.6 * separate_mb, f"the template and {JOBS} jobs hold {family_mb:.0f} MB, {JOBS} workers would hold about {separate_mb:.0f} MB"
        for job in held:
            job.send({"op": "stop"})
    finally:
        templates.stop()
