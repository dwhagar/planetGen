"""
ADM.22 in a real browser: a running job's page shows its output in a terminal
fed over Server-Sent Events, keeps it live (new output, state), resumes
without losing or repeating a line when the stream is closed and reopened
(the server ends each response after STREAM_SECONDS), and raises no
Content-Security-Policy violations. Skipped without Playwright, Chromium or
the MySQL test server, like `test_web_a11y.py`.
"""

import json
import os
import threading
import time

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.api.authz import SESSION_COOKIE_NAME  # noqa: E402
from planetgen.web import generate_page, jobs  # noqa: E402

from tests.test_web_a11y import (  # noqa: E402,F401 -- fixtures
    CSP_LISTENER,
    PAGE_ENDPOINTS,
    admin_token,
    base_url,
    browser,
    page_targets,
    sample_job,
    sample_job_tree,
    sample_params,
    site_app,
    site_db,
)


@pytest.fixture
def running_job(tmp_path, monkeypatch):
    root = str(tmp_path / "jobs")
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", root)
    monkeypatch.setattr(jobs, "_runner_alive", lambda *args, **kwargs: True)
    monkeypatch.setattr(generate_page, "STREAM_SECONDS", 1.5)
    job_id = jobs.start_job("reset", "Live job", [{"label": "Step one", "argv": ["true"]}],
                            admin="admin", spawn=False)
    path = os.path.join(root, job_id)

    def write_state(**fields):
        state = {"status": "running", "pid": os.getpid(), "step": 1, "started_at": time.time() - 3}
        state.update(fields)
        with open(os.path.join(path, "state.json"), "w", encoding="utf-8") as f:
            json.dump(state, f)

    def append(text):
        with open(os.path.join(path, "output.log"), "a", encoding="utf-8") as f:
            f.write(text)

    append("line one\n")
    write_state()
    return job_id, append, write_state


def test_the_terminal_follows_a_job_across_reconnects(browser, base_url, admin_token, running_job):
    job_id, append, write_state = running_job
    context = browser.new_context(viewport={"width": 1000, "height": 900}, reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    context.add_init_script(CSP_LISTENER)
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    try:
        page.goto(f"{base_url}/admin/generate/jobs/{job_id}", wait_until="load")
        rows = page.locator(".job-terminal .xterm-rows")
        rows.wait_for()
        _eventually(lambda: "line one" in rows.inner_text())

        # Output written after the first stream has closed arrives through the reconnect.
        time.sleep(2)
        append("line two\r\nbar 50%\rbar 100%\n")
        _eventually(lambda: "line two" in rows.inner_text())
        text = rows.inner_text()
        assert text.count("line one") == 1 and text.count("line two") == 1, text
        assert "bar 100%" in text and "bar 50%" not in text

        write_state(status="succeeded", finished_at=time.time(), exit_code=0)
        badge = page.locator("[data-job-status]")
        _eventually(lambda: badge.inner_text() == "succeeded")
        assert page.locator("[data-job-cancel]").count() == 0
        assert page.locator("a[download]").first.get_attribute("href").endswith(f"/jobs/{job_id}/log")
        assert page.evaluate("window.__cspViolations || []") == []
    finally:
        context.close()
    assert not errors, errors


def _eventually(check, timeout=15):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if check():
            return
        time.sleep(0.1)
    raise AssertionError("condition not reached in time")
