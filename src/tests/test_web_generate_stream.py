"""
ADM.22: a job's output and state stream to the browser over Server-Sent
Events (`/admin/generate/jobs/<id>/stream`), resume from a byte offset after
a dropped connection, and download in full (`.../log`).
"""

import json
import os
import time

from planetgen.web import jobs

from tests.test_web_generate import client, jobs_root, site  # noqa: F401 -- fixtures


def _job(jobs_root, text, status="succeeded"):
    job_id = jobs.start_job("reset", "Stream", [{"label": "x", "argv": ["true"]}], spawn=False)
    with open(os.path.join(jobs_root, job_id, "output.log"), "wb") as f:
        f.write(text)
    now = time.time()
    state = {"status": status, "step": 1, "started_at": now - 5}
    if status != "running":
        state.update(finished_at=now, exit_code=0 if status == "succeeded" else 1)
    with open(os.path.join(jobs_root, job_id, "state.json"), "w") as f:
        json.dump(state, f)
    return job_id


def _events(response):
    out = []
    for block in response.get_data(as_text=True).split("\n\n"):
        fields = {}
        for line in block.split("\n"):
            if ": " in line and not line.startswith(":"):
                key, value = line.split(": ", 1)
                fields[key] = value
        if "event" in fields:
            out.append(fields)
    return out


def test_a_finished_job_streams_its_log_then_state_then_done(site, client, jobs_root):
    job_id = _job(jobs_root, b"one\r\nbar 50%\rbar 100%\ncaf\xc3\xa9\n")
    resp = client.get(f"/admin/generate/jobs/{job_id}/stream")
    assert resp.mimetype == "text/event-stream"
    assert resp.headers["Cache-Control"] == "no-store"
    events = _events(resp)
    assert [e["event"] for e in events] == ["log", "state", "done"]
    assert json.loads(events[0]["data"])["text"] == "one\r\nbar 50%\rbar 100%\ncafé\n"
    assert events[0]["id"] == str(len(b"one\r\nbar 50%\rbar 100%\ncaf\xc3\xa9\n"))
    assert json.loads(events[1]["data"])["status"] == "succeeded"


def test_a_reconnect_resumes_after_the_last_event_id(site, client, jobs_root):
    job_id = _job(jobs_root, b"first\nsecond\n")
    resp = client.get(f"/admin/generate/jobs/{job_id}/stream", headers={"Last-Event-ID": "6"})
    assert json.loads(_events(resp)[0]["data"])["text"] == "second\n"
    resp = client.get(f"/admin/generate/jobs/{job_id}/stream?offset=13")
    assert [e["event"] for e in _events(resp)] == ["state", "done"]


def test_read_log_never_splits_a_character_and_restarts_past_the_end(jobs_root):
    job_id = _job(jobs_root, "abé".encode())  # 4 bytes: a b 0xC3 0xA9
    assert jobs.read_log(job_id, 0, max_bytes=3) == ("ab", 2)
    assert jobs.read_log(job_id, 2) == ("é", 4)
    assert jobs.read_log(job_id, 99) == ("abé", 4)
    assert jobs.read_log("nosuchjob") is None


def test_a_long_log_replays_only_its_end_from_a_line_start(site, client, jobs_root, monkeypatch):
    monkeypatch.setattr("planetgen.web.generate_page.STREAM_REPLAY_BYTES", 20)
    job_id = _job(jobs_root, b"aaaaaaaaaa\nbbbbbbbbbb\ncccccccccc\n")
    text = json.loads(_events(client.get(f"/admin/generate/jobs/{job_id}/stream"))[0]["data"])["text"]
    assert text == "cccccccccc\n"


def test_the_stream_and_download_need_an_admin_and_a_real_job(site, client, jobs_root):
    job_id = _job(jobs_root, b"x\n")
    site.admin = None
    assert client.get(f"/admin/generate/jobs/{job_id}/stream").status_code == 403
    assert client.get(f"/admin/generate/jobs/{job_id}/log").status_code == 403
    site.admin = {"username": "boss", "must_change_credentials": False}
    assert client.get("/admin/generate/jobs/nosuchjob/stream").status_code == 404
    assert client.get("/admin/generate/jobs/nosuchjob/log").status_code == 404


def test_the_log_downloads_whole_as_text(site, client, jobs_root):
    job_id = _job(jobs_root, b"line\n" * 3)
    resp = client.get(f"/admin/generate/jobs/{job_id}/log")
    assert resp.status_code == 200 and resp.mimetype == "text/plain"
    assert f'filename=job-{job_id}.log' in resp.headers["Content-Disposition"]
    assert resp.get_data() == b"line\n" * 3


def test_only_the_two_terminal_pages_allow_inline_styles(site, client, jobs_root):
    from planetgen import web

    job_id = _job(jobs_root, b"x\n")
    for path in ("/admin/generate", f"/admin/generate/jobs/{job_id}"):
        csp = client.get(path).headers["Content-Security-Policy"]
        assert csp == web.JOB_LOG_CONTENT_SECURITY_POLICY and "script-src" not in csp, path
    for path in ("/", "/admin/generate/status", f"/admin/generate/jobs/{job_id}/log"):
        assert "unsafe-inline" not in client.get(path).headers.get("Content-Security-Policy", ""), path
