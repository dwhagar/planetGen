"""OPS.299: a site that does not answer /api/health for a while makes the update say what it can see (the site's processes, Apache's error log)."""

import os
import subprocess

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _wait(tmp_path, curl_body, wait="4"):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    curl = bin_dir / "curl"
    curl.write_text("#!/bin/bash\n" + curl_body)
    curl.chmod(0o755)
    env = dict(os.environ, PATH=f"{bin_dir}:{os.environ['PATH']}", HEALTH_WAIT_SECONDS=wait, HEALTH_SETTLE_SECONDS="0", HEALTH_SLOW_NOTE_SECONDS="1")
    return subprocess.run(["bash", "-c", f'source "{REPO_DIR}/scripts/deploy-common.sh"; wait_for_health http://x/api/health "the API" "config.json api_base_url"'],
                          env=env, capture_output=True, text=True, timeout=60)


def test_a_site_that_never_answers_gets_one_note_of_what_runs_and_the_warning(tmp_path):
    proc = _wait(tmp_path, "exit 28\n")
    assert proc.returncode == 1
    assert proc.stdout.count("the site's processes:") == 1
    assert "the health check of the API failed" in proc.stderr + proc.stdout


def test_a_site_that_answers_late_says_how_late(tmp_path):
    counter = tmp_path / "count"
    body = f'n=$(cat {counter} 2>/dev/null || echo 0); echo $((n+1)) > {counter}\nif (( n < 3 )); then exit 28; fi\nprintf \'ok\\n200\'\n'
    proc = _wait(tmp_path, body, wait="30")
    assert proc.returncode == 0
    assert "the site answered" in proc.stdout and "the site's processes:" in proc.stdout


def test_a_credential_in_an_error_log_line_is_masked_in_the_update_log(tmp_path):
    error_log = tmp_path / "error.log"
    error_log.write_text("[wsgi:error] redis.exceptions.ConnectionError: redis://planetgen:SECRETPASSWORD123@127.0.0.1:6379/0\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    curl = bin_dir / "curl"
    curl.write_text("#!/bin/bash\nexit 28\n")
    curl.chmod(0o755)
    env = dict(os.environ, PATH=f"{bin_dir}:{os.environ['PATH']}", HEALTH_WAIT_SECONDS="3", HEALTH_SETTLE_SECONDS="0", HEALTH_SLOW_NOTE_SECONDS="1",
               HEALTH_ERROR_LOGS=str(error_log))
    proc = subprocess.run(["bash", "-c", f'source "{REPO_DIR}/scripts/deploy-common.sh"; wait_for_health http://x/api/health "the API" "config.json api_base_url"'],
                          env=env, capture_output=True, text=True, timeout=60)
    text = proc.stdout + proc.stderr
    assert "Apache's error log:" in text and "ConnectionError" in text
    assert "SECRETPASSWORD123" not in text
