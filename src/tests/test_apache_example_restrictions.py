# tests/test_apache_example_restrictions.py

"""
The directory restrictions in `examples/apache/planetgen.conf.example`: with
the checkout owned by the web user, Apache must serve nothing from disk but
`static/`, run nothing but the WSGI app, list no directory, and refuse
source, configuration, state, `.git` and dotfiles even if a later `Alias`
points at them.

Runs a real Apache (the example file itself, paths swapped for a throwaway
checkout) on a loopback port; skipped when Apache or mod_wsgi isn't
installed. The server around it is as permissive as it can be (every
directory granted, .htaccess honoured), so the example's own rules are what
is under test. The WSGI app is a stub that echoes the path it was given, so a
200 with `app:` in the body means the request reached the app and never
the file on disk.
"""

import os
import pathlib
import re
import shutil
import signal
import socket
import subprocess
import tempfile
import time
import urllib.error
import urllib.request

import pytest

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EXAMPLE = os.path.join(REPO_DIR, "examples", "apache", "planetgen.conf.example")
MODULES = "/usr/lib/apache2/modules"
APACHE = shutil.which("apache2")

pytestmark = pytest.mark.skipif(
    APACHE is None or not os.path.exists(f"{MODULES}/mod_wsgi.so"),
    reason="needs Apache 2 with mod_wsgi (Debian/Ubuntu layout)",
)

WSGI = (
    "def application(environ, start_response):\n"
    '    start_response("200 OK", [("Content-Type", "text/plain")])\n'
    '    return [("app:" + environ["PATH_INFO"]).encode()]\n'
)
FILES = {
    "src/html/wsgi.py": WSGI,
    "src/html/static/style.css": "css-body",
    "src/html/static/app.js": "js-body",
    "src/html/static/.hidden": "hidden",
    "src/html/static/leftover.py": "py-body",
    "src/html/lib/fmt.py": "source-body",
    "src/planetgen/secret.py": "source-body",
    "config.json": '{"password": "x"}',
    "state/maintenance.json": "{}",
    ".git/config": "[core]",
    ".env": "KEY=value",
}


def _modules(root):
    names = ["mpm_event", "authz_core", "alias", "headers", "mime", "dir", "autoindex",
             "wsgi", "filter", "deflate"]
    return "\n".join(f"LoadModule {n}_module {MODULES}/mod_{n}.so" for n in names)


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture
def apache():
    # Not pytest's tmp_path: /tmp/pytest-of-<user> is private, and the web
    # user has to be able to walk down to the files.
    base = pathlib.Path(tempfile.mkdtemp(prefix="planetgen-apache-"))
    os.chmod(base, 0o755)
    root = base / "planetGen"
    for rel, text in FILES.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    for rel in ("venv/bin", "jobs", "db", "src/html/static/vendor/data"):
        (root / rel).mkdir(parents=True, exist_ok=True)
    (root / "src/html/static/vendor/data/table.js").write_text("nested-data-folder")
    port = _free_port()

    def start(extra=""):
        example = open(EXAMPLE, encoding="utf-8").read().split("# --- HTTPS ---")[0]
        example = example.replace("/var/lib/planetGen", str(root))
        example = re.sub(r" python-home=\S+", "", example)
        example = example.replace("<VirtualHost *:80>", f"<VirtualHost *:{port}>")
        example = example.replace("    Alias /static/", extra.replace("@ROOT@", str(root)) + "    Alias /static/", 1)
        run_as = "User www-data\nGroup www-data\n" if os.geteuid() == 0 else ""
        (base / "httpd.conf").write_text(
            f"ServerRoot {base}\nPidFile {base}/pid\nServerName localhost\nListen 127.0.0.1:{port}\n"
            f"Define APACHE_LOG_DIR {base}\nTypesConfig /etc/mime.types\n{_modules(root)}\n"
            f"WSGISocketPrefix {base}/wsgi\nErrorLog {base}/error.log\nLogLevel warn authz_core:debug\n{run_as}"
            "<Directory />\n    AllowOverride All\n    Require all granted\n</Directory>\n" + example)
        if os.geteuid() == 0:
            subprocess.run(["chown", "-R", "www-data:www-data", str(base)], check=True)
        subprocess.run([APACHE, "-f", str(base / "httpd.conf"), "-k", "start"], check=True,
                       capture_output=True)
        for _ in range(100):
            try:
                socket.create_connection(("127.0.0.1", port), timeout=0.2).close()
                return
            except OSError:
                time.sleep(0.1)
        pytest.fail("Apache did not start: " + (base / "error.log").read_text())

    def get(path):
        request = urllib.request.Request(f"http://127.0.0.1:{port}{path}")
        try:
            with urllib.request.urlopen(request, timeout=10) as response:
                return response.status, response.read().decode()
        except urllib.error.HTTPError as e:
            return e.code, e.read().decode()

    yield start, get, base
    pid = base / "pid"
    if pid.exists():
        try:
            os.kill(int(pid.read_text()), signal.SIGTERM)
        except (OSError, ValueError):
            pass
        time.sleep(0.5)
    shutil.rmtree(base, ignore_errors=True)


def test_the_example_config_starts_apache(apache):
    start, _, _ = apache
    start()


def test_static_files_are_served_and_everything_else_is_the_app(apache):
    start, get, _ = apache
    start()
    assert get("/static/style.css") == (200, "css-body")
    assert get("/static/app.js") == (200, "js-body")
    # A folder in static/ that happens to be named like a data folder is still static.
    assert get("/static/vendor/data/table.js") == (200, "nested-data-folder")
    # Any other URL reaches the app (WSGIScriptAlias /), never a file on disk.
    for path in ("/", "/lib/fmt.py", "/config.json", "/state/maintenance.json", "/.git/config",
                 "/.env", "/src/planetgen/secret.py", "/static/../lib/fmt.py",
                 "/static/%2e%2e/lib/fmt.py"):
        status, body = get(path)
        assert status == 200 and body.startswith("app:"), (path, status, body)


def test_static_lists_no_directory_and_serves_no_source_or_dotfiles(apache):
    start, get, _ = apache
    start()
    for path in ("/static/", "/static/vendor/", "/static/.hidden", "/static/leftover.py"):
        assert get(path)[0] == 403, path


@pytest.mark.parametrize("alias, path", [
    ("/leak/state/", "/leak/state/maintenance.json"),
    ("/leak/git/", "/leak/git/config"),
    ("/leak/src/", "/leak/src/planetgen/secret.py"),
    ("/leak/lib/", "/leak/lib/fmt.py"),
    ("/leak/all/", "/leak/all/config.json"),
    ("/leak/all/", "/leak/all/.env"),
])
def test_a_later_alias_cannot_expose_the_install_tree(apache, alias, path):
    """Defence in depth: the deny rules hold even if someone adds an Alias."""
    start, get, base = apache
    folder = {"/leak/state/": "state/", "/leak/git/": ".git/", "/leak/src/": "src/planetgen/",
              "/leak/lib/": "src/html/lib/", "/leak/all/": ""}[alias]
    start(extra=f"    Alias {alias} @ROOT@/{folder}\n")
    assert get(path)[0] == 403
    # Refused by the configuration (AH01630), not by file permissions.
    assert "AH01630" in (base / "planetgen_error.log").read_text()
