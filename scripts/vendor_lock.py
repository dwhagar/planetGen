"""Write or check src/html/static/vendor/VENDORED.json (OPS.35), the lock of the libraries the site ships itself.

Per package: its version, the npm `dist.integrity` of the release it was taken from, the tool it was built with (the
bundles are made with esbuild), its licence, the libraries of other authors bundled into it, and the SHA-256 of
every file shipped for it. `tests/test_vendored_lock.py` hashes every file under `vendor/` against it, so an
accidental edit, a stale file or an unlisted file fails a test, and a reviewer sees in the lock's diff exactly what a
vendor bump changed.

Usage:

    python scripts/vendor_lock.py            # rewrite the lock (a vendor bump; run after the vendor script)
    python scripts/vendor_lock.py --check    # exit 1 and say what differs, change nothing

The integrity strings come from `npm view <package>@<version> dist.integrity` (needs npm and the network); a package
already in the lock at the same version keeps its recorded string, so a rewrite for one bump asks only for that
package. `scripts/vendor_shoelace.py` calls this when it has copied Shoelace.
"""

import hashlib
import json
import os
import subprocess
import sys

VENDOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "html", "static", "vendor")
LOCK = os.path.join(VENDOR, "VENDORED.json")
NOT_LOCKED = ("VENDORED.json", "THIRD_PARTY_NOTICES.txt")
"""tuple: Files under vendor/ that are not libraries."""

ESBUILD = "esbuild 0.28.2"
"""str: The bundler of the files that are bundled here; `npx esbuild@0.28.2 --bundle --minify --format=esm` makes
three.module.min.js and both TanStack files byte for byte (checked 2026-10-24)."""

COPIED = "not rebuilt: the files are copied from the npm package"

PACKAGES = (
    {"package": "three", "version": "0.186.0", "license": "MIT", "built_with": ESBUILD,
     "files": ("three.module.min.js", "three-LICENSE.txt")},
    {"package": "@shoelace-style/shoelace", "version": "2.20.1", "license": "MIT", "built_with": COPIED,
     "files": ("shoelace/",),
     "bundled": (
         {"package": "lit, lit-html, lit-element, @lit/reactive-element", "license": "BSD-3-Clause",
          "note": "Copyright Google LLC, in the license banners of the vendored chunks"},
         {"package": "@floating-ui/dom, @floating-ui/core", "license": "MIT", "note": "chunks/chunk.NDC2Y3JV.js"},
         {"package": "@shoelace-style/localize", "license": "MIT", "note": "the localize controller of the chunks"},
     )},
    {"package": "@xterm/xterm", "version": "6.0.0", "license": "MIT", "built_with": COPIED,
     "files": ("xterm/xterm.mjs", "xterm/xterm.css", "xterm/LICENSE.txt")},
    {"package": "@xterm/addon-fit", "version": "0.11.0", "license": "MIT", "built_with": COPIED,
     "files": ("xterm/addon-fit.mjs",)},
    {"package": "@tanstack/table-core", "version": "8.21.3", "license": "MIT", "built_with": ESBUILD,
     "files": ("tanstack/table-core.min.js", "tanstack/LICENSE-table-core.txt")},
    {"package": "@tanstack/virtual-core", "version": "3.13.12", "license": "MIT", "built_with": ESBUILD,
     "files": ("tanstack/virtual-core.min.js", "tanstack/LICENSE-virtual.txt")},
)
"""tuple: The libraries shipped, and the files of each (a name ending in `/` is a folder)."""


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def shipped_files():
    """Every library file under vendor/, as paths relative to it with `/`, sorted."""
    found = []
    for folder, _dirs, names in os.walk(VENDOR):
        for name in names:
            relative = os.path.relpath(os.path.join(folder, name), VENDOR).replace(os.sep, "/")
            if relative not in NOT_LOCKED:
                found.append(relative)
    return sorted(found)


def files_of(entry, shipped):
    """The shipped files an entry covers, `{path: sha256}`."""
    paths = [p for want in entry["files"] for p in (
        [q for q in shipped if q.startswith(want)] if want.endswith("/") else [want])]
    return {path: sha256_of(os.path.join(VENDOR, path)) for path in sorted(paths)}


def integrity_of(entry, old):
    """The npm integrity string for the release: the recorded one when the lock has this version, else the registry's."""
    for package in old.get("packages", []):
        if package["package"] == entry["package"] and package["version"] == entry["version"] and package.get("integrity"):
            return package["integrity"]
    spec = f"{entry['package']}@{entry['version']}"
    return subprocess.run(["npm", "view", spec, "dist.integrity"], check=True, capture_output=True,
                          text=True).stdout.strip()


def build(old=None):
    """The lock as a dict."""
    old = old or {}
    shipped = shipped_files()
    packages = []
    for entry in PACKAGES:
        item = {"package": entry["package"], "version": entry["version"], "integrity": integrity_of(entry, old),
                "built_with": entry["built_with"], "license": entry["license"]}
        if entry.get("bundled"):
            item["bundled"] = [dict(b) for b in entry["bundled"]]
        item["files"] = files_of(entry, shipped)
        packages.append(item)
    return {"note": "Written by scripts/vendor_lock.py (OPS.35); tests/test_vendored_lock.py checks it.",
            "packages": packages}


def read_lock():
    try:
        with open(LOCK, encoding="utf-8") as handle:
            return json.load(handle)
    except FileNotFoundError:
        return {}


def render(lock):
    return json.dumps(lock, indent=2, sort_keys=False) + "\n"


def refresh():
    """Rewrites the lock."""
    text = render(build(read_lock()))    # before the file is opened for writing: it reads the old lock
    with open(LOCK, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def differences():
    """What the files under vendor/ differ in from the lock: a list of sentences (empty when they match)."""
    lock = read_lock()
    locked = {path: digest for package in lock.get("packages", []) for path, digest in package["files"].items()}
    problems = []
    for path in shipped_files():
        if path not in locked:
            problems.append(f"{path} is shipped but is in no package of the lock")
        elif sha256_of(os.path.join(VENDOR, path)) != locked[path]:
            problems.append(f"{path} differs from the lock")
    for path in sorted(set(locked) - set(shipped_files())):
        problems.append(f"{path} is in the lock but not shipped")
    return problems


if __name__ == "__main__":
    if "--check" in sys.argv[1:]:
        found = differences()
        print("\n".join(found) if found else "the lock matches")
        sys.exit(1 if found else 0)
    refresh()
    print(f"wrote {os.path.normpath(LOCK)}")
