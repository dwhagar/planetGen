# tests/test_connectionnote.py

"""
The Generate page's "Connection lost" note (ADM.40): the job log stream
ends on purpose every 40 seconds and the browser reconnects at once, so the
note (`html/static/connectionnote.js`) waits out a grace period instead of
flashing at every reconnect, and shows only when the stream stays down.
Checked in Node, with a fake clock, when Node is installed.
"""

import json
import os
import shutil
import subprocess

import pytest

from planetgen.web import generate_page

_STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")

SCRIPT = """
import { quietReconnect } from %(module)s;
let now = 0, nextId = 1; const timers = new Map(); const shown = [];
const fake = {
  setTimeout(fn, ms) { const id = nextId++; timers.set(id, {fn, at: now + ms}); return id; },
  clearTimeout(id) { timers.delete(id); },
};
const advance = (ms) => { now += ms; for (const [id, t] of [...timers]) if (t.at <= now) { timers.delete(id); t.fn(); } };
const link = quietReconnect((text) => shown.push(text), 8000, fake);
const log = [];
// A normal end of stream: the error, then the reconnect a moment later.
link.lost(); advance(2000); link.restored(); advance(60000); log.push(shown.filter((t) => t).length);
// Many in a row, never down for long.
for (let i = 0; i < 5; i++) { link.lost(); advance(2500); link.restored(); }
log.push(shown.filter((t) => t).length);
// A stream that really is down.
link.lost(); link.lost(); advance(7999); log.push(shown.filter((t) => t).length);
advance(2); log.push(shown.filter((t) => t).length);
link.restored(); log.push(shown[shown.length - 1]);
console.log(JSON.stringify(log));
"""


@pytest.mark.skipif(shutil.which("node") is None, reason="Node is not installed")
def test_the_note_waits_out_the_expected_reconnect():
    module = json.dumps("file://" + os.path.abspath(os.path.join(_STATIC, "connectionnote.js")))
    result = subprocess.run(["node", "--input-type=module", "-e", SCRIPT % {"module": module}],
                            capture_output=True, text=True, timeout=60, check=True)
    assert json.loads(result.stdout) == [0, 0, 0, 1, ""]


def test_the_stream_ends_well_inside_the_request_timeout():
    # The grace period in generatejobs.js (8 s) must stay well under the stream length, or a reconnect would
    # still be a visible gap.
    source = open(os.path.join(_STATIC, "generatejobs.js"), encoding="utf-8").read()
    assert "quietReconnect(connection, RECONNECT_GRACE_MS)" in source
    assert generate_page.STREAM_SECONDS >= 20
