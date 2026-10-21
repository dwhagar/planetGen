// tests/js/adminjobstrip.test.mjs -- static/adminjobstrip.js, the admin job strip (ADM.119, ADM.123): how the status
// JSON's `strip` fills the running job and the one waiting behind it, and what shows when the job ends.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const { apply } = await import(new URL("adminjobstrip.js", STATIC_URL).href);

function parts() {
  const node = () => ({ hidden: false, textContent: "", value: 0, attrs: {},
    setAttribute(name, value) { this.attrs[name] = value; }, removeAttribute(name) { delete this.attrs[name]; } });
  const p = {};
  for (const name of ["root", "run", "title", "text", "bar", "elapsed", "cancel", "jobId", "queued", "queuedTitle", "queuedNote", "start"]) p[name] = node();
  return p;
}

const job = { id: "j1", title: "Add sectors", url: "/admin/generate/jobs/j1", text: "Stage 2 of 3", elapsed_text: "3 m 05 s", value: 0.4 };

test("a running job fills the strip and shows it", () => {
  const p = parts();
  p.root.hidden = true;
  assert.equal(apply({ job, queued: null }, p, null), job);
  assert.equal(p.root.hidden, false);
  assert.equal(p.run.hidden, false);
  assert.equal(p.queued.hidden, true);
  assert.equal(p.title.textContent, "Add sectors");
  assert.equal(p.title.attrs.href, "/admin/generate/jobs/j1");
  assert.equal(p.jobId.value, "j1");
  assert.equal(p.bar.value, 0.4);
});

test("a job that does not know how far it is gets an unmoving bar", () => {
  const p = parts();
  p.bar.attrs.value = "0.2";
  apply({ job: { ...job, value: null }, queued: null }, p, null);
  assert.equal("value" in p.bar.attrs, false);
});

test("the waiting job shows with its note, and Start now only when it is held", () => {
  const p = parts();
  apply({ job, queued: { title: "Redo scatters", admin: "boss", held: null } }, p, null);
  assert.equal(p.queued.hidden, false);
  assert.equal(p.queuedTitle.textContent, "Redo scatters");
  assert.equal(p.queuedNote.textContent, "Starts when the running job ends.");
  assert.equal(p.start.hidden, true);
  apply({ job: null, queued: { title: "Redo scatters", admin: "boss", held: "The job before it did not succeed (failed)." } }, p, job);
  assert.equal(p.start.hidden, false);
  assert.equal(p.queuedNote.textContent, "The job before it did not succeed (failed).");
});

test("nothing running or waiting hides the strip, but a job that just ended says it finished", () => {
  const p = parts();
  assert.equal(apply(null, p, null), null);
  assert.equal(p.root.hidden, true);
  const q = parts();
  assert.equal(apply(null, q, job), null);
  assert.equal(q.root.hidden, false);
  assert.equal(q.text.textContent, "Finished.");
  assert.equal(q.cancel.hidden, true);
});

test("only a waiting job keeps the strip open without a running job", () => {
  const p = parts();
  apply({ job: null, queued: { title: "x", admin: "boss", held: "held" } }, p, null);
  assert.equal(p.root.hidden, false);
  assert.equal(p.run.hidden, true);
});
