// tests/js/adminjobstrip.test.mjs -- static/adminjobstrip.js, the admin job strip (ADM.119, ADM.123, PERF.154): how the
// status JSON's `strip` fills the running job and the waiting list behind it, and what shows when the job ends.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL, installDom } from "./fakedom.mjs";

installDom("http://site.test/admin");
const { apply } = await import(new URL("adminjobstrip.js", STATIC_URL).href);

function parts() {
  const node = () => ({ hidden: false, textContent: "", value: 0, attrs: {},
    setAttribute(name, value) { this.attrs[name] = value; }, removeAttribute(name) { delete this.attrs[name]; } });
  const p = {};
  for (const name of ["run", "title", "text", "bar", "elapsed", "cancel", "jobId", "queued"]) p[name] = node();
  p.root = { ...node(), dataset: { actionUrl: "/admin/generate", csrfField: "csrf_token", csrfToken: "tok" } };
  p.line = document.createElement("ol");
  return p;
}

const job = { id: "j1", title: "Add sectors", url: "/admin/generate/jobs/j1", text: "Stage 2 of 3", elapsed_text: "3 m 05 s", value: 0.4 };
const waiting = { ticket: "t1", type: "job", title: "Redo scatters", who: "boss", held: null, position: 1,
  note: "Starts when the running job ends.", startable: false, discardable: true };

// The entry's list item as text, and the actions of its forms.
function items(p) {
  return p.line.children.map((li) => ({
    texts: li.children.filter((c) => c.tagName !== "FORM").map((c) => c.textContent),
    actions: li.children.filter((c) => c.tagName === "FORM").map((f) => f.children.find((i) => i.getAttribute("name") === "action").getAttribute("value")),
  }));
}

test("a running job fills the strip and shows it", () => {
  const p = parts();
  p.root.hidden = true;
  assert.equal(apply({ job, line: [] }, p, null), job);
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
  apply({ job: { ...job, value: null }, line: [] }, p, null);
  assert.equal("value" in p.bar.attrs, false);
});

test("the waiting list shows every entry in order, with Discard for a job and Start now only for a held one", () => {
  const p = parts();
  const held = { ...waiting, ticket: "t2", title: "Add sectors again", held: "The job before it did not succeed (failed).",
    note: "The job before it did not succeed (failed).", startable: true };
  const run = { ticket: "t3", type: "run", title: "Regenerate sector 7", who: null, held: null, position: 1,
    note: "Waiting for the machine, place 1.", startable: false, discardable: false };
  apply({ job, line: [waiting, held, run] }, p, null);
  assert.equal(p.queued.hidden, false);
  assert.deepEqual(items(p), [
    { texts: ["Redo scatters", "Starts when the running job ends."], actions: ["discard_queued"] },
    { texts: ["Add sectors again", "The job before it did not succeed (failed)."], actions: ["start_queued", "discard_queued"] },
    { texts: ["Regenerate sector 7", "Waiting for the machine, place 1."], actions: [] },
  ]);
  const first = p.line.children[0].children.find((c) => c.tagName === "FORM");
  assert.equal(first.getAttribute("action"), "/admin/generate");
  const sent = Object.fromEntries(first.children.filter((c) => c.tagName === "INPUT").map((i) => [i.getAttribute("name"), i.getAttribute("value")]));
  assert.deepEqual(sent, { csrf_token: "tok", action: "discard_queued", ticket: "t1" });
});

test("an emptied list clears the entries", () => {
  const p = parts();
  apply({ job, line: [waiting] }, p, null);
  apply({ job, line: [] }, p, job);
  assert.equal(p.line.children.length, 0);
  assert.equal(p.queued.hidden, true);
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
  apply({ job: null, line: [{ ...waiting, held: "held", note: "held", startable: true }] }, p, null);
  assert.equal(p.root.hidden, false);
  assert.equal(p.run.hidden, true);
});
