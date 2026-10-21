// tests/js/adminjobstrip.test.mjs -- static/adminjobstrip.js, the admin job strip (ADM.119): how the status JSON's
// `strip` fills the strip, and what it shows when the job ends.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const { apply } = await import(new URL("adminjobstrip.js", STATIC_URL).href);

function parts() {
  const node = () => ({ hidden: false, textContent: "", value: 0, attrs: {},
    setAttribute(name, value) { this.attrs[name] = value; }, removeAttribute(name) { delete this.attrs[name]; } });
  return { root: node(), title: node(), text: node(), bar: node(), elapsed: node(), cancel: node(), jobId: node() };
}

const strip = { id: "j1", title: "Add sectors", url: "/admin/generate/jobs/j1", text: "Stage 2 of 3", elapsed_text: "3 m 05 s", value: 0.4 };

test("a running job fills the strip and shows it", () => {
  const p = parts();
  p.root.hidden = true;
  assert.equal(apply(strip, p, null), strip);
  assert.equal(p.root.hidden, false);
  assert.equal(p.title.textContent, "Add sectors");
  assert.equal(p.title.attrs.href, "/admin/generate/jobs/j1");
  assert.equal(p.jobId.value, "j1");
  assert.equal(p.bar.value, 0.4);
  assert.equal(p.elapsed.textContent, "3 m 05 s");
});

test("a job that does not know how far it is gets an unmoving bar", () => {
  const p = parts();
  p.bar.attrs.value = "0.2";
  apply({ ...strip, value: null }, p, null);
  assert.equal("value" in p.bar.attrs, false);
});

test("nothing running hides the strip, but a job that just ended says it finished", () => {
  const p = parts();
  assert.equal(apply(null, p, null), null);
  assert.equal(p.root.hidden, true);
  const q = parts();
  assert.equal(apply(null, q, strip), null);
  assert.equal(q.root.hidden, false);
  assert.equal(q.text.textContent, "Finished.");
  assert.equal(q.cancel.hidden, true);
});
