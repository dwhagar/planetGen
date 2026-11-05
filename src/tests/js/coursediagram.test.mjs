// tests/js/coursediagram.test.mjs -- static/coursediagram.js, the link between the course diagram and the Galaxy Map's course (NAV.69).

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const D = await import(new URL("coursediagram.js", STATIC_URL).href);

// A stop of the diagram as the page has it: its place in the course, a class list and attributes.
function stop(index) {
  const classes = new Set();
  const attributes = {};
  return {
    dataset: { index: String(index) },
    classList: { toggle: (name, on) => (on ? classes.add(name) : classes.delete(name)), contains: (name) => classes.has(name) },
    setAttribute: (name, value) => { attributes[name] = value; },
    removeAttribute: (name) => { delete attributes[name]; },
    attributes,
    closest(selector) { return selector.includes("nav-diagram-stop") ? this : null; },
  };
}

function page(stops) {
  const listeners = {};
  return {
    stops,
    querySelectorAll: () => stops,
    addEventListener: (name, handler) => { listeners[name] = handler; },
    listeners,
  };
}

test("marking a stop marks it in every drawing and unmarks the others", () => {
  // Three drawings of one route: the same indexes three times over.
  const doc = page([stop(0), stop(1), stop(2), stop(0), stop(1), stop(2)]);
  D.markStop(doc, 1);
  assert.deepEqual(doc.stops.map((el) => el.classList.contains("is-selected")), [false, true, false, false, true, false]);
  assert.equal(doc.stops[1].attributes["aria-current"], "true");
  assert.equal(doc.stops[0].attributes["aria-current"], undefined);
  D.markStop(doc, 2);
  assert.deepEqual(doc.stops.map((el) => el.classList.contains("is-selected")), [false, false, true, false, false, true]);
  assert.equal(doc.stops[1].attributes["aria-current"], undefined);
});

test("only a plain left click selects a stop; a modified or middle click follows the link", () => {
  assert.equal(D.isPlainClick({ button: 0 }), true);
  for (const modifier of ["ctrlKey", "metaKey", "shiftKey", "altKey"]) {
    assert.equal(D.isPlainClick({ button: 0, [modifier]: true }), false, modifier);
  }
  assert.equal(D.isPlainClick({ button: 1 }), false);
});

test("a click on a stop selects it on the map and stays on the page; a click on the selected stop follows the link", () => {
  const doc = page([stop(0), stop(1)]);
  const selected = [];
  let current = null;
  const canvas = {
    galaxyCourseSelect: (index) => { current = index; selected.push(index); return true; },
    galaxyCourseSelected: () => current,
    addEventListener: () => {},
  };
  D.connectCourseDiagram(canvas, doc);
  let prevented = 0;
  const click = (index, extra = {}) => doc.listeners.click({ target: doc.stops[index], button: 0, preventDefault: () => { prevented += 1; }, ...extra });
  click(1);
  assert.deepEqual(selected, [1]);
  assert.equal(prevented, 1);
  click(1);                                   // the selected stop again: its page
  assert.deepEqual(selected, [1]);
  assert.equal(prevented, 1);
  click(0, { ctrlKey: true });                // a modified click: its page
  assert.deepEqual(selected, [1]);
  assert.equal(prevented, 1);
});

test("the map's event marks the stop in the diagram", () => {
  const doc = page([stop(0), stop(1)]);
  const handlers = {};
  const canvas = { addEventListener: (name, handler) => { handlers[name] = handler; } };
  D.connectCourseDiagram(canvas, doc);
  handlers.galaxycoursestop({ detail: { index: 0 } });
  assert.equal(doc.stops[0].classList.contains("is-selected"), true);
  assert.equal(doc.stops[1].classList.contains("is-selected"), false);
});
