// tests/js/viewmenusearch.test.mjs -- static/viewmenusearch.js (UX.279): the search over the View menu's entries.

import { test } from "node:test";
import assert from "node:assert/strict";

import { installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/galaxy");
const L = await import(new URL("maplayers.js", STATIC_URL).href);
const S = await import(new URL("viewmenusearch.js", STATIC_URL).href);

function entries() {
  const reg = L.createLayerRegistry("galaxy");
  const add = (id, label, heading, group, keywords) => reg.register({ id, label, heading, group, keywords });
  add("star", "Stars", "Show", "Show on the map", ["star", "sun"]);
  add("black-hole", "Black holes", "Show", "Show on the map", ["black hole", "singularity"]);
  add("nebula", "Nebulae", "Show", "Show on the map", ["nebula", "cloud"]);
  add("class-m", "Class M", "Show", "World classes", ["rocky", "red dwarf world"]);
  add("class-k", "Class K", "Show", "World classes", ["rocky", "orange"]);
  add("grid", "Cell grid", "Color and fill", null, ["grid", "cells", "lines"]);
  add("spin", "Spin", "Move and Settings", null, ["spin", "rotate", "turn"]);
  return reg.entries();
}

test("words are lower case runs of letters and digits", () => {
  assert.deepEqual(S.wordsOf("World classes, M-type"), ["world", "classes", "m", "type"]);
  assert.deepEqual(S.wordsOf("  "), []);
  assert.deepEqual(S.wordsOf(null), []);
  assert.deepEqual(S.wordsOf("Nébulæ"), ["nébulæ"]);
});

test("a blank search is not a search: every entry stays and nothing is forced open", () => {
  const all = entries();
  for (const blank of ["", "   ", "--", null, undefined]) {
    const found = S.searchEntries(all, blank);
    assert.equal(found.active, false);
    assert.equal(found.count, all.length);
    assert.deepEqual(found.headings, []);
    assert.deepEqual(found.groups, []);
    assert.equal(found.text, "");
  }
});

test("a word finds an entry by the start of its name, its group's name or a keyword, whatever the case", () => {
  const all = entries();
  assert.deepEqual(S.searchEntries(all, "BLA").ids, ["black-hole"]);
  assert.deepEqual(S.searchEntries(all, "singul").ids, ["black-hole"], "a keyword");
  assert.deepEqual(S.searchEntries(all, "world").ids, ["class-m", "class-k"], "the group name finds the family");
  assert.deepEqual(S.searchEntries(all, "rotate").ids, ["spin"]);
  assert.deepEqual(S.searchEntries(all, "ell").ids, [], "a word starts a word; it is not found in the middle of one");
});

test("every typed word must be found, in any order", () => {
  const all = entries();
  assert.deepEqual(S.searchEntries(all, "class rocky orange").ids, ["class-k"]);
  assert.deepEqual(S.searchEntries(all, "orange class").ids, ["class-k"]);
  assert.deepEqual(S.searchEntries(all, "class sun").ids, []);
});

test("the matches keep the menu's order, and the headings and groups to open are those with a match", () => {
  const all = entries();
  const found = S.searchEntries(all, "c");
  assert.deepEqual(found.ids, ["nebula", "class-m", "class-k", "grid"], "cloud, Class, Class, Cell");
  assert.deepEqual(found.headings, ["Show", "Color and fill"]);
  assert.deepEqual(found.groups, [{ heading: "Show", group: "Show on the map" }, { heading: "Show", group: "World classes" }],
    "a layer with no group opens only its heading");
});

test("the count is told in words for a screen reader", () => {
  const all = entries();
  assert.equal(S.searchEntries(all, "spin").text, "1 match");
  assert.equal(S.searchEntries(all, "rocky").text, "2 matches");
  assert.equal(S.searchEntries(all, "zzz").text, "No matches");
  assert.equal(S.searchEntries(all, "zzz").count, 0);
  assert.deepEqual(S.searchEntries(all, "zzz").headings, []);
});

test("a search does not change the entries it is given", () => {
  const all = entries();
  const before = JSON.stringify(all);
  S.searchEntries(all, "class");
  assert.equal(JSON.stringify(all), before);
});
