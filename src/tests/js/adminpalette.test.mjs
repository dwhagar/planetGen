// tests/js/adminpalette.test.mjs -- static/adminpalette.js, the admin palette and keyboard set (ADM.120, ADM.122):
// the ranking (the same as web/adminregistry.py `search`) and what each key means.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const { rank, meaning, typing, GO_KEYS } = await import(new URL("adminpalette.js", STATIC_URL).href);

const items = [
  { kind: "page", label: "Settings", summary: "Every site setting.", keywords: ["config"] },
  { kind: "setting", label: "Page answer cache", summary: "How long a page is kept.", keywords: ["cache"] },
  { kind: "page", label: "Jobs and queue", summary: "Background work.", keywords: ["cache", "warm"] },
  { kind: "task", label: "Add sectors", summary: "Make sectors.", keywords: [] },
];

test("a word that starts the label ranks above one inside it, which ranks above the summary", () => {
  assert.deepEqual(rank(items, "cache").map((i) => i.label), ["Page answer cache", "Jobs and queue"]);
  assert.deepEqual(rank(items, "set").map((i) => i.label)[0], "Settings");
});

test("every word must match, and an empty query lists the pages", () => {
  assert.deepEqual(rank(items, "cache warm").map((i) => i.label), ["Jobs and queue"]);
  assert.deepEqual(rank(items, "zzz"), []);
  assert.deepEqual(rank(items, "").map((i) => i.kind), ["page", "page"]);
  assert.equal(rank(items, "a", 2).length, 2);
});

const press = (key, extra = {}) => ({ key, target: { tagName: "BODY" }, ...extra });

test("keys: Ctrl+K and / open, ? helps, Esc closes, Ctrl+Enter reviews", () => {
  assert.equal(meaning(press("k", { ctrlKey: true }), false), "open");
  assert.equal(meaning(press("/"), false), "open");
  assert.equal(meaning(press("?"), false), "help");
  assert.equal(meaning(press("Escape"), false), "close");
  assert.equal(meaning(press("Enter", { ctrlKey: true }), false), "review");
});

test("g then g, s, j or d goes to a destination, and only right after g", () => {
  for (const [key, destination] of Object.entries(GO_KEYS)) {
    assert.equal(meaning(press(key), true), `go:${destination}`);
    assert.equal(meaning(press(key), false), "");
  }
});

test("single keys are off while typing, but Ctrl+K and Esc still work", () => {
  const field = { tagName: "INPUT" };
  assert.equal(typing(field), true);
  assert.equal(typing({ tagName: "DIV", isContentEditable: true }), true);
  assert.equal(meaning({ key: "/", target: field }, false), "");
  assert.equal(meaning({ key: "s", target: field }, true), "");
  assert.equal(meaning({ key: "k", ctrlKey: true, target: field }, false), "open");
  assert.equal(meaning({ key: "Escape", target: field }, false), "close");
});
