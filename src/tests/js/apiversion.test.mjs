// tests/js/apiversion.test.mjs -- static/apiversion.js (API.31): a request to the site's own /api/ carries the API
// version header, whatever else it sends; a request to anything else does not.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const { install, isApiAddress } = await import(new URL("apiversion.js", STATIC_URL).href);

const BASE = "https://planets.example/galaxy/3";
const HEADER = "X-PlanetGen-API-Version";

function target() {
  const calls = [];
  return { calls, fetch: async (input, init) => { calls.push({ input, init }); return "answer"; } };
}

test("the site's own /api/ is recognised, anything else is not", () => {
  assert.equal(isApiAddress("/api/sectors", BASE), true);
  assert.equal(isApiAddress("api/sectors", BASE), false);
  assert.equal(isApiAddress("https://planets.example/api/health?db=x", BASE), true);
  assert.equal(isApiAddress("https://elsewhere.example/api/health", BASE), false);
  assert.equal(isApiAddress("/galaxy", BASE), false);
  assert.equal(isApiAddress(new URL("/api/x", BASE), BASE), true);
});

test("a fetch of /api/ gains the version header and keeps its own", async () => {
  const t = target();
  install(t, BASE, HEADER, "6");
  assert.equal(await t.fetch("/api/sectors", { cache: "no-store", headers: { Accept: "application/json" } }), "answer");
  const sent = t.calls[0];
  assert.equal(sent.init.cache, "no-store");
  assert.equal(sent.init.headers.get(HEADER), "6");
  assert.equal(sent.init.headers.get("Accept"), "application/json");
});

test("a fetch with no options, and one given a Request, also carries it", async () => {
  const t = target();
  install(t, BASE, HEADER, "6");
  await t.fetch("/api/health");
  assert.equal(t.calls[0].init.headers.get(HEADER), "6");
  await t.fetch(new Request("https://planets.example/api/health", { headers: { "X-Mine": "1" } }));
  assert.equal(t.calls[1].init.headers.get(HEADER), "6");
  assert.equal(t.calls[1].init.headers.get("X-Mine"), "1");
});

test("a fetch of a page or another site is left as it was", async () => {
  const t = target();
  install(t, BASE, HEADER, "6");
  const options = { cache: "no-store" };
  await t.fetch("/galaxy/4", options);
  await t.fetch("https://elsewhere.example/api/health");
  assert.equal(t.calls[0].init, options);
  assert.equal(t.calls[1].init, undefined);
});

test("an XMLHttpRequest to /api/ sets the header too", () => {
  const sent = [];
  class Fake {
    open() { this.opened = true; }
    setRequestHeader(name, value) { sent.push([name, value]); }
  }
  const t = { XMLHttpRequest: Fake };
  install(t, BASE, HEADER, "6");
  new t.XMLHttpRequest().open("GET", "/api/sectors");
  new t.XMLHttpRequest().open("GET", "/galaxy");
  assert.deepEqual(sent, [[HEADER, "6"]]);
});

test("with no version nothing is wrapped", () => {
  const t = target();
  const original = t.fetch;
  install(t, BASE, HEADER, "");
  assert.equal(t.fetch, original);
});
