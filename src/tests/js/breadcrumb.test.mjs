// tests/js/breadcrumb.test.mjs -- static/breadcrumb.js, the one breadcrumb
// every page and map shares (NAV.14): links and buttons for the steps, the
// current one as text, and the steps a server-rendered list holds.

import { test } from "node:test";
import assert from "node:assert/strict";

import { h, installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/system/1");
const { createBreadcrumb, stepsOfList } = await import(new URL("breadcrumb.js", STATIC_URL).href);

test("set draws a link or a button per step and the current one as text", () => {
  installDom("http://localhost/system/1");
  const nav = h("nav", {});
  const picked = [];
  const crumbs = createBreadcrumb(nav, { onSelect: (step) => picked.push(step.label) });
  crumbs.set([
    { label: "Home", href: "/" },
    { label: "Galaxy" },
    { label: "Sector 5", last: true },
  ]);
  const items = Array.from(nav.querySelector("ol").children);
  assert.equal(items.length, 3);
  assert.equal(nav.querySelector("a").getAttribute("href"), "/");
  nav.querySelector("button").click();
  assert.deepEqual(picked, ["Galaxy"]);
  const here = nav.querySelector("[aria-current]");
  assert.equal(here.textContent, "Sector 5");
});

test("stepsOfList reads a server-rendered list: links, then the current page", () => {
  installDom("http://localhost/system/1");
  const nav = h("nav", {}, [h("ol", {}, [
    h("li", {}, [h("a", { href: "/" }, "Home")]),
    h("li", {}, [h("a", { href: "/sectors" }, "Sectors")]),
    h("li", {}, [h("span", { "aria-current": "page" }, "Voranthis")]),
  ])]);
  const steps = stepsOfList(nav);
  assert.deepEqual(steps.map((s) => [s.label, s.href, s.last]), [
    ["Home", "/", false], ["Sectors", "/sectors", false], ["Voranthis", null, true],
  ]);
});
