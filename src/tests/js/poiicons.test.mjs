// tests/js/poiicons.test.mjs -- static/poiicons.js (MAP.320): the point-of-interest symbols the map and the pages share.

import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";

import { STATIC_URL } from "./fakedom.mjs";

const P = await import(new URL("poiicons.js", STATIC_URL).href);
const M = await import(new URL("galaxymarkers.js", STATIC_URL).href);

const fileText = (id) => readFileSync(new URL(`poi/${id}.svg`, STATIC_URL), "utf8");

test("seven symbols, each with a file, a name and a token, the two the score colours having none", () => {
  assert.deepEqual(P.POI_ICONS.map((icon) => icon.id), ["habitable-world", "uninhabitable-world", "inhabited-world",
    "gas-giant", "asteroid-field", "sector", "region"]);
  assert.deepEqual(P.POI_ICONS.filter((icon) => icon.token === null).map((icon) => icon.id), ["habitable-world", "inhabited-world"]);
  assert.equal(new Set(P.POI_ICONS.map((icon) => icon.label)).size, 7);
  assert.equal(P.poiIcon("region").label, "Region");
  assert.equal(P.poiIcon("star"), null);
});

test("a symbol's token is its own, or the habitability tier's for a world the score colours", () => {
  assert.equal(P.tokenOf("uninhabitable-world", 0), "--poi-brown");     // brown whatever the tier: never red
  assert.equal(P.tokenOf("gas-giant", null), "--poi-gas");
  assert.equal(P.tokenOf("inhabited-world", 0), "--tier-blue");
  assert.equal(P.tokenOf("inhabited-world", 3), "--tier-red");
  assert.equal(P.tokenOf("habitable-world", null), "--tier-green");
  assert.throws(() => P.tokenOf("habitable-world", 4));
  assert.throws(() => P.tokenOf("habitable-world", 1.5));
  assert.throws(() => P.tokenOf("star", 0));
});

test("every file reads into square-viewBox shapes for a canvas", () => {
  for (const icon of P.POI_ICONS) {
    const shapes = P.shapesOf(fileText(icon.id));
    assert.equal(shapes.size, 24, icon.id);
    assert.ok(shapes.paths.length >= 1 && shapes.paths.every((path) => /^M[\d. -]/.test(path.d)), icon.id);
  }
  // A shape with a hole is drawn even-odd, as the file says.
  assert.ok(P.shapesOf(fileText("habitable-world")).paths[0].evenodd);
  assert.equal(P.shapesOf(fileText("asteroid-field")).paths.length, 4);
  assert.ok(P.shapesOf(fileText("asteroid-field")).paths.every((path) => !path.evenodd));
});

test("a file that is not a square icon is refused", () => {
  assert.throws(() => P.shapesOf('<svg viewBox="0 0 24 12"><path d="M0 0z"/></svg>'));
  assert.throws(() => P.shapesOf('<svg viewBox="0 0 24 24"></svg>'));
  assert.throws(() => P.shapesOf("<svg></svg>"));
});

test("a marker kind can be a symbol: the habitable worlds are, and the legend line carries the symbol and its token", () => {
  assert.equal(M.HABITABLE_KIND.icon, "habitable-world");
  assert.equal(M.HABITABLE_KIND.glyph, undefined);
  assert.deepEqual(M.legendOf(M.HABITABLE_KIND), { text: "Habitable worlds", icon: "habitable-world", token: "--tier-green" });
  assert.deepEqual(M.legendOf(M.BLACK_HOLE_KIND), { text: "Black holes", swatch: "#c58cff", swatchClass: "map-marker-swatch" });
  assert.equal(M.tokenOfItem(M.HABITABLE_KIND, { tier: 0 }), "--tier-blue");
  assert.equal(M.tokenOfItem(M.HABITABLE_KIND, {}), "--tier-green");
  assert.equal(M.tokenOfItem(M.BLACK_HOLE_KIND, {}), null);       // a glyph kind has a color, not a token
  const registry = M.createMarkerKinds();
  registry.register(M.HABITABLE_KIND);
  assert.throws(() => registry.register({ id: "x", label: "X", icon: "star" }));
  assert.throws(() => registry.register({ id: "y", label: "Y", glyph: "nope", icon: "region" }));    // an unknown glyph is not saved by an icon
});
