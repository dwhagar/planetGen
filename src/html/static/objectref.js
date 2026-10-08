// html/static/objectref.js
//
// The browser mirror of planetgen/galaxy/objectref.py: one reference form,
// `<kind>:<id>`, for every object (NAV.7). A bare number means a system.
// Keep the kinds in step with the Python module; tests/test_objectref.py
// checks both.

export const GALAXY = "galaxy";

export const KINDS = [
  "sector", "system", "star", "planet", "moon", "belt", "comet",
  "nebula", "asteroid_field", "black_hole", "neutron_star",
  "supernova_remnant", "rogue_planet", "interstellar_comet", "quasar",
];

const REF = /^(?:([a-z_]+):)?(\d+)$/;

// The reference for an object: format("moon", 7) is "moon:7".
export function format(kind, id) {
  if (!KINDS.includes(kind)) throw new Error(`unknown object kind: ${kind}`);
  const n = Number(id);
  if (!Number.isInteger(n) || n < 0) throw new Error(`bad object id: ${id}`);
  return `${kind}:${n}`;
}

// {kind, id} for a reference, or null when it is not one.
export function parse(raw) {
  const match = REF.exec(String(raw).trim());
  if (!match || (match[1] !== undefined && !KINDS.includes(match[1]))) return null;
  return { kind: match[1] || "system", id: Number(match[2]) };
}
