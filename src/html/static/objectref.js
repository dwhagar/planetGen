// html/static/objectref.js
//
// The browser mirror of planetgen/galaxy/objectref.py's public form
// (`parse_public`, `format_public`): one reference form, `<kind>:<id>`, for
// every object (NAV.7), where the id is the object's printed ID
// (`FE81000A2B-0000005-000`, API.23). A bare three-part ID means a system.
// Ids stay strings: a sector's printed ID can be all decimal digits and is
// never a number. Keep the kinds in step with the Python module;
// tests/test_objectref.py checks both.

export const GALAXY = "galaxy";

export const KINDS = [
  "sector", "system", "star", "planet", "moon", "belt", "comet",
  "nebula", "asteroid_field", "black_hole", "neutron_star",
  "supernova_remnant", "rogue_planet", "interstellar_comet", "quasar",
];

const ID = /^[0-9A-Fa-f]+(?:-[0-9A-Fa-f]+){0,2}$/;
const REF = /^(?:([a-z_]+):)?([0-9A-Fa-f]+(?:-[0-9A-Fa-f]+){0,2})$/;

// The reference for an object: format("moon", "FE81000A2B-0000005-003") is
// "moon:FE81000A2B-0000005-003".
export function format(kind, id) {
  if (!KINDS.includes(kind)) throw new Error(`unknown object kind: ${kind}`);
  const text = String(id);
  if (!ID.test(text)) throw new Error(`bad object id: ${id}`);
  return `${kind}:${text.toUpperCase()}`;
}

// {kind, id} for a reference (id as upper-case text), or null when it is not one.
export function parse(raw) {
  const match = REF.exec(String(raw).trim());
  if (!match || (match[1] !== undefined && !KINDS.includes(match[1]))) return null;
  if (match[1] === undefined && !match[2].includes("-")) return null;
  return { kind: match[1] || "system", id: match[2].toUpperCase() };
}
