// html/static/findplace.js
//
// MAP.336: the choosing half of the Galaxy Map's find-a-place box (the box itself is findplacebox.js): what a typed line can be,
// as groups of results each holding a target for the map's one `goTo`. A line can be coordinates ("x y z" in parsecs), a sector
// designation or shell/band/wedge, an object ID (the printed 96-bit form, whose first group names the sector), the name of a charted
// region, a name the server finds (`/galaxy/locate`: sectors, systems, planets, moons), or a place visited before. No DOM and no
// three.js: tests/js/findplace.test.mjs holds it to this description.

export const RECENT_LIMIT = 8;
export const RECENT_KEY = "planetgen.recentPlaces";
const REGION_LIMIT = 6;
const RECENT_SHOWN = 5;
const CHARTED_SHOWN = 5;

// A result group's order on the list and its heading.
export const GROUP_TITLES = {
  here: "Go to this place",
  recent: "Recent",
  charted: "Charted places",
  region: "Regions",
  name: "Names",
  home: "Places",
};

const BODY_KINDS = ["star", "planet", "moon", "belt", "comet"];
const OBJECT_ID = /^[0-9A-Fa-f]+(?:-[0-9A-Fa-f]+){1,2}$/;

// {sector: {shell, band, wedge}, system, object, id} for a printed object ID ("FE81000A2B-0000005-000", optionally with its kind,
// "moon:FE81000A2B-0000005-003"): its first group is the sector's designation, the first two with the system's "-000" name the
// system, and `object` is the body's reference when the kind says it is a body and the last group is not zero (else null). null
// when `text` is not an object ID. `parseSector(designation)` is galaxystages' parseSectorDesignation.
export function objectIdOf(text, parseSector) {
  const raw = String(text == null ? "" : text).trim();
  const kind = (/^([a-z_]+):/i.exec(raw) || [])[1];
  const bare = raw.replace(/^[a-z_]+:/i, "");
  if (!OBJECT_ID.test(bare)) return null;
  const groups = bare.toUpperCase().split("-");
  const sector = parseSector(groups[0]);
  if (!sector) return null;
  const last = groups[2];
  const isBody = !!kind && BODY_KINDS.indexOf(kind.toLowerCase()) >= 0 && last !== undefined && /[^0]/.test(last);
  return {
    sector: sector, system: groups[0] + "-" + groups[1] + "-" + (last !== undefined ? "0".repeat(last.length) : "000"),
    object: isBody ? kind.toLowerCase() + ":" + groups.join("-") : null, id: groups.join("-"),
  };
}

function lower(text) {
  return String(text == null ? "" : text).trim().toLowerCase();
}

// The charted regions whose name contains `text` (the ones named for it first), as results.
export function regionResults(blobs, text) {
  const wanted = lower(text);
  if (!wanted) return [];
  const found = [];
  (blobs || []).forEach(function (blob) {
    const name = lower(blob.name || "Region " + blob.key);
    const at = name.indexOf(wanted);
    if (at < 0) return;
    found.push({ at: at, name: name, blob: blob });
  });
  found.sort(function (a, b) { return a.at - b.at || (a.name < b.name ? -1 : a.name > b.name ? 1 : 0); });
  return found.slice(0, REGION_LIMIT).map(function (f) {
    const label = f.blob.name || "Region " + f.blob.key;
    return {
      kind: "region", label: label, detail: f.blob.parent ? "region in " + f.blob.parent : "region",
      target: { region: { block: f.blob.block, center: [f.blob.x, f.blob.y, f.blob.z], radius: f.blob.radiusPc }, label: label },
    };
  });
}

// MAP.357: the charted regions with the most sectors generated, biggest first, for a blank line: the places worth going to.
export function chartedResults(blobs) {
  const list = (blobs || []).filter(function (blob) { return blob.sectorCount > 0; }).slice();
  list.sort(function (a, b) { return (b.sectorCount - a.sectorCount) || (a.key < b.key ? -1 : a.key > b.key ? 1 : 0); });
  return list.slice(0, CHARTED_SHOWN).map(function (blob) {
    const label = blob.name || "Region " + blob.key;
    return {
      kind: "region", label: label, detail: blob.sectorCount + " sectors generated",
      target: { region: { block: blob.block, center: [blob.x, blob.y, blob.z], radius: blob.radiusPc }, label: label },
    };
  });
}

// The places visited before whose label contains `text`; all of them, newest first, for a blank line.
export function recentResults(recents, text) {
  const wanted = lower(text);
  return (recents || []).filter(function (r) { return !wanted || lower(r.label).indexOf(wanted) >= 0; })
    .slice(0, RECENT_SHOWN).map(function (r) {
      return { kind: "recent", label: r.label, detail: "visited", target: Object.assign({}, r.target, { label: r.label }) };
    });
}

// The visit log: the last RECENT_LIMIT places gone to by name, newest first, in `storage` (localStorage, which can be missing or
// throw, so every access is guarded and the log is then only for this page).
export function createRecents(storage) {
  let memory = [];
  function read() {
    try {
      const parsed = JSON.parse(storage.getItem(RECENT_KEY) || "[]");
      return Array.isArray(parsed) ? parsed.filter(function (r) { return r && typeof r.label === "string" && r.target; }) : [];
    } catch (error) {
      return memory;
    }
  }
  return {
    list: read,
    // Records that `label` was gone to, `target` being what goTo takes to go there again (plain data only).
    add: function (label, target) {
      const text = String(label || "").trim();
      if (!text || !target) return;
      memory = [{ label: text, target: target }].concat(read().filter(function (r) { return lower(r.label) !== lower(text); })).slice(0, RECENT_LIMIT);
      try { storage.setItem(RECENT_KEY, JSON.stringify(memory)); } catch (error) { /* this page only */ }
    },
    clear: function () {
      memory = [];
      try { storage.removeItem(RECENT_KEY); } catch (error) { /* this page only */ }
    },
  };
}

// The printed ID of the system an object's printed ID is in or is: its first two groups, then the system's own "000".
function systemIdOf(id) {
  const groups = String(id).toUpperCase().split("-");
  return groups.length < 2 ? null : groups[0] + "-" + groups[1] + "-" + "0".repeat(groups[2] ? groups[2].length : 3);
}

function isBodyKind(kind) {
  return kind === "planet" || kind === "moon";
}

// What a server name match (`/galaxy/locate`) goes to: a sector is gone to; a system opens its sector in place and its own scene;
// a planet or moon opens its system with the body picked, as an object ID does.
export function nameResult(match) {
  const sector = { shell: match.shell_index, band: match.band_index, wedge: match.wedge_index };
  const label = match.name || "(unnamed)";
  const chain = (match.parents || []).filter(function (p) { return p.kind === "planet" || p.kind === "system"; })
    .reverse().map(function (p) { return p.name; }).filter(Boolean).join(", ");
  return {
    kind: match.kind, label: label, detail: match.kind + (chain ? " of " + chain : ""),
    target: match.kind === "sector" ? { sector: sector, label: label } : {
      sector: sector, open: true, system: systemIdOf(match.id), object: isBodyKind(match.kind) ? match.kind + ":" + String(match.id).toUpperCase() : null, label: label,
    },
    match: match,
  };
}

// The groups for a typed line, in the order shown: [{key, title, items: [{kind, label, detail, target}]}], empty groups left out.
// `parsed` is galaxystages' parseAddress of the line ({sector, point}, {name} or {problem}); `idOf` is objectIdOf for the line (or
// null); `regions` the charted regions' blobs; `recents` the visit log; `names` the server's name matches (null while asked);
// `home` the result for the whole galaxy.
export function buildGroups(options) {
  const text = String(options.text == null ? "" : options.text).trim();
  const groups = [];
  const add = function (key, items) {
    if (items.length) groups.push({ key: key, title: GROUP_TITLES[key], items: items });
  };
  const parsed = options.parsed || {};
  const here = [];
  if (text && parsed.point) {
    here.push({
      kind: "point", label: parsed.point.map(function (v) { return Math.round(v * 100) / 100; }).join(", ") + " pc",
      detail: "coordinates", target: { center: parsed.point.slice(), radius: 1, label: "Point " + parsed.point.join(", ") },
    });
  } else if (text && options.id) {
    const id = options.id;
    here.push({
      kind: "object", label: id.id, detail: "object ID",
      target: { sector: id.sector, open: true, system: id.system, object: id.object, label: id.id },
    });
  } else if (text && parsed.sector) {
    here.push({ kind: "sector", label: "Sector " + sectorText(parsed.sector), detail: "sector", target: { sector: parsed.sector, label: "Sector " + sectorText(parsed.sector) } });
  }
  add("here", here);
  add("recent", recentResults(options.recents, text));
  if (text) add("region", regionResults(options.regions, text));
  else add("charted", chartedResults(options.regions));
  if (text && options.names) add("name", options.names.map(nameResult));
  const homeLabel = options.home ? options.home.label : "";
  if (options.home && (!text || lower(homeLabel).indexOf(lower(text)) >= 0 || /^(home|galaxy|whole galaxy)$/i.test(text))) {
    add("home", [{ kind: "home", label: homeLabel, detail: "the whole galaxy", target: options.home.target }]);
  }
  return groups;
}

function sectorText(sector) {
  return [sector.shell, sector.band, sector.wedge].join("/");
}

// Whether a typed line should be sent to the server as a name: it is not a problem line, not coordinates or a designation, and has at
// least two characters.
export function asksServer(parsed, text) {
  return !!(parsed && parsed.name && String(text).trim().length >= 2);
}
