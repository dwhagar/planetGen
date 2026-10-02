// html/static/galaxymap3d.js
// (Why three.js and not another renderer: docs/html-interface.md.)
//
// Renders the interactive 3D Galaxy Map (lib/galaxymap3d.py) as a real
// WebGL scene (three.js, vendored at ./vendor/three.module.min.js -- see
// sectormap.js's own docstring for why it's vendored rather than loaded
// from a CDN). Unlike sectormap.js's scene (every star/cloud baked into
// one JSON block at page load, camera always orbiting a fixed origin),
// this camera's own orbit TARGET moves freely through the galaxy -- most
// of what's drawn is fetched live from /galaxy/tiles (the JSON's
// fetchPath) as the camera moves, one fixed cube of space ("tile") at a time,
// never baked into the page beyond the very first frame
// (#galaxymap3d-data's own "initial" payload) -- see the "Cube tiles"
// section below.
//
// Coordinate convention: every position here is the design doc's own
// galaxy-frame parsecs (docs/design/galaxy-coordinate-system.md), passed
// straight through as three.js world units with NO axis flip/remap --
// x/y span the galactic plane, +z is galactic north. `camera.up` is set
// to (0, 0, 1) once, at scene setup, specifically so that convention
// reads as "up" on screen without needing any per-position sign flip the
// way this map's own former flat SVG projection (a legacy "+y is down on
// screen" convention) needed one.
//
// The galaxy's sector grid is drawn as blocks of whole sectors, one
// drill-down stage at a time (./galaxystageview.js, rules in
// ./galaxystages.js): the visitor picks a quarter, a layer, an arc, a
// layer, an arc, ... down to a sector. The whole galaxy and its quarters
// are seen from straight above; below them the view can be turned, moved
// and zoomed. Density is computed right here from the galaxy's own
// analytic shape, not fetched, and colors each block.
// Unfilled space is mostly see-through; a block holding generated
// sectors is amber and grows more solid with their share (MAP.37).
//
// Nebulae and supernova remnants are drawn over the blocks as soft
// translucent spheres their real size (each tile lists the ones reaching
// into it); clicking one shows it and links to its page -- see "Clouds".
//
// Stars are points of light a few pixels across with a soft glow, the
// same size at every zoom: the pre-placed bright stars (500 L☉ or more)
// everywhere, and the stars of generated systems fainter and fainter as
// the view closes in (MAP.51); clicking one shows it -- see "Stars".

// Sibling modules are imported with this module's own `?v=<version>`
// query (html/lib/fmt.py's `static_url`), so they are cached and
// refreshed with the page's script. A plain static `import "./x.js"`
// would drop the query: an update could then leave a stale copy cached,
// and a page that also loaded the same file by its versioned URL would
// get a second, separate instance of it.
const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const {
  blockSectorCount, blockSectorRanges, blockSlotRange, cellCoordinates, cellVertices, drillBlockBounds, drillSlotRange,
  sectorAddressAt, wedgeLines,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);
const { createStageView } = await import(`./galaxystageview.js${VERSION_QUERY}`);
const { createBlockScene } = await import(`./galaxyblocks.js${VERSION_QUERY}`);
const { formatDistancePc, LIGHTYEAR_M, PARSEC_M } = await import(`./distance.js${VERSION_QUERY}`);
const { formatNumber } = await import(`./numberformat.js${VERSION_QUERY}`);
const { boostLight, starLightBoost } = await import(`./starlight.js${VERSION_QUERY}`);
const { blockGenerateButtons, generateButtons } = await import(`./generatebuttons.js${VERSION_QUERY}`);

var canvas = document.getElementById("galaxymap3d-canvas");
var dataEl = document.getElementById("galaxymap3d-data");

function readSceneData() {
  if (!dataEl) {
    return null;
  }
  try {
    return JSON.parse(dataEl.textContent);
  } catch (err) {
    return null;
  }
}

var sceneData = readSceneData();

function cssVar(name, fallback) {
  var value = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return value || fallback;
}

// Whether the map's background (--bg-subtle, the light or dark theme) is
// light, from its computed color.
function isLightBackground() {
  var probe = document.createElement("span");
  probe.style.color = cssVar("--bg-subtle", "#000");
  document.body.appendChild(probe);
  var rgb = (getComputedStyle(probe).color.match(/[\d.]+/g) || [0, 0, 0]).map(Number);
  probe.remove();
  return 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2] > 140;
}

function addField(dl, label, value) {
  if (!value && value !== 0) {
    return;
  }
  var dt = document.createElement("dt");
  dt.textContent = label;
  var dd = document.createElement("dd");
  dd.textContent = value;
  dl.appendChild(dt);
  dl.appendChild(dd);
}

// A plain <a href> link (GET, bookmarkable, opens in a new tab).
function pageLink(href, label) {
  var link = document.createElement("a");
  link.href = href;
  link.className = "btn";
  link.textContent = label;
  return link;
}

// The sector page's URL: the server's template (sceneData.sectorUrl,
// built with page_url so it follows the sector page wherever it lives)
// with the id filled in.
function sectorUrl(id) {
  return String(sceneData.sectorUrl || "").replace("{id}", encodeURIComponent(id));
}

export function formatAddress(ringIndex, layerIndex, slotIndex) {
  return "ring " + ringIndex + " layer " + layerIndex + " slot " + slotIndex;
}

// The map's control buttons (lib/galaxymap3d.py's panel), by their
// data-action: what each does, given the map's own parts in `ctx`
// ({stageView, setTerritories(on, button), territoriesWanted(),
// wedgeGroup}). A button whose action isn't here does nothing.
export function mapControlHandlers(ctx) {
  return {
    "territories": function (button) {
      ctx.setTerritories(button.getAttribute("aria-pressed") !== "true", button);
      button.setAttribute("aria-pressed", String(ctx.territoriesWanted()));
    },
    "generated-only": function (button) {
      var on = button.getAttribute("aria-pressed") !== "true";
      button.setAttribute("aria-pressed", String(on));
      ctx.stageView.setGeneratedOnly(on);
    },
    "back": function () { ctx.stageView.travel(-1); },
    "forward": function () { ctx.stageView.travel(1); },
    "up": function () { ctx.stageView.up(); },
    "reset": function () { ctx.stageView.home(); },
    "reset-view": function () { ctx.stageView.resetView(); },
    "wedges": function (button) {
      ctx.wedgeGroup.visible = !ctx.wedgeGroup.visible;
      button.setAttribute("aria-pressed", String(ctx.wedgeGroup.visible));
    },
  };
}

// Sends each [data-action] button's clicks to its handler.
export function wireMapControls(controlsEl, handlers) {
  controlsEl.querySelectorAll("[data-action]").forEach(function (button) {
    button.addEventListener("click", function () {
      var handler = handlers[button.dataset.action];
      if (handler) {
        handler(button);
      }
    });
  });
}

// --- Info panel ----------------------------------------------------------

function showPlacedInfo(entry) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";

  var heading = document.createElement("h3");
  heading.textContent = entry.name || "Unnamed sector";
  panel.appendChild(heading);

  var relativeDensity = placedRelativeDensity(entry, sceneData.referenceDensityPerLy3);

  var dl = document.createElement("dl");
  addField(dl, "Systems", entry.system_count != null ? entry.system_count : 0);
  addField(dl, "Density", relativeDensity != null ? relativeDensity.toFixed(2) + "× local average" : null);
  addField(dl, "Distance from core", entry.galactic_radius_pc != null ? Math.round(entry.galactic_radius_pc) + " pc" : null);
  addField(dl, "Address", entry.ring_index != null ? formatAddress(entry.ring_index, entry.layer_index, entry.ring_slot_index) : null);
  addField(dl, "Designation", entry.designation);
  panel.appendChild(dl);
  if (sceneData.sectorUrl && entry.id != null) {
    panel.appendChild(pageLink(sectorUrl(entry.id), "View sector →"));
  }
}

// A share as a percentage, never rounding a nonzero one down to 0%.
function formatShare(share) {
  if (share > 0 && share < 0.001) {
    return "< 0.1%";
  }
  return (100 * share).toFixed(share < 0.1 ? 1 : 0) + "%";
}

function formatDeg(rad) {
  return ((rad * 180) / Math.PI).toFixed(3) + "°";
}

// provisional_sector_designation's hex code: ring << 33 | (layer + 4096)
// << 20 | slot (BigInt, since it runs past 32 bits).
export function sectorDesignation(ring, layer, slot) {
  var packed = (BigInt(ring) << 33n) | (BigInt(layer + 4096) << 20n) | BigInt(slot);
  return packed.toString(16).toUpperCase();
}

// One block's info: a single sector (m == 1) or a block of sectors.
// `cell` is {m, bounds, address} for a sector or {m, bounds, ring, seg,
// slab, density, edgePc, shape} for a block, plus `filled` (generated
// sectors in it) when it came from the drawn blocks. A generated sector's
// own panel is showPlacedInfo.
function showCellInfo(cell) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";
  var b = cell.bounds;
  var single = cell.m === 1;

  var heading = document.createElement("h3");
  heading.textContent = single ? "Sector cell" : "Sector block (" + cell.m + " sectors a side)";
  panel.appendChild(heading);

  var dl = document.createElement("dl");
  if (single) {
    var a = cell.address;
    addField(dl, "Address", formatAddress(a.ring, a.layer, a.slot));
    addField(dl, "Designation", sectorDesignation(a.ring, a.layer, a.slot));
    if (cell.filled != null) {
      addField(dl, "Generated", cell.filled > 0 ? "Yes" : "Not yet");
    }
  } else {
    var ranges = blockSectorRanges(cell.ring, cell.slab, cell.m);
    addField(dl, "Rings", ranges.ringFirst + "–" + ranges.ringLast);
    addField(dl, "Layers", ranges.layerFirst + "–" + ranges.layerLast);
    // Slot numbers restart in every ring, so give the innermost and
    // outermost member rings' ranges.
    [ranges.ringFirst, ranges.ringLast].forEach(function (ring) {
      var slots = blockSlotRange(cell.ring, cell.seg, cell.m, ring);
      addField(dl, "Slots in ring " + ring, slots.first + "–" + slots.last);
    });
    var total = blockSectorCount(cell.ring, cell.seg, cell.slab, cell.m, cell.edgePc, cell.shape);
    addField(dl, "Sectors", formatNumber(total));
    if (cell.filled != null) {
      addField(dl, "Generated", formatNumber(cell.filled) + (cell.filled > 0 && total > 0 ? " (" + formatShare(cell.filled / total) + ")" : ""));
    }
  }
  if (cell.density != null) {
    addField(dl, "Predicted density", cell.density.toFixed(2) + "× local average");
  }
  var coords = cellCoordinates(b);
  var c = coords.cartesian;
  addField(dl, "Center x, y, z", c.map(function (v) { return v.toFixed(1); }).join(", ") + " pc");
  addField(dl, "Cylindrical R, θ, z", coords.cylindrical[0].toFixed(1) + " pc, " + formatDeg(coords.cylindrical[1]) + ", " + coords.cylindrical[2].toFixed(1) + " pc");
  addField(dl, "Spherical r, θ, φ", coords.spherical[0].toFixed(1) + " pc, " + formatDeg(coords.spherical[1]) + ", " + formatDeg(coords.spherical[2]));
  addField(dl, "Radial width", formatDistancePc(b.r1 - b.r0));
  addField(dl, "Height", formatDistancePc(b.z1 - b.z0));
  addField(dl, "Mean arc length", formatDistancePc(((b.r0 + b.r1) / 2) * (b.t1 - b.t0)));
  panel.appendChild(dl);

  var corners = document.createElement("details");
  var summary = document.createElement("summary");
  summary.textContent = "8 corners (x, y, z pc)";
  corners.appendChild(summary);
  var list = document.createElement("ol");
  list.start = 0;
  cellVertices(b).forEach(function (v) {
    var item = document.createElement("li");
    item.textContent = v.map(function (n) { return n.toFixed(2); }).join(", ");
    list.appendChild(item);
  });
  corners.appendChild(list);
  panel.appendChild(corners);

  // A sector that isn't generated yet: for a logged-in admin
  // (sceneData.generate, set by lib/galaxymap3d.py only then), the same
  // Generate buttons the Sector Map gives a neighbor.
  if (single && !(cell.filled > 0) && sceneData.generate) {
    panel.appendChild(generateButtons(sceneData.generate, cell.address.ring, cell.address.layer, cell.address.slot,
      sceneData.edgeLy));
  }
}

// A drill-down block's info (galaxystageview.js): its rings, layers and
// slots, how many sectors it can hold and how many are generated, where
// it is, and (`enter`) a button that flies into it; `info.generate`
// (stages 7-8, admins only) adds the block's Generate buttons.
// `info.hint` goes under it.
function showBlockInfo(info, edgePc) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";
  var block = info.block;
  var b = drillBlockBounds(block, edgePc);
  var heading = document.createElement("h3");
  heading.textContent = "Block " + block.ring + "·" + block.wedge + " (" + block.m + " sectors a side)";
  panel.appendChild(heading);
  var half = (block.m - 1) / 2;
  var ringFirst = block.ring * block.m;
  var ringLast = ringFirst + block.m - 1;
  var dl = document.createElement("dl");
  addField(dl, "Rings", ringFirst + "–" + ringLast);
  addField(dl, "Layers", (block.slab * block.m - half) + "–" + (block.slab * block.m + half));
  [ringFirst, ringLast].forEach(function (ring) {
    var slots = drillSlotRange(block, ring);
    addField(dl, "Slots in ring " + ring, slots.first + "–" + slots.last);
  });
  if (info.total != null) {
    addField(dl, "Sectors", formatNumber(info.total));
  }
  if (info.generated != null) {
    addField(dl, "Generated", formatNumber(info.generated)
      + (info.generated > 0 && info.total > 0 ? " (" + formatShare(info.generated / info.total) + ")" : ""));
  }
  var coords = cellCoordinates(b);
  addField(dl, "Center x, y, z", coords.cartesian.map(function (v) { return v.toFixed(1); }).join(", ") + " pc");
  addField(dl, "Distance from core", formatDistancePc((b.r0 + b.r1) / 2));
  panel.appendChild(dl);
  if (info.enter) {
    var enter = document.createElement("button");
    enter.type = "button";
    enter.className = "btn";
    enter.textContent = "Fly into this block →";
    enter.addEventListener("click", info.enter);
    panel.appendChild(enter);
  }
  // Stages 7-8 for an admin: generate the block, or the layer shown.
  if (info.generate && sceneData.generate) {
    panel.appendChild(blockGenerateButtons(sceneData.generate, info.generate));
  }
  if (info.hint) {
    showHint(info.hint, true);
  }
}

// A hint paragraph in the info panel: replacing what's there, or (keep)
// under it.
function showHint(text, keep) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  if (!keep) {
    panel.textContent = "";
  }
  var hint = document.createElement("p");
  hint.className = "hint";
  hint.textContent = text;
  panel.appendChild(hint);
}

// --- Clouds: nebulae and supernova remnants --------------------------------

// A nebula's color by its type, the same hues the Sector Map uses
// (lib/starmap.py's _NEBULA_TYPE_COLORS and _NEBULA_TYPE_ALPHA, as core
// and edge opacity). A dark nebula is a near-black silhouette.
var NEBULA_LOOKS = {
  diffuse: ["#e3a6c8", 0.47, 0.14],
  emission: ["#ff6f91", 0.69, 0.25],
  reflection: ["#6fa8ff", 0.63, 0.22],
  planetary: ["#5be8c9", 0.66, 0.24],
  dark: ["#1c1c24", 0.91, 0.56],
};
var DEFAULT_NEBULA_LOOK = ["#c9a8e0", 0.56, 0.22];

function capitalize(text) {
  text = String(text || "");
  return text.charAt(0).toUpperCase() + text.slice(1);
}

// "Emission Nebula", "Supernova Remnant (Shell)".
function cloudTypeLabel(cloud) {
  if (cloud.type === "nebula") {
    return cloud.descriptor ? capitalize(cloud.descriptor) + " Nebula" : "Nebula";
  }
  return "Supernova Remnant" + (cloud.descriptor ? " (" + capitalize(cloud.descriptor) + ")" : "");
}

// The phenomenon page's URL: the server's template (sceneData.phenomenonUrl)
// with the type and id filled in.
function phenomenonUrl(cloud) {
  return String(sceneData.phenomenonUrl || "")
    .replace("{type}", encodeURIComponent(cloud.type))
    .replace("{id}", encodeURIComponent(cloud.id));
}

function showCloudInfo(cloud) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";
  var heading = document.createElement("h3");
  heading.textContent = cloud.name || cloudTypeLabel(cloud);
  panel.appendChild(heading);
  var dl = document.createElement("dl");
  addField(dl, "Type", cloudTypeLabel(cloud));
  addField(dl, "Class", cloud.class);
  addField(dl, "Radius", formatDistancePc(cloud.radius_pc));
  addField(dl, "Center x, y, z", [cloud.x, cloud.y, cloud.z].map(function (v) { return v.toFixed(1); }).join(", ") + " pc");
  addField(dl, "Distance from core", formatDistancePc(Math.hypot(cloud.x, cloud.y, cloud.z)));
  panel.appendChild(dl);
  if (sceneData.phenomenonUrl) {
    panel.appendChild(pageLink(phenomenonUrl(cloud), "View phenomenon →"));
  }
}

// A star's color from its surface temperature: Tanner Helland's fit to
// blackbody colors, as [r, g, b] in 0..1 (red giants orange, O and B
// stars blue-white).
function starColor(temperatureK) {
  var t = Math.min(40000, Math.max(1000, temperatureK || 5800)) / 100;
  var r = t <= 66 ? 255 : 329.698727446 * Math.pow(t - 60, -0.1332047592);
  var g = t <= 66 ? 99.4708025861 * Math.log(t) - 161.1195681661 : 288.1221695283 * Math.pow(t - 60, -0.0755148492);
  var b = t >= 66 ? 255 : t <= 19 ? 0 : 138.5177312231 * Math.log(t - 10) - 305.0447927307;
  return [r, g, b].map(function (v) { return Math.min(255, Math.max(0, v)) / 255; });
}

// "9,300 L☉", "1.23 × 10⁴ L☉" (UX.20), "0.0031 L☉".
function formatLuminosity(sol) {
  if (sol < 100) {
    return Number(sol.toPrecision(2)).toLocaleString("en-US", { maximumSignificantDigits: 2 }) + " L☉";
  }
  return formatNumber(sol) + " L☉";
}

function systemUrl(id) {
  return String(sceneData.systemUrl || "").replace("{id}", encodeURIComponent(id));
}

// A star on the map: a pre-placed bright star
// (queryDb.galaxy_bright_stars_in_box) or a generated system's star
// (queryDb.galaxy_generated_stars_in_box, `generated` set): what it is,
// and its system's page once its sector is filled.
function showStarInfo(star) {
  var panel = document.getElementById("galaxymap3d-info");
  if (!panel) {
    return;
  }
  panel.textContent = "";
  var heading = document.createElement("h3");
  heading.textContent = star.generated ? star.name || "Star" : "Bright star";
  panel.appendChild(heading);
  var dl = document.createElement("dl");
  addField(dl, "Type", [star.star_type, star.yerkes_class].filter(Boolean).join(" "));
  addField(dl, "Luminosity", formatLuminosity(star.luminosity_sol));
  addField(dl, "Temperature", formatNumber(star.temperature_k) + " K");
  if (star.radius_sol != null) {
    addField(dl, "Radius", formatNumber(Number(star.radius_sol.toPrecision(2)), 6, 0) + " R☉");
  }
  addField(dl, "Distance from core", formatDistancePc(Math.hypot(star.x, star.y, star.z)));
  addField(dl, "Sector", sectorDesignation(star.ring_index, star.layer_index, star.ring_slot_index));
  addField(dl, "Address", formatAddress(star.ring_index, star.layer_index, star.ring_slot_index));
  addField(dl, "System", star.system_id != null ? null : "Not generated yet (its sector isn't filled)");
  panel.appendChild(dl);
  if (sceneData.systemUrl && star.system_id != null) {
    panel.appendChild(pageLink(systemUrl(star.system_id), "View system →"));
  }
}

function rgba(hex, alpha) {
  var n = parseInt(hex.slice(1), 16);
  return "rgba(" + (n >> 16) + "," + ((n >> 8) & 255) + "," + (n & 255) + "," + alpha + ")";
}

// A cloud's sprite texture, drawn so the sphere's edge is the texture's
// edge (the sprite is scaled to the cloud's diameter): a nebula glows
// from a denser core out to nothing; a remnant is a thin bright shell
// around a faint interior, as Cassiopeia A or the Veil look.
var cloudTextures = new Map();

function cloudTexture(cloud) {
  var key = cloud.type === "nebula" ? "nebula:" + (NEBULA_LOOKS[cloud.descriptor] ? cloud.descriptor : "") : "remnant";
  var texture = cloudTextures.get(key);
  if (texture) {
    return texture;
  }
  var size = 128;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = size;
  canvasEl.height = size;
  var ctx = canvasEl.getContext("2d");
  var gradient = ctx.createRadialGradient(size / 2, size / 2, 0, size / 2, size / 2, size / 2);
  if (cloud.type === "nebula") {
    var look = NEBULA_LOOKS[cloud.descriptor] || DEFAULT_NEBULA_LOOK;
    gradient.addColorStop(0, rgba(look[0], look[1]));
    gradient.addColorStop(0.35, rgba(look[0], look[1] * 0.8));
    gradient.addColorStop(0.75, rgba(look[0], look[2]));
    gradient.addColorStop(1, rgba(look[0], 0));
  } else {
    gradient.addColorStop(0, rgba("#ffb07a", 0.12));
    gradient.addColorStop(0.6, rgba("#ffb07a", 0.2));
    gradient.addColorStop(0.82, rgba("#ff8a5c", 0.7));
    gradient.addColorStop(0.9, rgba("#8fd6ff", 0.6));
    gradient.addColorStop(1, rgba("#8fd6ff", 0));
  }
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, size, size);
  texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  cloudTextures.set(key, texture);
  return texture;
}

// --- Sprite textures -------------------------------------------------------

// A text label for a sprite, with its width / height in userData.aspect.
function makeLabelTexture(text, color) {
  var height = 64;
  var canvasEl = document.createElement("canvas");
  var ctx = canvasEl.getContext("2d");
  var font = "600 " + Math.round(height * 0.7) + "px system-ui, sans-serif";
  ctx.font = font;
  canvasEl.width = Math.ceil(ctx.measureText(text).width + height * 0.3);
  canvasEl.height = height;
  ctx.font = font;
  ctx.fillStyle = color;
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillText(text, canvasEl.width / 2, height / 2);
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  texture.userData.aspect = canvasEl.width / height;
  return texture;
}

function makeRingTexture(color) {
  var size = 64;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = canvasEl.height = size;
  var ctx = canvasEl.getContext("2d");
  var r = size / 2;
  ctx.beginPath();
  ctx.arc(r, r, r - 4, 0, Math.PI * 2);
  ctx.lineWidth = 3;
  ctx.strokeStyle = color;
  ctx.stroke();
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

// --- Filled-sector colors ---------------------------------------------
//
// Unlike lib/starmap.py (size/color baked server-side, once), these run
// client-side for EVERY fetch (initial payload and every live re-fetch
// alike) -- see galaxymap3d.py's own module docstring for why: almost
// everything drawn here arrives through a live fetch that never passes
// through that Python module again, so a server-computed style would
// only ever apply to the first frame.

// A placed sector's own REAL stellar density (system_count / edge_ly^3),
// relative to physical_constants.LOCAL_STELLAR_DENSITY_LY3 (the real
// local-neighborhood average this whole generator already calibrates
// against -- see lib/galaxymap3d.py's own referenceDensityPerLy3
// comment) -- 1.0 means exactly average, >1 denser, <1 sparser. `null`
// when edge_ly isn't available (a sector placed before per-sector edge
// tracking existed) rather than a false 0, so a generated sector's block
// (galaxyblocks.js) can fall back to a neutral mid-tone instead of
// reading as "empty". Generated sectors are colored on a warm
// bronze-to-gold range from PLACED_LOW_DENSITY_COLOR (sparse) to
// PLACED_HIGH_DENSITY_COLOR (dense), apart from the density shading's
// cooler dim-to-accent one, so they stand out at a glance.
var PLACED_LOW_DENSITY_COLOR = "#4a3f2e";
var PLACED_HIGH_DENSITY_COLOR = "#fff6df";

function placedRelativeDensity(entry, referenceDensityPerLy3) {
  if (!entry.edge_ly || !referenceDensityPerLy3) {
    return null;
  }
  var densityPerLy3 = (entry.system_count || 0) / Math.pow(entry.edge_ly, 3);
  return densityPerLy3 / referenceDensityPerLy3;
}

// The selection ring, drawn this many screen pixels across.
var HIGHLIGHT_PX = 9;

var highlightTexture = null;

function ensureTextures(accentColor) {
  if (!highlightTexture) {
    highlightTexture = makeRingTexture(accentColor);
  }
}

// --- Scene setup -----------------------------------------------------------

function initGalaxyMap3d(canvasEl, data) {
  var viewport = canvasEl.closest(".starmap-viewport");
  var accentColor = cssVar("--accent", "#4f5fe8");
  ensureTextures(accentColor);

  var renderer = new THREE.WebGLRenderer({ canvas: canvasEl, antialias: true, alpha: true, logarithmicDepthBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.setClearColor(0x000000, 0);
  if (THREE.SRGBColorSpace) {
    renderer.outputColorSpace = THREE.SRGBColorSpace;
  }

  var scene = new THREE.Scene();

  var FOV_DEG = data.fovDeg || 50;
  var farPlane = Math.max(data.maxViewRadiusPc * 6, 1000);
  var camera = new THREE.PerspectiveCamera(FOV_DEG, 1, Math.max(data.minViewRadiusPc / 50, 0.001), farPlane);
  camera.up.set(0, 0, 1);

  var MIN_RADIUS = data.minViewRadiusPc;
  var MAX_RADIUS = data.maxViewRadiusPc;
  var GALAXY_RADIUS = data.galaxyRadiusPc || data.maxViewRadiusPc;
  // The galaxy's real edge (its outermost ring, unpadded): the wedge
  // lines stop here (MAP.43).
  var GALAXY_EDGE = data.galaxyEdgePc || GALAXY_RADIUS;

  var target = new THREE.Vector3(data.initialCenter[0], data.initialCenter[1], data.initialCenter[2]);
  var initialRadius = data.initialRadiusPc;

  // theta: azimuth from +x in the xy-plane; phi: polar angle from +z --
  // the same (r, theta, phi) convention docs/design/galaxy-coordinate-
  // system.md and stellarObjects.galaxyGeometry.sector_position_pc use,
  // deliberately not THREE.Spherical (which assumes a +y-up world).
  var orbit = { radius: initialRadius, theta: THREE.MathUtils.degToRad(-32), phi: THREE.MathUtils.degToRad(60) };

  function offsetFromOrbit(o) {
    var sinPhi = Math.sin(o.phi);
    return new THREE.Vector3(
      o.radius * sinPhi * Math.cos(o.theta),
      o.radius * sinPhi * Math.sin(o.theta),
      o.radius * Math.cos(o.phi)
    );
  }

  function applyCamera() {
    var offset = offsetFromOrbit(orbit);
    camera.position.set(target.x + offset.x, target.y + offset.y, target.z + offset.z);
    camera.lookAt(target);
    // Raycasts read matrixWorld, which otherwise only updates at the next
    // render -- a click right after a move would pick against the old view.
    camera.updateMatrixWorld();
  }
  applyCamera();

  // --- Wedge lines ---------------------------------------------------------
  //
  // The master lines of the sector grid (galaxyprisms.wedgeLines): lines in
  // the galactic plane out to the galaxy's edge (its outermost ring, not
  // past it: MAP.43), 3 from the core, then 6, 12, ... each starting where
  // its zone does. The coarsest are labelled with their
  // bearing (degrees counterclockwise from +X, the zero meridian), so a
  // view can be placed around the galaxy at a glance. Each zone's lines
  // show only while they are at least WEDGE_MIN_GAP_PX apart on screen
  // (updateWedgeLevels). Drawn over everything (no depth test) and never
  // picked. The Wedges button hides them.
  var wedgeColor = new THREE.Color(cssVar("--text", "#e6e8f0"));
  // Guide lines only: mostly see-through, just visible enough to follow
  // (Boss, 2026-10-01). The bearing labels stay fully readable.
  var WEDGE_LINE_OPACITY = 0.22;
  var wedgeGroup = new THREE.Group();
  wedgeGroup.renderOrder = 2;
  scene.add(wedgeGroup);
  var WEDGE_LABEL_PX = 16;
  // Zones with more master lines than this (15 degrees apart) go unlabelled.
  var WEDGE_LABEL_MAX_MASTERS = 24;
  var WEDGE_MIN_GAP_PX = 24;
  var wedgeLabels = [];
  // [{masters, r0, objects: [LineSegments, label sprites...]}], coarsest first.
  var wedgeLevels = [];
  (function buildWedgeLines() {
    var reach = GALAXY_EDGE;
    var byMasters = new Map();
    wedgeLines(data.edgePc || 1, GALAXY_EDGE).forEach(function (line) {
      // A zone starting at (or past) the edge has nothing to draw.
      if (line.r0 >= reach) {
        return;
      }
      if (!byMasters.has(line.masters)) {
        byMasters.set(line.masters, []);
      }
      byMasters.get(line.masters).push(line);
    });
    var material = new THREE.LineBasicMaterial({
      color: wedgeColor, transparent: true, opacity: WEDGE_LINE_OPACITY, depthTest: false, depthWrite: false,
    });
    byMasters.forEach(function (lines, masters) {
      var level = { masters: masters, r0: lines[0].r0, lines: lines, objects: [], segments: null };
      var points = new Float32Array(lines.length * 6);
      lines.forEach(function (line, n) {
        var cos = Math.cos(line.angleRad);
        var sin = Math.sin(line.angleRad);
        points.set([line.r0 * cos, line.r0 * sin, 0, reach * cos, reach * sin, 0], 6 * n);
        if (masters > WEDGE_LABEL_MAX_MASTERS) {
          return;
        }
        var label = new THREE.Sprite(new THREE.SpriteMaterial({
          map: makeLabelTexture(String(line.bearingDeg).padStart(3, "0"), "#" + wedgeColor.getHexString()),
          transparent: true, depthTest: false, depthWrite: false, sizeAttenuation: false,
        }));
        label.position.set(reach * 1.06 * cos, reach * 1.06 * sin, 0);
        label.renderOrder = 2;
        wedgeLabels.push(label);
        level.objects.push(label);
      });
      var geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(points, 3));
      var segments = new THREE.LineSegments(geometry, material);
      segments.renderOrder = 2;
      segments.frustumCulled = false;
      level.objects.push(segments);
      level.segments = segments;
      level.objects.forEach(function (object) { wedgeGroup.add(object); });
      wedgeLevels.push(level);
    });
  })();

  // The part of the galaxy the drill-down shows ({r0, r1, a0, a1, z0,
  // z1}: radii and heights in pc, bearings in radians, and optionally
  // `cells`, the bounds {r0, r1, t0, t1, z0, z1} of the blocks in view),
  // or null for the whole galaxy. Zoomed in, only that wedge shows (Boss,
  // 2026-10-01): the wedge lines, stars and clouds are kept to its blocks
  // (MAP.44); zoomed out to the galaxy they run to its edge (MAP.43).
  var wedgeClip = null;
  var WEDGE_CLIP_EPSILON = 1e-6;

  function setWedgeClip(clip) {
    wedgeClip = clip;
    markClippedStars();
    // Zoomed in, the lines lie on the floor of the layers in view, not
    // the galactic plane, so seen at a slant they run along the blocks.
    wedgeGroup.position.z = clip && clip.z0 != null ? clip.z0 : 0;
    updateWedgeLevels(pcPerPixelAtTarget());
    updateClouds();
  }

  // How far round from bearing a0 bearing `angle` is, 0 to a full turn.
  function bearingFrom(angle, a0) {
    return (((angle - a0) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  }

  // Whether a point is inside the wedge shown (always, for the whole
  // galaxy): inside one of its blocks, or its bounds when it lists none.
  function inWedgeClip(x, y, z) {
    var c = wedgeClip;
    if (!c) {
      return true;
    }
    var r = Math.hypot(x, y);
    var angle = Math.atan2(y, x);
    var inside = function (b, t0, t1) {
      return r >= b.r0 && r <= b.r1 && bearingFrom(angle, t0) <= t1 - t0
        && (b.z0 == null || z >= b.z0) && (b.z1 == null || z <= b.z1);
    };
    if (c.cells) {
      return c.cells.some(function (b) { return inside(b, b.t0, b.t1); });
    }
    return inside(c, c.a0, c.a1);
  }

  // Shows each zone's lines while neighbours in that zone are at least
  // WEDGE_MIN_GAP_PX apart on screen, measured where the view reaches
  // farthest out (the focus's radius plus half the screen's height), or at
  // the zone's start. The bearing labels only show over the whole galaxy.
  function updateWedgeLevels(pcPerPixel) {
    var focusR = Math.hypot(target.x, target.y);
    var reachR = Math.min(GALAXY_EDGE, focusR + 0.5 * (canvasEl.clientHeight || 1) * pcPerPixel);
    wedgeLevels.forEach(function (level) {
      var gapPx = (2 * Math.PI * Math.max(level.r0, reachR)) / level.masters / pcPerPixel;
      var show = level.masters === wedgeLevels[0].masters || gapPx >= WEDGE_MIN_GAP_PX;
      level.objects.forEach(function (object) {
        object.visible = show && (object === level.segments || !wedgeClip);
      });
      if (show) {
        clipWedgeLevel(level);
      }
    });
  }

  // Each of a zone's lines from where it starts to the galaxy's edge, or
  // with a clip, only the part across the clip (none for a line whose
  // bearing is outside it).
  function clipWedgeLevel(level) {
    var attribute = level.segments.geometry.getAttribute("position");
    var points = attribute.array;
    var c = wedgeClip;
    var marginR = 0;
    var marginA = c ? WEDGE_CLIP_EPSILON : 0;
    level.lines.forEach(function (line, n) {
      var cos = Math.cos(line.angleRad);
      var sin = Math.sin(line.angleRad);
      var near = line.r0;
      var far = GALAXY_EDGE;
      if (c && c.cells) {
        // Only across the blocks the line's bearing passes through.
        var lo = Infinity;
        var hi = -Infinity;
        c.cells.forEach(function (b) {
          if (bearingFrom(line.angleRad, b.t0 - marginA) <= b.t1 - b.t0 + 2 * marginA) {
            lo = Math.min(lo, b.r0);
            hi = Math.max(hi, b.r1);
          }
        });
        if (lo > hi) {
          far = near;
        } else {
          near = Math.max(near, lo);
          far = Math.min(far, hi);
        }
      } else if (c) {
        var off = bearingFrom(line.angleRad, c.a0 - marginA);
        if (off > c.a1 - c.a0 + 2 * marginA) {
          far = near;
        } else {
          near = Math.max(near, c.r0 - marginR);
          far = Math.min(far, c.r1 + marginR);
        }
      }
      far = Math.max(near, far);
      points.set([near * cos, near * sin, 0, far * cos, far * sin, 0], 6 * n);
    });
    attribute.needsUpdate = true;
  }

  // Keeps the labels WEDGE_LABEL_PX tall on screen: a sprite without size
  // attenuation spans scale / tan(fov / 2) half-heights of the view.
  function updateWedgeLabels() {
    var heightPx = canvasEl.clientHeight || 1;
    var h = (WEDGE_LABEL_PX * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / heightPx;
    wedgeLabels.forEach(function (label) {
      label.scale.set(h * label.material.map.userData.aspect, h, 1);
    });
  }

  // --- The block solid -------------------------------------------------------
  //
  // A stage's blocks (galaxyprisms.js) fill their whole cells; a thin
  // brighter outline along each face's own edges (found from the face's
  // 0..1 uv, about a screen pixel wide) keeps neighboring blocks apart.
  // Each part is two meshes sharing one shader:
  // - solid: blocks whose every sector is generated, opaque;
  // - glass: everything else, translucent (no depth writes), its blocks
  //   sorted back to front from the camera when built. Unfilled space is
  //   10% to 30% opaque by density, and a block grows more solid with its
  //   filled share (galaxyblocks.blockOpacity).
  // `fade` scales a mesh's opacity (the drill-down's fades and dimming).
  // The logdepthbuf chunks match the renderer's logarithmic depth buffer.
  // Filled blocks (MAP.37: "a much higher contrast"): a saturated amber
  // over most of the face and on the edges, nothing like the blue-white
  // density ramp, the same for one generated sector as for many (their
  // opacity still tells them apart). A deeper shade on the light theme's
  // pale background.
  var FILLED_TINT = new THREE.Color(isLightBackground() ? "#d06a00" : "#ffb02e");
  var FILLED_FACE_MIX = 0.85;
  // How far an unfilled block's edges brighten toward white: the ring,
  // wedge and layer boundaries of the grid, kept faint like the wedge
  // lines (a filled block's amber edges add up to 0.7 more).
  var GRID_EDGE_MIX = 0.1;

  function makeBlockMaterial(translucent) {
    return new THREE.ShaderMaterial({
      uniforms: {
        filledTint: { value: FILLED_TINT }, fade: { value: 1 },
      },
      transparent: translucent,
      depthWrite: !translucent,
      vertexShader: [
        "#include <common>",
        "#include <logdepthbuf_pars_vertex>",
        "attribute vec3 prismColor;",
        "attribute float prismAlpha;",
        "attribute float prismFill;",
        "attribute vec2 faceUv;",
        "varying vec3 vColor;",
        "varying float vAlpha;",
        "varying float vFill;",
        "varying vec2 vUv;",
        "void main() {",
        "  vColor = prismColor;",
        "  vAlpha = prismAlpha;",
        "  vFill = prismFill;",
        "  vUv = faceUv;",
        "  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);",
        "  gl_Position = projectionMatrix * mvPosition;",
        "  #include <logdepthbuf_vertex>",
        "}",
      ].join("\n"),
      fragmentShader: [
        "#include <common>",
        "#include <logdepthbuf_pars_fragment>",
        "varying vec3 vColor;",
        "uniform vec3 filledTint;",
        "uniform float fade;",
        "varying float vAlpha;",
        "varying float vFill;",
        "varying vec2 vUv;",
        "void main() {",
        "  #include <logdepthbuf_fragment>",
        "  vec2 toEdge = min(vUv, 1.0 - vUv) / max(fwidth(vUv), vec2(1e-6));",
        "  float edge = 1.0 - smoothstep(0.5, 1.5, min(toEdge.x, toEdge.y));",
        "  vec3 face = mix(vColor, filledTint, step(0.001, vFill) * " + FILLED_FACE_MIX.toFixed(3) + ");",
        "  vec3 edgeColor = mix(vec3(1.0), filledTint, step(0.001, vFill));",
        "  gl_FragColor = vec4(mix(face, edgeColor, (" + GRID_EDGE_MIX.toFixed(3) + " + 0.7 * vFill) * edge), vAlpha * fade);",
        "  #include <colorspace_fragment>",
        "}",
      ].join("\n"),
    });
  }

  var highlightSprite = new THREE.Sprite(
    new THREE.SpriteMaterial({ map: highlightTexture, transparent: true, depthWrite: false, depthTest: false })
  );
  highlightSprite.visible = false;
  highlightSprite.renderOrder = 5;
  scene.add(highlightSprite);

  function highlightPosition(x, y, z) {
    highlightSprite.position.set(x, y, z);
    highlightSprite.visible = true;
  }

  // Keeps the selection ring HIGHLIGHT_PX across on screen whatever its
  // distance: a sprite spans its scale in world units, so the scale is the
  // world size of that many pixels at the sprite's own distance.
  function updateHighlightScale() {
    if (!highlightSprite.visible) {
      return;
    }
    var heightPx = canvasEl.clientHeight || 1;
    var worldPerPixel = (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2) * camera.position.distanceTo(highlightSprite.position)) / heightPx;
    var size = HIGHLIGHT_PX * 2 * worldPerPixel;
    highlightSprite.scale.set(size, size, 1);
  }

  // Density shading's dim-to-accent-to-white ramp (galaxyblocks.js).
  var PRISM_DIM = new THREE.Color(0x1d2340);
  var PRISM_ACCENT = new THREE.Color(accentColor);
  var PRISM_HOT = new THREE.Color(0xeef0ff);
  var galaxyShape = data.densityShape || null;
  var edgePc = data.edgePc || 1;
  // How many sectors a side the drawn blocks are.
  var drawnSectorsPerPrism = 1;

  // Parsecs per screen pixel at the view's focus (the orbit target), for
  // the camera at orbit radius `radius` (default: where it is). Per CSS
  // pixel, on purpose, not per device pixel: the block size (a block at
  // least blockMinPx across) and the scale readout's "1 px" are about what
  // a person can see and click, which is the same on a 2x phone screen as
  // on a 1x monitor. The renderer draws at up to 2 device pixels per CSS
  // pixel (setPixelRatio) for sharpness only, and the bright stars' sizes
  // are CSS pixels too (scaled by pixelRatio in their shader).
  function pcPerPixelAtTarget(radius) {
    var heightPx = canvasEl.clientHeight || 1;
    return ((radius || orbit.radius) * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / heightPx;
  }

  // --- Blocks ------------------------------------------------------------------
  //
  // The drill-down (galaxystageview.js) draws each stage's blocks with
  // galaxyblocks.js's scene and the block shader above.
  var reducedMotion = typeof window.matchMedia === "function"
    && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  var blockSceneConfig = {
    edgePc: edgePc, galaxyRadius: GALAXY_RADIUS, shape: galaxyShape, minPx: data.blockMinPx, budget: data.blockBudget,
    palette: {
      dim: PRISM_DIM.toArray(), accent: PRISM_ACCENT.toArray(), hot: PRISM_HOT.toArray(),
      placedLow: new THREE.Color(PLACED_LOW_DENSITY_COLOR).toArray(),
      placedHigh: new THREE.Color(PLACED_HIGH_DENSITY_COLOR).toArray(),
    },
  };
  var localBlockScene = createBlockScene(blockSceneConfig);

  // The view ball for the camera at orbit radius `radius` around `center`
  // (the fetched tiles' reach).
  function viewRadiusFor(radius) {
    return Math.max(MIN_RADIUS, Math.min(MAX_RADIUS * FETCH_RADIUS_FACTOR, radius * FETCH_RADIUS_FACTOR));
  }

  // One packed part (galaxyblocks' pack) as a mesh with the block shader,
  // keeping the part for picking.
  function makeBlockMesh(part, translucent) {
    var geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(part.positions, 3));
    geometry.setAttribute("prismCenter", new THREE.BufferAttribute(part.centers, 3));
    geometry.setAttribute("prismColor", new THREE.BufferAttribute(part.colors, 3, true));
    geometry.setAttribute("faceUv", new THREE.BufferAttribute(part.uvs, 2, true));
    geometry.setAttribute("prismAlpha", new THREE.BufferAttribute(part.alphas, 1, true));
    geometry.setAttribute("prismFill", new THREE.BufferAttribute(part.fills, 1, true));
    geometry.setIndex(new THREE.BufferAttribute(part.indices, 1));
    var mesh = new THREE.Mesh(geometry, makeBlockMaterial(translucent));
    mesh.frustumCulled = false;
    mesh.visible = part.vertexCount > 0;
    mesh.userData.part = part;
    return mesh;
  }

  // --- Cube tiles ----------------------------------------------------------
  //
  // The map asks for fixed cubes of space ("tiles") rather than "everything
  // within R of the camera target" -- see stellarObjects.galaxyViewport's
  // "Cube tiles" section. Space is an octree: level 0 is one cube
  // tileRootEdgePc on a side centered on the galactic origin, each level
  // halves the edge, and a tile's key is "level/ix/iy/iz". neededTiles()
  // picks the smallest level whose tiles are at least the view radius
  // across (so at most 27 tiles cover the view). lib/galaxymap3d.py's
  // initial_tile_request does the same in Python for the first frame.
  // Once the view's own tiles are in, the tiles one click-zoom step in and
  // out are fetched too (prefetchTiles), so zooming never waits on them.
  //
  // A tile's contents depend only on its key and the database's content
  // stamp, so tiles are cached three ways: in memory here, in this
  // browser's localStorage (so a revisit or reload doesn't refetch them),
  // and on the server's disk (lib/tilecache.py) -- only tiles in none of
  // those reach the API and database. Every response carries the current
  // stamp; when it differs from ours, the response also lists which tiles
  // changed since (history), and only those are dropped and refetched. A
  // change the history can't account for (a deleted sector, a new
  // release, a stamp older than the history) drops every cached tile.

  var TILE_ROOT = data.tileRootEdgePc || 65536;
  var TILE_MAX_LEVEL = data.tileMaxLevel != null ? data.tileMaxLevel : 12;
  var FETCH_RADIUS_FACTOR = data.fetchRadiusFactor || 1.6;
  var MAX_TILES_PER_REQUEST = data.maxTilesPerRequest || 128;

  function tileLevelForRadius(radius) {
    if (!(radius > 0)) {
      return TILE_MAX_LEVEL;
    }
    var level = Math.floor(Math.log2(TILE_ROOT / radius));
    return Math.max(0, Math.min(TILE_MAX_LEVEL, level));
  }

  function tileEdge(level) {
    return TILE_ROOT / Math.pow(2, level);
  }

  // Nearest first, matching galaxyViewport.tiles_intersecting_sphere.
  function tilesIntersectingSphere(level, center, radius) {
    var edge = tileEdge(level);
    var origin = -TILE_ROOT / 2;
    var span = Math.pow(2, level);
    var c = [center.x, center.y, center.z];
    var lo = [];
    var hi = [];
    for (var axis = 0; axis < 3; axis++) {
      lo.push(Math.max(0, Math.floor((c[axis] - radius - origin) / edge)));
      hi.push(Math.min(span - 1, Math.floor((c[axis] + radius - origin) / edge)));
    }
    var found = [];
    for (var ix = lo[0]; ix <= hi[0]; ix++) {
      for (var iy = lo[1]; iy <= hi[1]; iy++) {
        for (var iz = lo[2]; iz <= hi[2]; iz++) {
          var index = [ix, iy, iz];
          var distanceSq = 0;
          for (var a = 0; a < 3; a++) {
            var boxLo = origin + index[a] * edge;
            var nearest = Math.min(Math.max(c[a], boxLo), boxLo + edge);
            distanceSq += (c[a] - nearest) * (c[a] - nearest);
          }
          if (distanceSq <= radius * radius) {
            found.push({ d: distanceSq, key: level + "/" + ix + "/" + iy + "/" + iz });
          }
        }
      }
    }
    found.sort(function (p, q) { return p.d - q.d; });
    return found.map(function (f) { return f.key; });
  }

  // The tiles the camera at orbit radius `radius` (default: where it is)
  // needs.
  function neededTiles(radius) {
    var viewRadius = viewRadiusFor(radius || orbit.radius);
    var level = tileLevelForRadius(viewRadius);
    var keys = tilesIntersectingSphere(level, target, viewRadius);
    return { keys: keys, viewRadius: viewRadius };
  }

  // --- Tile caches ---------------------------------------------------------

  var TILE_MEMORY_MAX = 2000;
  // Keyed by database name, the same as before the page moved to
  // /galaxy, so a visitor's stored tiles carry over.
  var STORAGE_PREFIX = "planetgen:tile:" + data.storageKey + ":";
  // Where the stamp our stored tiles are at is kept, and the generation
  // (the label our stored tiles are filed under since they last all went
  // stale).
  var STAMP_RECORD = STORAGE_PREFIX + "@stamp";
  var remembered = readStampRecord();
  var currentStamp = remembered.stamp || "";
  var currentGeneration = remembered.generation || "";
  var tileMemory = new Map();

  // Map iteration order is insertion order, so re-inserting on every hit
  // and evicting from the front makes these LRU caches.
  function memoryGet(cache, key) {
    if (!cache.has(key)) {
      return undefined;
    }
    var value = cache.get(key);
    cache.delete(key);
    cache.set(key, value);
    return value;
  }

  function memorySet(cache, key, value, max) {
    cache.delete(key);
    cache.set(key, value);
    while (cache.size > max) {
      cache.delete(cache.keys().next().value);
    }
  }

  // localStorage can be missing, full, or throw on any access (private
  // browsing, blocked site data) -- every use is wrapped, and the map
  // works the same without it, just refetching from the server cache.
  function storage() {
    try {
      return window.localStorage || null;
    } catch (err) {
      return null;
    }
  }

  function readStampRecord() {
    var store = storage();
    try {
      var record = store ? JSON.parse(store.getItem(STAMP_RECORD) || "{}") : {};
      return record && typeof record === "object" ? record : {};
    } catch (err) {
      return {};
    }
  }

  function rememberStamp() {
    var store = storage();
    if (!store) {
      return;
    }
    try {
      store.setItem(STAMP_RECORD, JSON.stringify({ stamp: currentStamp, generation: currentGeneration }));
    } catch (err) {
      // Next visit just starts afresh.
    }
  }

  // Drops this database's stored tiles, except the current generation's
  // when keepCurrent is set.
  function purgeStoredTiles(keepCurrent) {
    var store = storage();
    if (!store) {
      return;
    }
    try {
      var doomed = [];
      for (var i = 0; i < store.length; i++) {
        var name = store.key(i);
        if (name && name.indexOf(STORAGE_PREFIX) === 0 && name !== STAMP_RECORD) {
          if (!keepCurrent || name.indexOf(STORAGE_PREFIX + currentGeneration + ":") !== 0) {
            doomed.push(name);
          }
        }
      }
      doomed.forEach(function (name) { store.removeItem(name); });
    } catch (err) {
      // Nothing more to do -- storage just stays as it was.
    }
  }

  function storageName(key) {
    return STORAGE_PREFIX + currentGeneration + ":" + key;
  }

  function removeStoredTile(key) {
    var store = storage();
    if (!store || !currentGeneration) {
      return;
    }
    try {
      store.removeItem(storageName(key));
    } catch (err) {
      // Nothing more to do.
    }
  }

  function storedTile(key) {
    var store = storage();
    if (!store || !currentGeneration) {
      return undefined;
    }
    try {
      var raw = store.getItem(storageName(key));
      return raw ? JSON.parse(raw) : undefined;
    } catch (err) {
      return undefined;
    }
  }

  function storeTile(key, tile) {
    var store = storage();
    if (!store || !currentGeneration) {
      return;
    }
    var raw = JSON.stringify(tile);
    try {
      store.setItem(storageName(key), raw);
    } catch (err) {
      // Probably full: drop other stamps' leftovers and try once more,
      // then drop this database's tiles entirely and try a last time.
      purgeStoredTiles(true);
      try {
        store.setItem(storageName(key), raw);
      } catch (err2) {
        purgeStoredTiles(false);
        try {
          store.setItem(storageName(key), raw);
        } catch (err3) {
          // Leave it memory-only.
        }
      }
    }
  }

  function getTile(key) {
    var tile = memoryGet(tileMemory, key);
    if (tile === undefined) {
      tile = storedTile(key);
      if (tile !== undefined) {
        memorySet(tileMemory, key, tile, TILE_MEMORY_MAX);
      }
    }
    return tile;
  }

  function putTile(key, tile) {
    memorySet(tileMemory, key, tile, TILE_MEMORY_MAX);
    storeTile(key, tile);
  }

  // The tile keys that changed between fromStamp and toStamp, from a
  // response's history (checks oldest first, each {from, to, tiles}), or
  // null when the history doesn't reach back to fromStamp unbroken.
  function staleKeys(fromStamp, toStamp, history) {
    if (!fromStamp || !Array.isArray(history)) {
      return null;
    }
    var start = -1;
    for (var i = history.length - 1; i >= 0; i--) {
      if (history[i] && history[i].from === fromStamp) {
        start = i;
        break;
      }
    }
    if (start < 0) {
      return null;
    }
    var keys = [];
    var at = fromStamp;
    for (var j = start; j < history.length; j++) {
      var entry = history[j];
      if (!entry || entry.from !== at || !Array.isArray(entry.tiles)) {
        return null;
      }
      keys = keys.concat(entry.tiles);
      at = entry.to;
    }
    return at === toStamp ? keys : null;
  }

  function adoptStamp(payload) {
    var stamp = payload.stamp;
    if (!stamp || stamp === currentStamp) {
      return false;
    }
    var stale = currentGeneration ? staleKeys(currentStamp, stamp, payload.history) : null;
    var hadStamp = !!currentStamp;
    currentStamp = stamp;
    if (hadStamp && stageView) {
      stageView.invalidate();
    }
    if (stale) {
      stale.forEach(function (key) {
        tileMemory.delete(key);
        removeStoredTile(key);
      });
    } else {
      currentGeneration = payload.generation || stamp;
      tileMemory.clear();
      purgeStoredTiles(true);
    }
    rememberStamp();
    return true;
  }

  function absorb(payload) {
    if (!payload) {
      return false;
    }
    var stampChanged = adoptStamp(payload);
    Object.keys(payload.tiles || {}).forEach(function (key) {
      putTile(key, payload.tiles[key]);
    });
    return stampChanged;
  }

  absorb(data.initial);
  purgeStoredTiles(true);

  // --- Clouds ----------------------------------------------------------------
  //
  // Each tile lists the nebulae and supernova remnants that reach into it
  // (queryDb.galaxy_clouds_in_box); a cloud is a sprite the size of its
  // sphere, drawn over the blocks without a depth test (a flat sprite
  // through a cloud's middle would otherwise be cut in half by the blocks
  // around it, though the cloud fills that space). One that is only a pixel or two across
  // is hidden, and one the camera comes close to fades out so it never
  // fills the screen.
  var CLOUD_MIN_PX = 2;
  var CLOUD_OPACITY = 1;
  var CLOUD_NEAR_FADE = [1.2, 2.0];
  // Clicking picks a cloud only while it is small enough on screen to aim
  // at; up close, clicks go to the sectors inside it.
  var CLOUD_PICK_MAX_PX = 150;
  var cloudGroup = new THREE.Group();
  cloudGroup.renderOrder = 4;
  scene.add(cloudGroup);
  var cloudSprites = new Map();
  var cloudSignature = "";

  function cloudKey(cloud) {
    return cloud.type + ":" + cloud.id;
  }

  // Shows exactly `clouds` (one entry per type and id).
  function setClouds(clouds) {
    var wanted = new Map();
    clouds.forEach(function (cloud) {
      wanted.set(cloudKey(cloud), cloud);
    });
    cloudSprites.forEach(function (sprite, key) {
      if (!wanted.has(key)) {
        cloudGroup.remove(sprite);
        sprite.material.dispose();
        cloudSprites.delete(key);
      }
    });
    wanted.forEach(function (cloud, key) {
      if (cloudSprites.has(key)) {
        return;
      }
      var sprite = new THREE.Sprite(new THREE.SpriteMaterial({
        map: cloudTexture(cloud), transparent: true, depthWrite: false, depthTest: false, opacity: CLOUD_OPACITY,
      }));
      sprite.position.set(cloud.x, cloud.y, cloud.z);
      sprite.scale.set(2 * cloud.radius_pc, 2 * cloud.radius_pc, 1);
      sprite.renderOrder = 4;
      sprite.userData.cloud = cloud;
      cloudGroup.add(sprite);
      cloudSprites.set(key, sprite);
    });
  }

  // A cloud's radius in pixels from where the camera is.
  function cloudPixels(cloud, distance) {
    var heightPx = canvasEl.clientHeight || 1;
    return (cloud.radius_pc * heightPx) / (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2) * Math.max(distance, 1e-6));
  }

  function updateClouds() {
    cloudSprites.forEach(function (sprite) {
      var cloud = sprite.userData.cloud;
      var distance = camera.position.distanceTo(sprite.position);
      var near = THREE.MathUtils.smoothstep(distance / cloud.radius_pc, CLOUD_NEAR_FADE[0], CLOUD_NEAR_FADE[1]);
      sprite.visible = near > 0 && cloudPixels(cloud, distance) >= CLOUD_MIN_PX
        && inWedgeClip(sprite.position.x, sprite.position.y, sprite.position.z);
      sprite.material.opacity = CLOUD_OPACITY * near;
    });
  }

  // The smallest visible cloud under a screen point that is still small
  // enough to aim at, as {cloud, core} (core: the point is within the
  // inner CLOUD_CORE of its radius), or null.
  var CLOUD_CORE = 0.5;

  function cloudAtClientPoint(clientX, clientY) {
    var ndc = ndcFromClientPoint(clientX, clientY);
    if (!ndc) {
      return null;
    }
    raycaster.setFromCamera(ndc, camera);
    var ray = raycaster.ray;
    var best = null;
    var bestOffset = 0;
    cloudSprites.forEach(function (sprite) {
      var cloud = sprite.userData.cloud;
      if (!sprite.visible || (best && best.radius_pc <= cloud.radius_pc)) {
        return;
      }
      var distance = camera.position.distanceTo(sprite.position);
      if (cloudPixels(cloud, distance) > CLOUD_PICK_MAX_PX) {
        return;
      }
      var offset = Math.sqrt(ray.distanceSqToPoint(sprite.position));
      if (offset <= cloud.radius_pc && ray.direction.dot(new THREE.Vector3().subVectors(sprite.position, ray.origin)) > 0) {
        best = cloud;
        bestOffset = offset;
      }
    });
    return best ? { cloud: best, core: bestOffset <= CLOUD_CORE * best.radius_pc } : null;
  }

  // --- Stars ---------------------------------------------------------------
  //
  // Each tile lists its most luminous pre-placed stars (bright_stars, v43:
  // every star of 500 L☉ or more, placed before any sector is filled), so
  // the arms show before anything is generated, and (MAP.51) the most
  // luminous stars of its generated systems down to a floor that drops
  // fourfold with each finer tile level (queryDb.generated_star_floor_sol):
  // closing in on filled sectors brings out fainter and fainter stars,
  // down to the red dwarfs at sector depth. Boss: "no matter how far the
  // user zooms in they should be very small with a big glow", and "rough
  // sizes relative to the size of the star, brightness relative to
  // luminosity, and color relative to temperature" -- each is a point of
  // light a fixed number of pixels across at every zoom (never sized by
  // distance): a core sized by the star's radius (STAR_CORE_PX, a red
  // dwarf one pixel, a supergiant four or five) in a soft halo whose width
  // and strength grow with luminosity (STAR_MIN_PX to STAR_MAX_PX), the
  // core dimmer for a fainter star, all in the star's blackbody color;
  // the faint end is drawn brighter by the Sector Map's own curve (MAP.87,
  // static/starlight.js). Halos blend normally rather than adding up, so a crowded arm zoomed
  // out glows in its stars' colors instead of burning to white; stars are
  // drawn over the blocks, never hidden by them.
  var STAR_MIN_PX = 6;
  var STAR_MAX_PX = 30;
  var STAR_CORE_PX = [1.1, 4.5];
  var STAR_LOG_RADIUS = [-1, 3];
  var STAR_GLOW = [0.2, 0.6];
  var STAR_CORE_ALPHA = [0.75, 1];
  var STAR_LOG_LUMINOSITY = [-4, 6];
  // A click within this many pixels of a star's center picks it.
  var STAR_PICK_PX = 7;
  // Stars new to the view fade in over this long (MAP.48: nothing pops
  // in); with reduced motion they just appear.
  var STAR_FADE_IN_MS = reducedMotion ? 0 : 300;

  var starMaterial = new THREE.ShaderMaterial({
    uniforms: {
      pixelRatio: { value: renderer.getPixelRatio() }, now: { value: 0 },
      fadeIn: { value: Math.max(STAR_FADE_IN_MS, 1) / 1000 },
    },
    vertexShader: [
      "#include <common>",
      "#include <logdepthbuf_pars_vertex>",
      "attribute float starBorn;",
      "uniform float now;",
      "uniform float fadeIn;",
      "varying float vShown;",
      "attribute float starSize;",
      "attribute float starCore;",
      "attribute float starGlow;",
      "attribute float starBright;",
      "attribute vec3 starColor;",
      "uniform float pixelRatio;",
      "attribute float starClipped;",
      "varying vec3 vColor;",
      "varying float vCore;",
      "varying float vGlow;",
      "varying float vBright;",
      "void main() {",
      "  vColor = starColor;",
      "  vBright = starBright;",
      "  vCore = starCore / starSize;",
      "  vGlow = starGlow;",
      "  vShown = clamp((now - starBorn) / fadeIn, 0.0, 1.0);",
      "  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);",
      "  gl_PointSize = starSize * pixelRatio;",
      // Outside the wedge shown (setWedgeClip): dropped.
      "  if (starClipped > 0.5) {",
      "    gl_PointSize = 0.0;",
      "    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);",
      "  }",
      "  #include <logdepthbuf_vertex>",
      "}",
    ].join("\n"),
    fragmentShader: [
      "#include <common>",
      "#include <logdepthbuf_pars_fragment>",
      "varying vec3 vColor;",
      "varying float vCore;",
      "varying float vGlow;",
      "varying float vShown;",
      "varying float vBright;",
      "void main() {",
      "  #include <logdepthbuf_fragment>",
      "  float r = length(gl_PointCoord * 2.0 - 1.0);",
      "  if (r > 1.0) discard;",
      "  float core = 1.0 - smoothstep(vCore * 0.5, vCore, r);",
      "  float halo = vGlow * exp(-r * r * 4.0) * (1.0 - r);",
      "  gl_FragColor = vec4(mix(vColor, vec3(1.0), core * 0.6 * vBright), clamp(core * vBright + halo, 0.0, 1.0) * vShown);",
      "}",
    ].join("\n"),
    transparent: true,
    depthWrite: false,
    // Drawn over the blocks (MAP.51): a generated sector's block is
    // solid, and depth-tested stars inside it never showed.
    depthTest: false,
  });
  var starPoints = new THREE.Points(new THREE.BufferGeometry(), starMaterial);
  starPoints.renderOrder = 5;
  starPoints.frustumCulled = false;
  scene.add(starPoints);
  var starList = [];
  // When each drawn star first showed (seconds, performance.now's clock),
  // by starKey, so a star already on screen never fades in again.
  var starBornAt = new Map();

  function starClock() {
    return performance.now() / 1000;
  }

  // A star's key among the drawn ones: bright_stars and stars rows are
  // numbered separately.
  function starKey(star) {
    return (star.generated ? "s" : "b") + star.id;
  }

  // Where `value`'s log10 falls in `range`, 0..1.
  function logShare(value, range) {
    return THREE.MathUtils.clamp((Math.log10(Math.max(value, 1e-12)) - range[0]) / (range[1] - range[0]), 0, 1);
  }

  // Draws exactly `stars` (one entry per starKey).
  function setStars(stars) {
    starList = stars;
    var n = stars.length;
    var now = starClock();
    var born = new Float32Array(n);
    var bornAt = new Map();
    stars.forEach(function (star, i) {
      var key = starKey(star);
      var at = starBornAt.has(key) ? starBornAt.get(key) : now;
      bornAt.set(key, at);
      born[i] = at;
    });
    starBornAt = bornAt;
    var positions = new Float32Array(3 * n);
    var colors = new Float32Array(3 * n);
    var sizes = new Float32Array(n);
    var cores = new Float32Array(n);
    var glows = new Float32Array(n);
    var brights = new Float32Array(n);
    stars.forEach(function (star, i) {
      var t = logShare(star.luminosity_sol, STAR_LOG_LUMINOSITY);
      // Without a stored radius, guess one from the luminosity.
      var r = logShare(star.radius_sol != null ? star.radius_sol : Math.pow(star.luminosity_sol, 0.35), STAR_LOG_RADIUS);
      positions.set([star.x, star.y, star.z], 3 * i);
      colors.set(starColor(star.temperature_k), 3 * i);
      cores[i] = THREE.MathUtils.lerp(STAR_CORE_PX[0], STAR_CORE_PX[1], r);
      // MAP.87: the faint end drawn brighter (static/starlight.js).
      var halo = boostLight({
        sizePx: THREE.MathUtils.lerp(STAR_MIN_PX, STAR_MAX_PX, t * t),
        glow: THREE.MathUtils.lerp(STAR_GLOW[0], STAR_GLOW[1], t),
      }, starLightBoost(star.luminosity_sol));
      sizes[i] = Math.max(halo.sizePx, 2 * cores[i] + 2);
      glows[i] = halo.glow;
      brights[i] = THREE.MathUtils.lerp(STAR_CORE_ALPHA[0], STAR_CORE_ALPHA[1], t);
    });
    var geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute("starColor", new THREE.BufferAttribute(colors, 3));
    geometry.setAttribute("starSize", new THREE.BufferAttribute(sizes, 1));
    geometry.setAttribute("starCore", new THREE.BufferAttribute(cores, 1));
    geometry.setAttribute("starGlow", new THREE.BufferAttribute(glows, 1));
    geometry.setAttribute("starBright", new THREE.BufferAttribute(brights, 1));
    geometry.setAttribute("starBorn", new THREE.BufferAttribute(born, 1));
    geometry.setAttribute("starClipped", new THREE.BufferAttribute(new Float32Array(n), 1));
    starPoints.geometry.dispose();
    starPoints.geometry = geometry;
    markClippedStars();
  }

  // Marks the stars outside the wedge shown, which the shader drops.
  function markClippedStars() {
    var attribute = starPoints.geometry.getAttribute("starClipped");
    if (!attribute) {
      return;
    }
    starList.forEach(function (star, i) {
      attribute.array[i] = inWedgeClip(star.x, star.y, star.z) ? 0 : 1;
    });
    attribute.needsUpdate = true;
  }

  // The star whose center is nearest a screen point, within
  // STAR_PICK_PX, or null.
  function starAtClientPoint(clientX, clientY) {
    var rect = canvasEl.getBoundingClientRect();
    if (!starList.length || !rect.width || !rect.height) {
      return null;
    }
    var projected = new THREE.Vector3();
    var best = null;
    var bestPx = STAR_PICK_PX;
    starList.forEach(function (star) {
      if (!inWedgeClip(star.x, star.y, star.z)) {
        return;
      }
      projected.set(star.x, star.y, star.z).project(camera);
      if (projected.z < -1 || projected.z > 1) {
        return;
      }
      var px = Math.hypot(
        rect.left + ((projected.x + 1) / 2) * rect.width - clientX,
        rect.top + ((1 - projected.y) / 2) * rect.height - clientY);
      if (px <= bestPx) {
        best = star;
        bestPx = px;
      }
    });
    return best;
  }

  // --- Drawing from tiles --------------------------------------------------

  var starSignature = "";

  // A tile's box, [lo, hi] in parsecs.
  function tileBox(key) {
    var parts = key.split("/").map(Number);
    var edge = tileEdge(parts[0]);
    var origin = -TILE_ROOT / 2;
    var lo = [origin + parts[1] * edge, origin + parts[2] * edge, origin + parts[3] * edge];
    return [lo, [lo[0] + edge, lo[1] + edge, lo[2] + edge]];
  }

  function inBox(box, star) {
    return star.x >= box[0][0] && star.x < box[1][0] && star.y >= box[0][1] && star.y < box[1][1]
      && star.z >= box[0][2] && star.z < box[1][2];
  }

  // The nearest cached tile holding `key`'s box (a coarser level), or
  // undefined. Looks in memory only, without touching its LRU order.
  function cachedAncestor(key) {
    var parts = key.split("/").map(Number);
    for (var level = parts[0] - 1; level >= 0; level--) {
      var shift = Math.pow(2, parts[0] - level);
      var tile = tileMemory.get(level + "/" + Math.floor(parts[1] / shift) + "/" + Math.floor(parts[2] / shift)
        + "/" + Math.floor(parts[3] / shift));
      if (tile !== undefined) {
        return tile;
      }
    }
    return undefined;
  }

  // Stars to show in a tile that hasn't arrived yet (MAP.48): the ones
  // already drawn there, and a cached coarser tile's there, so a zoom
  // keeps what was on screen instead of blanking until the new tiles come.
  function carriedStars(key, into) {
    var box = tileBox(key);
    starList.forEach(function (star) {
      if (!into.has(starKey(star)) && inBox(box, star)) {
        into.set(starKey(star), star);
      }
    });
    tileStars(cachedAncestor(key)).forEach(function (star) {
      if (!into.has(starKey(star)) && inBox(box, star)) {
        into.set(starKey(star), star);
      }
    });
  }

  // A tile's stars, bright and generated (each generated one marked),
  // or none for a missing tile.
  function tileStars(tile) {
    if (!tile) {
      return [];
    }
    if (!tile.allStars) {
      Object.defineProperty(tile, "allStars", {
        value: (tile.stars || []).concat((tile.generated || []).map(function (star) {
          return Object.assign({ generated: true }, star);
        })),
      });
    }
    return tile.allStars;
  }

  // Takes the filled sectors from whatever of the needed tiles is already
  // cached (the blocks are rebuilt from them when they change); returns
  // what's still missing.
  function renderFromCache(need) {
    var missing = [];
    var clouds = new Map();
    var stars = new Map();
    need.keys.forEach(function (key) {
      var tile = getTile(key);
      if (tile === undefined) {
        missing.push(key);
        return;
      }
      (tile.clouds || []).forEach(function (cloud) {
        clouds.set(cloudKey(cloud), cloud);
      });
      tileStars(tile).forEach(function (star) {
        stars.set(starKey(star), star);
      });
    });
    missing.forEach(function (key) {
      carriedStars(key, stars);
    });
    var starKeys = currentStamp + "|" + Array.from(stars.keys()).sort().join(",");
    if (starKeys !== starSignature) {
      starSignature = starKeys;
      setStars(Array.from(stars.values()));
    }
    var cloudKeys = Array.from(clouds.keys()).sort().join(",");
    if (cloudKeys !== cloudSignature) {
      cloudSignature = cloudKeys;
      setClouds(Array.from(clouds.values()));
    }
    return missing;
  }

  // --- Live tile fetching --------------------------------------------------

  var FETCH_DEBOUNCE_MS = 300;
  var fetchTimer = null;
  var activeAbort = null;

  // Caps how often a click/double-click/zoom-button interaction (every
  // caller that passes scheduleFetch(true) -- see each one below) can
  // actually trigger an IMMEDIATE live fetch: at most MAX_CLICKS_PER_SECOND
  // per second. `doFetch`'s own activeAbort.abort() only stops the
  // BROWSER from waiting on a superseded response -- it doesn't reliably
  // stop the server from finishing a request it already started, so
  // rapid clicking still costs the server work per click even when every
  // earlier response gets thrown away client-side. A click/double-click's
  // own visual effect (the camera recentering/zooming, via centerOn/
  // zoomInOnTarget) is never throttled here, only the network fetch that
  // follows it -- clicking faster than the cap still feels instant, it
  // just falls back to the standard debounced delay below instead of
  // firing right away, so a rapid burst still settles on exactly one
  // fetch shortly after it stops. (With tiles, an aborted request's work
  // isn't wasted either: the server caches whatever it computed.)
  var MAX_CLICKS_PER_SECOND = 4;
  var MIN_MS_BETWEEN_IMMEDIATE_FETCHES = 1000 / MAX_CLICKS_PER_SECOND;
  var lastImmediateFetchAt = 0;

  function scheduleFetch(immediate) {
    if (fetchTimer) {
      clearTimeout(fetchTimer);
    }
    var delay = FETCH_DEBOUNCE_MS;
    if (immediate) {
      var now = Date.now();
      if (now - lastImmediateFetchAt >= MIN_MS_BETWEEN_IMMEDIATE_FETCHES) {
        lastImmediateFetchAt = now;
        delay = 0;
      }
      // else: rate-limited -- falls through to the debounced delay above
      // instead of a bare no-op, so state still eventually syncs.
    }
    fetchTimer = setTimeout(doFetch, delay);
  }

  function doFetch() {
    fetchTimer = null;
    var need = neededTiles();
    var missing = renderFromCache(need);
    if (!missing.length) {
      prefetchTiles();
      return;
    }
    if (activeAbort) {
      activeAbort.abort();
    }
    var controller = typeof AbortController !== "undefined" ? new AbortController() : null;
    activeAbort = controller;
    var params = new URLSearchParams({
      tiles: missing.slice(0, MAX_TILES_PER_REQUEST).join(","),
    });
    if (currentStamp) {
      params.set("stamp", currentStamp);
    }
    fetch(data.fetchPath + "?" + params.toString(), controller ? { signal: controller.signal } : undefined)
      .then(function (response) {
        if (!response.ok) {
          throw new Error(data.fetchPath + " returned " + response.status);
        }
        return response.json();
      })
      .then(function (payload) {
        if (activeAbort === controller) {
          activeAbort = null;
        }
        var stampChanged = absorb(payload);
        var stillMissing = renderFromCache(neededTiles());
        // A new stamp dropped changed (or all) cached tiles, and a view
        // needing more than one request's worth of tiles has more to
        // fetch -- either way, go again (tiles already fetched are
        // cached by now).
        if (stampChanged || stillMissing.length) {
          scheduleFetch(false);
        } else {
          prefetchTiles();
        }
      })
      .catch(function (err) {
        if (err && err.name === "AbortError") {
          return;
        }
        // A transient fetch failure just leaves the currently-drawn
        // content in place -- the next camera move retries automatically,
        // and there's no useful place to surface a network error inside
        // this canvas.
      });
  }

  // The tiles one pick in and out (each about PREFETCH_ZOOM times closer
  // or farther), fetched at low priority once
  // the view's own are in. One prefetch runs at a time and a real fetch
  // never waits on or cancels it (or the other way round); its tiles just
  // land in the cache.
  var prefetching = false;
  var PREFETCH_ZOOM = 3;

  function prefetchTiles() {
    if (prefetching || activeAbort) {
      return;
    }
    var factor = PREFETCH_ZOOM;
    var missing = [];
    [orbit.radius / factor, orbit.radius * factor].forEach(function (radius) {
      neededTiles(Math.max(MIN_RADIUS, Math.min(MAX_RADIUS, radius))).keys.forEach(function (key) {
        if (missing.indexOf(key) < 0 && getTile(key) === undefined) {
          missing.push(key);
        }
      });
    });
    if (!missing.length) {
      return;
    }
    var params = new URLSearchParams({ tiles: missing.slice(0, MAX_TILES_PER_REQUEST).join(",") });
    if (currentStamp) {
      params.set("stamp", currentStamp);
    }
    prefetching = true;
    fetch(data.fetchPath + "?" + params.toString(), { priority: "low" })
      .then(function (response) {
        if (!response.ok) {
          throw new Error(data.fetchPath + " returned " + response.status);
        }
        return response.json();
      })
      .then(function (payload) {
        prefetching = false;
        // A new stamp means the view's own tiles may be stale too.
        if (absorb(payload)) {
          scheduleFetch(false);
        } else if (missing.length > MAX_TILES_PER_REQUEST) {
          prefetchTiles();
        }
      })
      .catch(function () {
        // Only a head start: the real fetch gets these when they're needed.
        prefetching = false;
      });
  }

  // Normally every tile is already cached (the page embeds them), so this
  // only fetches if the page's own tile fetch came up short.
  // First frame's tiles, then (normally nothing is missing) the prefetch.
  renderFromCache(neededTiles());
  scheduleFetch(false);

  // --- Territories -----------------------------------------------------------
  //
  // The Territories button draws who holds what (docs/design/population-
  // and-politics.md): each polity's reach as a translucent ball around
  // its capital, and the systems it owns as small dots in its own color
  // (/galaxy/territories, at most 20,000 nearest their capitals). Both
  // are cheap at galaxy scale -- a handful of balls, one points object --
  // and the dots keep a fixed size on screen, so zooming out never turns
  // them into a smear. Fetched once, the first time the button is
  // pressed, and listed under the map with each polity's system count.
  var TERRITORY_DOT_PX = 3.5;
  var TERRITORY_BALL_OPACITY = 0.1;
  var territoryGroup = null;
  var territoryLoading = false;
  var territoryWanted = false;
  var territoryDots = null;

  function lyToPc(ly) {
    return (ly * LIGHTYEAR_M) / PARSEC_M;
  }

  function territoryColor(hex) {
    return new THREE.Color(hex || accentColor);
  }

  // The polities' balls and the owned systems' dots, from the payload.
  function buildTerritories(payload) {
    var group = new THREE.Group();
    group.renderOrder = 3;
    (payload.polities || []).forEach(function (polity) {
      if (!polity.capital_pc || !(polity.reach_ly > 0)) {
        return;
      }
      var ball = new THREE.Mesh(
        new THREE.SphereGeometry(lyToPc(polity.reach_ly), 32, 16),
        new THREE.MeshBasicMaterial({
          color: territoryColor(polity.color), transparent: true, opacity: TERRITORY_BALL_OPACITY,
          depthWrite: false, side: THREE.BackSide,
        })
      );
      ball.position.set(polity.capital_pc[0], polity.capital_pc[1], polity.capital_pc[2]);
      ball.frustumCulled = false;
      group.add(ball);
    });
    var points = payload.points || [];
    if (points.length) {
      var positions = new Float32Array(points.length * 3);
      var colors = new Float32Array(points.length * 3);
      points.forEach(function (point, n) {
        positions[3 * n] = point.x;
        positions[3 * n + 1] = point.y;
        positions[3 * n + 2] = point.z;
        var color = territoryColor(point.color);
        colors[3 * n] = color.r;
        colors[3 * n + 1] = color.g;
        colors[3 * n + 2] = color.b;
      });
      var geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
      geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
      territoryDots = new THREE.Points(geometry, new THREE.PointsMaterial({
        size: TERRITORY_DOT_PX * renderer.getPixelRatio(), sizeAttenuation: false, vertexColors: true,
        transparent: true, depthWrite: false,
      }));
      territoryDots.frustumCulled = false;
      group.add(territoryDots);
    }
    return group;
  }

  // The legend under the map: one row per polity, with its color, its
  // government and how many systems it holds.
  function showTerritoryLegend(payload) {
    var box = document.getElementById("galaxymap3d-territories");
    if (!box) {
      return;
    }
    box.textContent = "";
    var named = (payload.polities || []).filter(function (polity) { return polity.name; });
    var heading = document.createElement("h3");
    heading.textContent = "Territories";
    box.appendChild(heading);
    if (!named.length) {
      var empty = document.createElement("p");
      empty.className = "hint";
      empty.textContent = (payload.polities || []).length
        ? "No polity has a name yet."
        : "No polity holds any system yet.";
      box.appendChild(empty);
      return;
    }
    var list = document.createElement("ul");
    named.sort(function (p, q) { return (q.system_count || 0) - (p.system_count || 0); });
    named.forEach(function (polity) {
      var item = document.createElement("li");
      var swatch = document.createElement("span");
      swatch.className = "galaxy-territory-swatch";
      swatch.style.background = polity.color || accentColor;
      var name = document.createElement("span");
      name.textContent = polity.name + (polity.government ? " · " + polity.government : "");
      var count = document.createElement("span");
      count.className = "galaxy-territory-count";
      count.textContent = formatNumber(polity.system_count || 0)
        + (polity.system_count === 1 ? " system" : " systems");
      item.appendChild(swatch);
      item.appendChild(name);
      item.appendChild(count);
      list.appendChild(item);
    });
    box.appendChild(list);
  }

  function setTerritories(on, button) {
    territoryWanted = on;
    var box = document.getElementById("galaxymap3d-territories");
    if (box) {
      box.hidden = !on;
    }
    if (territoryGroup) {
      territoryGroup.visible = on;
      return;
    }
    if (!on || territoryLoading) {
      return;
    }
    territoryLoading = true;
    if (button) {
      button.disabled = true;
    }
    fetch(data.territoryPath || "/galaxy/territories", { headers: { Accept: "application/json" } })
      .then(function (response) {
        if (!response.ok) {
          throw new Error("HTTP " + response.status);
        }
        return response.json();
      })
      .then(function (payload) {
        territoryGroup = buildTerritories(payload);
        territoryGroup.visible = territoryWanted;
        scene.add(territoryGroup);
        showTerritoryLegend(payload);
      })
      .catch(function () {
        showHint("The territories could not be loaded. Please try again shortly.", true);
        if (button) {
          button.setAttribute("aria-pressed", "false");
        }
        territoryWanted = false;
        if (box) {
          box.hidden = true;
        }
      })
      .then(function () {
        territoryLoading = false;
        if (button) {
          button.disabled = false;
        }
      });
  }

  // --- NAV course ---------------------------------------------------------------
  //
  // With ?course=<from>,<to> (the NAV result's "Show on Galaxy Map", see
  // the drill-down design's section 9.4), the page embeds the course's
  // waypoints in galaxy-frame parsecs (web/nav_page.galaxy_course). They
  // are drawn as one line through every stop with a ring and a name at
  // each, over everything else and never picked, in the stages and in
  // free look alike; the map opens on the smallest stage holding them
  // all (galaxystages.courseStage), unless the URL names a stage itself.
  // A course inside one sector has no line to draw: the map opens that
  // sector instead.
  var course = data.course && (data.course.points || []).length >= 2 ? data.course : null;
  var courseSectors = [];
  var courseMarkers = [];
  var COURSE_RING_PX = 11;
  var COURSE_LABEL_PX = 13;

  if (data.course && data.course.sector) {
    courseSectors = [data.course.sector];
  }
  if (course) {
    var courseGroup = new THREE.Group();
    courseGroup.renderOrder = 4;
    scene.add(courseGroup);
    var linePoints = course.points.map(function (point) { return new THREE.Vector3(point.x, point.y, point.z); });
    courseGroup.add(new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(linePoints),
      new THREE.LineBasicMaterial({ color: new THREE.Color(accentColor), depthTest: false, transparent: true })
    ));
    course.points.forEach(function (point, n) {
      var ends = n === 0 || n === course.points.length - 1;
      var ring = new THREE.Sprite(new THREE.SpriteMaterial({
        map: highlightTexture, transparent: true, depthWrite: false, depthTest: false,
        sizeAttenuation: false, opacity: ends ? 1 : 0.7,
      }));
      ring.position.copy(linePoints[n]);
      courseGroup.add(ring);
      courseMarkers.push({ sprite: ring, px: ends ? COURSE_RING_PX : COURSE_RING_PX * 0.6 });
      if (!ends) {
        return;
      }
      var label = new THREE.Sprite(new THREE.SpriteMaterial({
        map: makeLabelTexture(point.name, accentColor), transparent: true, depthWrite: false, depthTest: false,
        sizeAttenuation: false,
      }));
      label.position.copy(linePoints[n]);
      courseGroup.add(label);
      courseMarkers.push({ sprite: label, px: COURSE_LABEL_PX, label: true });
    });
    courseSectors = course.points.map(function (point) {
      return sectorAddressAt(point.x, point.y, point.z, edgePc);
    });
  }

  // Keeps every course marker the same size on screen, like the wedge
  // labels (a sprite without size attenuation spans scale / tan(fov / 2)
  // half-heights of the view).
  function updateCourseMarkers() {
    if (!courseMarkers.length) {
      return;
    }
    var heightPx = canvasEl.clientHeight || 1;
    var per = (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / heightPx;
    courseMarkers.forEach(function (marker) {
      var h = marker.px * per;
      var aspect = marker.label ? marker.sprite.material.map.userData.aspect : 1;
      marker.sprite.scale.set(h * aspect, h, 1);
      // A name sits above its ring instead of on it.
      if (marker.label) {
        marker.sprite.center.set(0.5, -0.4);
      }
    });
  }

  // --- Drill-down stages ------------------------------------------------------
  //
  // The map opens on the drill-down (galaxystageview.js): the galaxy in
  // level-243 blocks, a slab of them from above, a block's level-27
  // children, and so on down to a sector, always from straight above.
  // Pointer and key input goes to the stage view; tiles (stars, clouds)
  // and the wedge lines follow the camera.
  var stageView = createStageView({
    THREE: THREE, scene: scene, camera: camera, canvasEl: canvasEl,
    edgePc: edgePc, shape: galaxyShape, galaxyRadius: GALAXY_RADIUS, reducedMotion: reducedMotion,
    accentColor: accentColor, canGenerate: !!data.generate, courseSectors: courseSectors,
    blockScene: localBlockScene,
    makeBlockMesh: makeBlockMesh,
    setCamera: function (v) {
      target.set(v.target[0], v.target[1], v.target[2]);
      orbit.radius = v.dist;
      orbit.theta = v.theta;
      orbit.phi = v.phi;
      applyCamera();
      updateScaleBar();
    },
    fetchStage: function (query) {
      return fetch((data.stagePath || "/galaxy/stage") + query, { headers: { Accept: "application/json" } })
        .then(function (response) {
          if (!response.ok) {
            throw new Error("HTTP " + response.status);
          }
          return response.json();
        });
    },
    setBlockSize: function (m) {
      drawnSectorsPerPrism = m;
      updateScaleBar();
    },
    showBlockInfo: function (info) { showBlockInfo(info, edgePc); },
    showPlacedInfo: showPlacedInfo,
    showCellInfo: showCellInfo,
    showHint: function (text) { showHint(text, false); },
    showPointAt: function (x, y) { return showPointAt(x, y); },
    setWedgeClip: function (clip) { setWedgeClip(clip); },
    sectorUrl: function (id) { return sceneData.sectorUrl ? sectorUrl(id) : null; },
    locate: function (name) {
      return fetch((data.locatePath || "/galaxy/locate") + "?q=" + encodeURIComponent(name), {
        headers: { Accept: "application/json" },
      })
        .then(function (response) {
          if (!response.ok) {
            throw new Error("HTTP " + response.status);
          }
          return response.json();
        })
        .then(function (payload) { return payload.matches || []; });
    },
    els: {
      crumbs: document.getElementById("galaxymap3d-crumbs"),
      slabs: document.getElementById("galaxymap3d-slabs"),
      tooltip: document.getElementById("galaxymap3d-tooltip"),
      notice: document.getElementById("galaxymap3d-notice"),
      address: document.getElementById("galaxymap3d-address"),
      matches: document.getElementById("galaxymap3d-matches"),
      controls: document.getElementById("galaxymap3d-controls"),
    },
  });

  // --- Pointer and keys ------------------------------------------------------
  //
  // Everything goes to the drill-down, which turns and zooms the view
  // where it may.
  // A click on a bright star or a cloud small enough to aim at shows it
  // instead of picking what's under it.
  canvasEl.addEventListener("pointerdown", function (event) { stageView.onPointerDown(event); });
  canvasEl.addEventListener("pointermove", function (event) { stageView.onPointerMove(event); });
  canvasEl.addEventListener("pointerup", function (event) { stageView.onPointerUp(event); });
  canvasEl.addEventListener("pointercancel", function (event) { stageView.onPointerUp(event); });
  canvasEl.addEventListener("pointerleave", function (event) { stageView.onPointerLeave(event); });
  canvasEl.addEventListener("keydown", function (event) { stageView.onKey(event); });
  // Below the galaxy and its quarters the wheel zooms; above them it is
  // left to scroll the page.
  canvasEl.addEventListener("wheel", function (event) {
    if (stageView.onWheel(event)) {
      event.preventDefault();
    }
  }, { passive: false });

  var raycaster = new THREE.Raycaster();

  function ndcFromClientPoint(clientX, clientY) {
    var rect = canvasEl.getBoundingClientRect();
    if (rect.width === 0 || rect.height === 0) {
      return null;
    }
    return new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
  }

  // The star or cloud a click at a screen point means, shown in the panel
  // and ringed: true when there was one.
  function showPointAt(clientX, clientY) {
    var star = starAtClientPoint(clientX, clientY);
    if (star) {
      highlightPosition(star.x, star.y, star.z);
      showStarInfo(star);
      return true;
    }
    var found = cloudAtClientPoint(clientX, clientY);
    if (found && found.core) {
      highlightPosition(found.cloud.x, found.cloud.y, found.cloud.z);
      showCloudInfo(found.cloud);
      return true;
    }
    highlightSprite.visible = false;
    return false;
  }

  // No right-click action any more (see this file's own module
  // docstring) -- the browser's own default context menu is left alone,
  // rather than preventDefault()-ing it for nothing.

  var controlsEl = document.getElementById("galaxymap3d-controls");
  if (controlsEl) {
    // Several buttons: let them wrap rather than run off the side panel.
    controlsEl.style.flexWrap = "wrap";
    wireMapControls(controlsEl, mapControlHandlers({
      stageView: stageView,
      setTerritories: setTerritories,
      territoriesWanted: function () { return territoryWanted; },
      wedgeGroup: wedgeGroup,
    }));
  }

  // --- Scale bar -----------------------------------------------------------

  var scaleEl = document.getElementById("galaxymap3d-scale");

  function niceScaleValue(raw) {
    if (!isFinite(raw) || raw <= 0) {
      return 0;
    }
    var magnitude = Math.pow(10, Math.floor(Math.log10(raw)));
    var mantissa = raw / magnitude;
    var niceMantissa;
    if (mantissa < 1.5) niceMantissa = 1;
    else if (mantissa < 3.5) niceMantissa = 2;
    else if (mantissa < 7.5) niceMantissa = 5;
    else niceMantissa = 10;
    return niceMantissa * magnitude;
  }

  // Up to 3 significant figures, grouped: 0.0512, 3.4, 1,280.
  function formatCount(value) {
    if (!(value > 0)) return "0";
    if (value >= 100) return formatNumber(value);
    return String(Number(value.toPrecision(value >= 1 ? 3 : 2)));
  }

  function plural(count, word) {
    return count === 1 ? word : word + "s";
  }

  // One line of the readout: "<lead> <sectors> · <distance>", the distance
  // on the shared ladder (static/distance.js), which adds ly in
  // parentheses to parsec values.
  function scaleLine(lead, pc, suffix) {
    var sectors = pc / edgePc;
    var line = document.createElement("span");
    line.className = "starmap-scale-label";
    // Wraps on a phone instead of running off the canvas.
    line.style.whiteSpace = "normal";
    line.textContent = lead + " " + formatCount(sectors) + " " + plural(sectors, "sector") + (suffix || "")
      + " · " + formatDistancePc(pc);
    return line;
  }

  // Three lines: what one screen pixel spans at the focus, how big one
  // prism (block) is, and a bar of about SCALE_BAR_PX rounded to a nice
  // number of parsecs. All per CSS pixel, like the block size itself:
  // they're about what a person can see and click, not device pixels.
  var SCALE_BAR_PX = 70;
  if (scaleEl) {
    // A column instead of the other maps' one-line row (set here, not in
    // style.css: only this map's readout stacks).
    scaleEl.style.flexDirection = "column";
    scaleEl.style.alignItems = "flex-start";
    scaleEl.style.gap = "0.1rem";
    scaleEl.style.maxWidth = "calc(100% - 1.5rem)";
  }

  function updateScaleBar() {
    if (!scaleEl) {
      return;
    }
    var pcPerScreenPx = pcPerPixelAtTarget();
    var nicePc = niceScaleValue(SCALE_BAR_PX * pcPerScreenPx);
    if (!nicePc) {
      return;
    }
    var m = drawnSectorsPerPrism || 1;
    var block = scaleLine("1 block =", m * edgePc, m === 1 ? "" : " across (≈ " + formatCount(m * m * m) + ")");
    var barRow = document.createElement("span");
    barRow.style.display = "flex";
    barRow.style.alignItems = "center";
    barRow.style.gap = "0.4rem";
    var bar = document.createElement("span");
    bar.className = "starmap-scale-bar";
    bar.style.width = (nicePc / pcPerScreenPx).toFixed(1) + "px";
    barRow.appendChild(bar);
    barRow.appendChild(scaleLine("≈", nicePc));
    scaleEl.replaceChildren(scaleLine("1 px ≈", pcPerScreenPx), block, barRow);
  }

  // --- Resize/render loop ----------------------------------------------

  function resize() {
    var width = canvasEl.clientWidth || 1;
    var height = canvasEl.clientHeight || 1;
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    updateScaleBar();
  }

  if (typeof ResizeObserver !== "undefined" && viewport) {
    new ResizeObserver(resize).observe(viewport);
  }
  resize();
  window.addEventListener("resize", resize);

  // NAV's "Pick on Galaxy Map": endpoints live only in generated
  // sectors, so "Generated only" stays on (design doc section 9).
  if (data.pick) {
    stageView.setGeneratedOnly(true);
    var onlyButton = document.querySelector('#galaxymap3d-controls [data-action="generated-only"]');
    if (onlyButton) {
      onlyButton.setAttribute("aria-pressed", "true");
      onlyButton.disabled = true;
      onlyButton.title = "Always on while choosing a NAV start or destination";
    }
  }

  // The drill-down, at the stage the URL names.
  stageView.openFromLocation();
  stageView.setActive(true);

  // Every frame: the drill-down's camera move, then (when the view moved)
  // the stars and clouds from cached tiles, a fetch for missing ones and
  // the wedge lines.
  var lastViewSignature = "";

  (function animate(now) {
    requestAnimationFrame(animate);
    now = now || performance.now();
    stageView.step(now);
    var viewSignature = [target.x, target.y, target.z, orbit.radius, canvasEl.clientHeight].join(",");
    if (viewSignature !== lastViewSignature) {
      lastViewSignature = viewSignature;
      var missing = renderFromCache(neededTiles());
      if (missing.length && !fetchTimer && !activeAbort) {
        scheduleFetch(false);
      }
      updateWedgeLevels(pcPerPixelAtTarget());
    }
    starMaterial.uniforms.now.value = starClock();
    updateClouds();
    updateHighlightScale();
    updateWedgeLabels();
    updateCourseMarkers();
    renderer.render(scene, camera);
  })();
}

if (canvas && sceneData) {
  initGalaxyMap3d(canvas, sceneData);
}
