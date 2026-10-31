// html/static/galaxymap3d.js
// (Why three.js and not another renderer: docs/reference/html-interface.md.)
//
// Renders the interactive 3D Galaxy Map (planetgen/web/maps/galaxymap3d.py) as a real
// WebGL scene (three.js, vendored at ./vendor/three.module.min.js -- see
// vendor/THIRD_PARTY_NOTICES.txt for why it's vendored rather than loaded
// from a CDN). Unlike a sector's scene (every star/cloud in one JSON the
// sector's own route answers, camera always orbiting a fixed origin),
// this camera's own orbit TARGET moves freely through the galaxy -- most
// of what's drawn is fetched live from /galaxy/nodes (the JSON's
// fetchPath) as the camera moves, one sector block ("node") at a time,
// never baked into the page beyond the very first frame
// (#galaxymap3d-data's own "initial" payload) -- see the "Sector-block nodes"
// section below.
//
// Coordinate convention: every position here is the design doc's own
// galaxy-frame parsecs (docs/design/galaxy/galaxy-coordinate-system.md), passed
// straight through as three.js world units with NO axis flip/remap --
// x/y span the galactic plane, +z is galactic north. `camera.up` is set
// to (0, 0, 1) once, at scene setup, specifically so that convention
// reads as "up" on screen without needing any per-position sign flip the
// way this map's own former flat SVG projection (a legacy "+y is down on
// screen" convention) needed one.
//
// The galaxy's sector grid is drawn as blocks of whole sectors, one
// drill-down stage at a time (./galaxystageview.js, rules in
// ./galaxystages.js): the whole galaxy's blocks, then for each block
// clicked the window of blocks one size down centred on it, ... down to a
// window of sectors (MAP.150, MAP.151). Every view can be
// turned and zoomed within limits; the whole galaxy has no grid lines
// drawn on it, so its stars and spiral show through. Density is computed right here from the galaxy's own
// analytic shape, not fetched, and colors each block.
// Unfilled space is mostly see-through; a block holding generated
// sectors is amber and grows more solid with their share (MAP.37).
//
// Nebulae and supernova remnants are drawn over the blocks as soft
// translucent spheres their real size (each node lists the ones reaching
// into it); clicking one shows it and links to its page -- see "Clouds".
//
// Stars are points of light a few pixels across with a soft glow, the
// same size at every zoom: the pre-placed bright stars (500 L☉ or more)
// everywhere, and the stars of generated systems fainter and fainter as
// the view closes in (MAP.51). They can't be clicked (MAP.101): a click
// on one picks what's under it, and a star's details are on its sector's
// page -- see "Stars".

// Sibling modules are imported with this module's own `?v=<version>`
// query (planetgen/web/lib/fmt.py's `static_url`), so they are cached and
// refreshed with the page's script. A plain static `import "./x.js"`
// would drop the query: an update could then leave a stale copy cached,
// and a page that also loaded the same file by its versioned URL would
// get a second, separate instance of it.
const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const {
  cellAt, cellCoordinates, cellVertices, drillBandCount, drillBlockBounds, drillParent, drillWedgeCount, formatDrillKey, isFullTurn,
  neighbourBlock, sectorAddressAt, sectorDesignation,
} = await import(`./galaxyprisms.js${VERSION_QUERY}`);
const { createStageView } = await import(`./galaxystageview.js${VERSION_QUERY}`);
const { createSectorFrame } = await import(`./sectorframe.js${VERSION_QUERY}`);
const { cellEdgePoints, parseSectorDesignation, sectorStage: stageOfSector } = await import(`./galaxystages.js${VERSION_QUERY}`);
const { objectIdOf } = await import(`./findplace.js${VERSION_QUERY}`);
const { createNavPad } = await import(`./mapnavpad.js${VERSION_QUERY}`);
const superbubble = await import(`./superbubble.js${VERSION_QUERY}`);
const { createBlockScene, COLOR_MODES } = await import(`./galaxyblocks.js${VERSION_QUERY}`);
const { formatDistancePc, formatDistanceKm, formatLuminositySol, formatMassSolar } = await import(`./unitladder.js${VERSION_QUERY}`);
const { convert } = await import(`./unitregistry.js${VERSION_QUERY}`);
const {
  BLACK_HOLE_KIND, INTEREST_KINDS, MAP_INTEREST_KINDS, interestLines, MARKER_PX, createMarkerKinds, legendOf: markerLegendOf, mergeMarks, nearFieldRule, selectMarkers,
} = await import(`./galaxymarkers.js${VERSION_QUERY}`);
const { createMarkerLayer } = await import(`./galaxymarkerlayer.js${VERSION_QUERY}`);
const { applyPrefs: applyPoiPrefs, savePrefs: savePoiPrefs } = await import(`./poiprefs.js${VERSION_QUERY}`);
const { bodyMarkers, createHandoff } = await import(`./bodymarkers.js${VERSION_QUERY}`);
const { createGalaxyLabels, LABEL_KINDS } = await import(`./galaxylabels.js${VERSION_QUERY}`);
const { regionLabels } = await import(`./regionlabels.js${VERSION_QUERY}`);
const { formatNumber, threeFigures } = await import(`./numberformat.js${VERSION_QUERY}`);
const { starLight } = await import(`./starlight.js${VERSION_QUERY}`);
const { nebulaLooks } = await import(`./nebulalook.js${VERSION_QUERY}`);
const { drawnSignature, nodeRevision } = await import(`./nodeview.js${VERSION_QUERY}`);
const {
  INSIDE_FACE_ALPHA, NEAR_FIELD_GLSL, nearFieldAtWorld, visibleEnough,
} = await import(`./nearfield.js${VERSION_QUERY}`);
const {
  POINT_FADE, cloudSizeShare, opacitySum, starShare,
} = await import(`./starfade.js${VERSION_QUERY}`);
const {
  GROWN_WHITEN, POINT_PX, POINT_WHITEN, STAR_POINT_GLSL, UNCHARTED_FAR_GLSL, UNCHARTED_FAR_KEEP, massRampShare,
  unchartedKeep, phenomenonLuminositySol, pointBrightness, pointHaloGlow, pointHaloPx,
} = await import(`./starpoint.js${VERSION_QUERY}`);
const { nodeCut, parseNodeKey, settledExpanded, weightCalculator } = await import(`./galaxynodecut.js${VERSION_QUERY}`);
const { SlotStore, relativeOffset, runsOf, splitPair } = await import(`./starbuffer.js${VERSION_QUERY}`);
const { blobsOf } = await import(`./galaxyglow.js${VERSION_QUERY}`);
const { createGlowLayer } = await import(`./galaxyglowlayer.js${VERSION_QUERY}`);
const { regionBlobs, regionLook, chartedText, regionMarkers, inhabitedText } = await import(`./galaxyregions.js${VERSION_QUERY}`);
const { blackHoleDensityPx, labelBufferPx, regionLanding } = await import(`./mappreferences.js${VERSION_QUERY}`);
const {
  regionBlobs: nebulaRegionBlobs, regionLook: nebulaRegionLook, cutToBudget: cutNebulaRegions, FAMILIES: NEBULA_FAMILIES,
  FAMILY_LABELS: NEBULA_FAMILY_LABELS, familyColour: nebulaFamilyColour,
} = await import(`./galaxynebulae.js${VERSION_QUERY}`);
const { createRegionList } = await import(`./galaxyregionlist.js${VERSION_QUERY}`);
const { unpackResponse, unpackNode } = await import(`./tilepack.js${VERSION_QUERY}`);
const {
  STAR_MAGNITUDE_GLSL, TARGET_STARS, TARGET_STARS_PHONE, absoluteMagnitude, apparentMagnitude, tierMagnitudePenalty,
  approachLimit, calibrationLimit, dimmestShown, histogramLimit, limitingMagnitude, magnitudeShare, RAMP,
} = await import(`./starmagnitude.js${VERSION_QUERY}`);
const {
  cssVar, entriesNearScreenPoint, fitRendererToCanvas, formatAddress, isLightBackground, makeRingTexture,
  niceScaleValue, readSceneData, watchResize, worldUnitsPerPixel,
} = await import(`./mapcore.js${VERSION_QUERY}`);
const { watchWindowFit } = await import(`./mapfit.js${VERSION_QUERY}`);
const { extentOf, groupByCell, linePieces } = await import(`./camerarelative.js${VERSION_QUERY}`);
const { blockGenerateButtons, generateButtons } = await import(`./generatebuttons.js${VERSION_QUERY}`);
const {
  NEBULA_MESH_MIN_PX, buildNebulaMesh, createShapeLoader, disposeNebulaMesh,
} = await import(`./nebulamesh.js${VERSION_QUERY}`);
const { createPicker, createRing, createHoverPreview, endpointBookmark, infoPanelOf } = await import(`./mappick.js${VERSION_QUERY}`);
const { createNavPick } = await import(`./navpick.js${VERSION_QUERY}`);

var canvas = document.getElementById("galaxymap3d-canvas");
var dataEl = document.getElementById("galaxymap3d-data");

var sceneData = readSceneData(dataEl);

// The NAV course being picked (navpick.js, NAV.29): the page opened with the
// NAV page's pick, or none; "Start Here" and "End Here" move it on in place.
var pickStart = sceneData || {};
var navPick = createNavPick({
  navUrl: pickStart.navUrl, pick: pickStart.pick, other: pickStart.pickOther, cancelUrl: pickStart.pickCancel,
  nearbyUrl: pickStart.nearbyUrl,
  onChange: function () { if (onPickChanged) onPickChanged(); },
});
var onPickChanged = null; // set once the map is built

// The sector page's URL: the server's template (sceneData.sectorUrl,
// built with page_url so it follows the sector page wherever it lives)
// with the id filled in and the NAV pick carried on.
function sectorUrl(id) {
  return String(sceneData.sectorUrl || "").replace("{id}", encodeURIComponent(id)) + navPick.query();
}

export { formatAddress };

// The map's control buttons (planetgen/web/maps/galaxymap3d.py's panel), by their
// data-action: what each does, given the map's own parts in `ctx`
// ({stageView, setTerritories(on, button), territoriesWanted(),
// setChartedOnly(on)}). A button whose action isn't here does nothing.
export function mapControlHandlers(ctx) {
  return {
    "territories": function (button) {
      ctx.setTerritories(button.getAttribute("aria-pressed") !== "true", button);
      button.setAttribute("aria-pressed", String(ctx.territoriesWanted()));
    },
    "cell-grid": function (button) {
      var on = button.getAttribute("aria-pressed") !== "true";
      button.setAttribute("aria-pressed", String(on));
      ctx.stageView.setCellGrid(on);
    },
    "charted-only": function (button) {
      var on = button.getAttribute("aria-pressed") !== "true";
      button.setAttribute("aria-pressed", String(on));
      ctx.setChartedOnly(on);
    },
    "back": function () { ctx.stageView.travel(-1); },
    "forward": function () { ctx.stageView.travel(1); },
    "current": function () { ctx.stageView.travelToCurrent(); },
    "up": function () { ctx.stageView.up(); },
    "reset": function () { ctx.stageView.home(); },
    "reset-view": function () { ctx.stageView.resetView(); },
    "step": function (button) { ctx.stageView.stepToNeighbour(button.dataset.step); },
    "map-help": function (button) {
      var dialog = document.getElementById("galaxymap3d-help");
      var menu = button && button.closest("details");
      if (menu) menu.open = false;
      if (dialog && typeof dialog.show === "function") dialog.show();
    },
    "center-on": function () { ctx.stageView.recenter(); },
    "zoom-in": function () { ctx.stageView.zoomIn(); },
    "zoom-out": function () { ctx.stageView.zoomOut(); },
  };
}

// The controls' Menu (MAP.55, a <details> dropping down over the page):
// Escape closes it and puts the focus back on its button, and a press
// anywhere outside it closes it.
export function wireMapMenu(menuEl) {
  menuEl.addEventListener("keydown", function (event) {
    if (event.key === "Escape" && menuEl.open) {
      event.preventDefault();
      menuEl.open = false;
      menuEl.querySelector("summary").focus();
    }
  });
  menuEl.ownerDocument.addEventListener("pointerdown", function (event) {
    if (menuEl.open && !menuEl.contains(event.target)) {
      menuEl.open = false;
    }
  });
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
//
// The panel is mappick.js's, shared with the Sector Map (MAP.65); these
// say what goes in it.

function infoPanel() {
  return infoPanelOf(document.getElementById("galaxymap3d-info"));
}

// What the panel says before anything is picked (a sector page's, which
// the panel returns to when a pick is cleared; MAP.113).
var restingHint = (infoPanel() && infoPanel().element && infoPanel().element.textContent.trim()) || "";

// The largest block whose info names a region: a bigger one spans many regions.
var REGION_BLOCK_MAX = 9;

// MAP.315: the charted region a point lies in, set once the map has its nodes: regionLookup(x, y, z) resolves to
// {name, parent, block, extent, parentBlock} (the 9-block with its extent {centre, radius} in pc, and its 81-block family,
// with their names) or null, and regionOpener.region(region) frames the region's extent (a block's window when the node
// has none) and regionOpener.block(block) flies the camera to a block.
var regionLookup = null;
var regionOpener = null;

// MAP.349: a Go to button for a line of a card that names a container or a place: it flies the camera there as the find box, the crumbs
// and the results list do (the map's one goTo), and is named for what it goes to.
function goToButton(name, go) {
  var button = document.createElement("button");
  button.type = "button";
  button.className = "btn btn-small btn-secondary map-info-go";
  button.textContent = "Go to";
  button.setAttribute("aria-label", "Go to " + name);
  button.addEventListener("click", go);
  return button;
}

// A line of a card, "<label> <name>" with its Go to button.
function placeLine(label, name, go) {
  var line = document.createElement("p");
  line.className = "map-info-place";
  line.append(label + " ", name, " ", goToButton(name, go));
  return line;
}

// MAP.315, MAP.349: the lines "Region: <name>" and "In: <parent>" for the point (x, y, z) pc, each with its own Go to that flies to its
// block (the region framed by its extent, MAP.319; its family by its 81-block). They stay hidden until the region is known, and for a
// point in no charted region (a 9-block with no generated sector has no name on the wire).
function regionNode(position) {
  if (!regionLookup || !position) {
    return null;
  }
  var box = document.createElement("div");
  box.className = "map-info-region";
  box.hidden = true;
  regionLookup(position[0], position[1], position[2]).then(function (region) {
    if (!region || !region.name) {
      return;
    }
    box.appendChild(placeLine("Region:", region.name, function () { regionOpener.region(region); }));
    if (region.parent) {
      box.appendChild(placeLine("In:", region.parent, function () { regionOpener.block(region.parentBlock); }));
    }
    box.hidden = false;
    if (box.isConnected && infoPanel()) {
      infoPanel().refit();
    }
  }, function () {
    // A region that could not be looked up is a line left out.
  });
  return box;
}

// `spec` with the region line (MAP.315) under its rows, for the point `position` ([x, y, z] pc).
function withRegion(spec, position) {
  var line = regionNode(position);
  if (line) {
    spec.after = [line].concat(spec.after || []);
  }
  return spec;
}

function showInfo(spec) {
  var panel = infoPanel();
  if (panel) {
    panel.show(spec);
  }
}

// A generated sector's ☆ entry: its designation, opening its page.
function sectorBookmark(entry) {
  if (!entry.designation) {
    return null;
  }
  return {
    kind: "sector", value: entry.designation, name: entry.name || "Sector " + entry.designation,
    url: entry.id != null && sceneData.sectorUrl ? sectorUrl(entry.id).split("?")[0] : null,
    sectorId: entry.id != null ? entry.id : null,
  };
}

function showPlacedInfo(entry) {
  var relativeDensity = placedRelativeDensity(entry, sceneData.referenceDensityPerLy3);
  showInfo(withRegion({
    title: entry.name || "Unnamed sector",
    kind: "Sector",
    fields: [
      ["Systems", entry.system_count != null ? entry.system_count : 0],
      ["Density", relativeDensity != null ? relativeDensity.toFixed(2) + "× local average" : null],
      ["Distance from core", entry.galactic_radius_pc != null ? Math.round(entry.galactic_radius_pc) + " pc" : null],
      ["Designation", entry.designation],
    ],
    more: [["Address", entry.shell_index != null ? formatAddress(entry.shell_index, entry.band_index, entry.wedge_index) : null]],
    bookmark: sectorBookmark(entry),
    view: sceneData.sectorUrl && entry.id != null ? { href: sectorUrl(entry.id), label: "View sector →" } : null,
    links: entry.id != null ? navPick.nearbyLinks("sector:" + entry.id) : [],
  }, entry.position));
}

// A share as a percentage, never rounding a nonzero one down to 0%.
function formatShare(share) {
  if (share > 0 && share < 0.001) {
    return "< 0.1%";
  }
  return (100 * share).toFixed(share < 0.1 ? 1 : 0) + "%";
}

export { sectorDesignation };

// Degrees counterclockwise from +X, 0 to 360 (the longitude a cell's side is at).
function longitudeDeg(rad) {
  return (((rad * 180) / Math.PI) % 360 + 360) % 360;
}

// The fields every cell's panel gives for its place: where its middle is, how far from the core, and the
// ranges its edges span (radii, latitudes north of the plane, longitudes).
function placeFields(b) {
  return [["Distance from core", formatDistancePc(cellCoordinates(b).spherical[0])]];
}

// The rest of a cell's place, behind the card's More (MAP.325): the middle, the depth and the edges' ranges.
function placeMoreFields(b) {
  var c = cellCoordinates(b).cartesian;
  var fields = [];
  fields.push(["Center x, y, z", c.map(function (v) { return v.toFixed(1); }).join(", ") + " pc"]);
  fields.push(["Shell depth", formatDistancePc(b.r1 - b.r0)]);
  fields.push(["Latitude", (90 - (b.c1 * 180) / Math.PI).toFixed(2) + "° to " + (90 - (b.c0 * 180) / Math.PI).toFixed(2) + "°"]);
  fields.push(["Longitude", isFullTurn(b) ? "all the way round"
    : longitudeDeg(b.t0).toFixed(2) + "° to " + (longitudeDeg(b.t1) || 360).toFixed(2) + "°"]);
  return fields;
}

// A sector cell's info (a single sector). `cell` is {m: 1, bounds, address:
// {shell, band, wedge}}, plus `filled` (generated sectors in it) when it
// came from the drawn blocks. A generated sector's own panel is
// showPlacedInfo.
function showCellInfo(cell) {
  var b = cell.bounds;
  var a = cell.address;
  var fields = [];
  fields.push(["Designation", sectorDesignation(a.shell, a.band, a.wedge)]);
  if (cell.filled != null) {
    fields.push(["Charted", cell.filled > 0 ? "Yes" : "Not yet"]);
  }
  if (cell.uncharted) {
    // Opened in place with what the scatters left in it (MAP.162).
    fields.unshift(["Status", "Uncharted: not charted yet, but scattered objects are placed here"]);
    fields.push(["Stars waiting", formatNumber(cell.uncharted.stars)]);
    fields.push(["Other objects", formatNumber(cell.uncharted.objects)]);
  }
  if (cell.density != null) {
    fields.push(["Predicted density", cell.density.toFixed(2) + "× local average"]);
  }
  fields = fields.concat(placeFields(b));
  var more = [["Address", formatAddress(a.shell, a.band, a.wedge)]].concat(placeMoreFields(b));

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

  // A sector that isn't generated yet: for a logged-in admin
  // (sceneData.generate, set by planetgen/web/maps/galaxymap3d.py only then), the same
  // Generate buttons the Sector Map gives a neighbor.
  var generate = !(cell.filled > 0) && sceneData.generate
    ? generateButtons(sceneData.generate, a.shell, a.band, a.wedge, sceneData.edgeLy)
    : null;
  showInfo(withRegion({
    title: cell.uncharted ? "Uncharted sector" : "Sector cell",
    fields: fields, more: more, after: isFullTurn(b) ? [] : [corners], generate: generate,
    // NAV.45: the cell's middle, so a sector nothing is generated in can still be asked about.
    links: navPick.nearbyLinks(cellCoordinates(b).cartesian),
  }, cellCoordinates(b).cartesian));
}

// A drill-down block's info (galaxystageview.js): its sectors, how many are
// charted and its distance from the core up front; its shells, row, wedge
// and where it is under More; and (`enter`) a button that flies into it.
function showBlockInfo(info, edgePc) {
  var block = info.block;
  var b = drillBlockBounds(block, edgePc);
  var first = block.shell * block.m;
  var fields = [];
  // MAP.325: the card says how big it is, how much is charted and how far it is; where it sits in the shells is behind More.
  var more = [["Shells", first + "–" + (first + block.m - 1)]];
  var rows = drillBandCount(block.m, block.shell);
  more.push(["Row", block.band + " of " + rows + " in its shell (north +)"]);
  more.push(["Wedge", block.wedge + " of " + drillWedgeCount(block.m, block.shell, block.band) + " in its row"]);
  if (info.total != null) {
    fields.push(["Sectors", formatNumber(info.total)]);
  }
  if (info.generated != null) {
    fields.push(["Charted", formatNumber(info.generated)
      + (info.generated > 0 && info.total > 0 ? " (" + formatShare(info.generated / info.total) + ")" : "")]);
  }
  fields = fields.concat(placeFields(b));
  more = more.concat(placeMoreFields(b));
  showInfo(withRegion({
    title: "Block " + block.shell + "·" + block.band + "·" + block.wedge + " (" + block.m + " sectors a side)",
    fields: fields,
    more: more,
    buttons: info.enter ? [{ label: "Fly into this block →", onClick: info.enter }] : [],
    // UX.127: an admin's Generate for the block, the same job the command's --block runs.
    generate: sceneData.generate ? blockGenerateButtons(sceneData.generate, formatDrillKey(block), info.total) : null,
  }, block.m <= REGION_BLOCK_MAX ? cellCoordinates(b).cartesian : null));
}

// A hint paragraph in the info panel: replacing what's there, or (keep)
// under it.
function showHint(text, keep) {
  var panel = infoPanel();
  if (panel) {
    panel.hint(text, keep);
  }
}

// --- Clouds: nebulae and supernova remnants --------------------------------

// A nebula's color by its type, per theme (static/nebulalook.js, MAP.175): the
// same hues the Sector Map uses (planetgen/web/maps/starmap.py's
// _NEBULA_TYPE_COLORS and _NEBULA_TYPE_ALPHA, as core and edge opacity),
// each reaching 3:1 against the map's background on its theme. A dark
// nebula is a dust-coloured silhouette on the dark theme and a near-black
// one on the light.
var NEBULA_LOOKS = nebulaLooks(isLightBackground());
var DEFAULT_NEBULA_LOOK = NEBULA_LOOKS.default;

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

// A phenomenon's ☆ entry (its NAV endpoint) and links: its page, and
// "Start Here" and "End Here", or the one button for the end being picked.
function phenomenonActions(item) {
  var endpoint = item.type + ":" + item.id;
  var spec = { bookmark: endpointBookmark(endpoint, item.name, sceneData.phenomenonUrl ? phenomenonUrl(item) : null) };
  // Not while choosing a NAV endpoint (NAV.30): it would leave the course.
  if (sceneData.phenomenonUrl && !navPick.active()) {
    // A scattered point that nothing has built has no record of its own: its page is the scatter row's
    // (MAP.211), which sends a visitor on to the object once its sector is generated.
    var href = item.scattered && sceneData.unchartedPhenomenonUrl
      ? String(sceneData.unchartedPhenomenonUrl).replace("{id}", encodeURIComponent(String(item.id).replace(/^s/, "")))
      : phenomenonUrl(item);
    spec.view = { href: href, label: "View phenomenon →" };
  }
  if (sceneData.navUrl) {
    spec.nav = navPick.actionsFor(endpoint, item.name);
  }
  // NAV.45: not for a scattered point nothing has built yet, which has no record to ask from.
  if (!item.scattered) {
    spec.links = (spec.links || []).concat(navPick.nearbyLinks(endpoint));
  }
  return spec;
}

// MAP.316: a charted region picked on the map: how much of it is charted, its family, and a button that flies to it.
// MAP.348: what a card says first of each kind of thing (its name, kind and values), also what the hover preview shows of it.
function regionFacts(blob) {
  return {
    title: blob.name || "Charted region",
    kind: "Region",
    fields: [
      ["Charted", chartedText(blob.share)],
      ["Sectors generated", blob.sectorCount + " of " + blob.capacity],
      ["Inhabited", blob.inhabited ? inhabitedText(blob.inhabited) : null],
      ["Family", blob.parent],
      ["Distance from core", formatDistancePc(Math.hypot(blob.x, blob.y, blob.z))],
    ],
  };
}

function showRegionInfo(blob) {
  showInfo(withRegion(Object.assign(regionFacts(blob), {
    more: [["Radius", formatDistancePc(blob.radiusPc)], ["Region cell", blob.key]],
    buttons: [{ label: "Go to region", onClick: function () { openRegionBlob(blob); } }],
  }), [blob.x, blob.y, blob.z]));
}

// Flies to a charted region of the regions layer: its window, the camera framing its extent (MAP.319).
function openRegionBlob(blob) {
  if (!regionOpener) return;
  regionOpener.region({ block: blob.block, name: blob.name, extent: { centre: [blob.x, blob.y, blob.z], radius: blob.radiusPc } });
}

function cloudFacts(cloud) {
  return {
    title: cloud.name || cloudTypeLabel(cloud),
    fields: [
      ["Type", cloudTypeLabel(cloud)],
      ["Class", cloud.class],
      ["Radius", formatDistancePc(cloud.radius_pc)],
      ["Distance from core", formatDistancePc(Math.hypot(cloud.x, cloud.y, cloud.z))],
    ],
  };
}

function showCloudInfo(cloud) {
  showInfo(withRegion(Object.assign(cloudFacts(cloud), {
    more: [["Center x, y, z", [cloud.x, cloud.y, cloud.z].map(function (v) { return v.toFixed(1); }).join(", ") + " pc"]],
  }, phenomenonActions(cloud)), [cloud.x, cloud.y, cloud.z]));
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

var POINT_TYPE_LABELS = { black_hole: "Black Hole", neutron_star: "Neutron Star", quasar: "Quasar" };

// "Black Hole (Accreting)", "Neutron Star (Millisecond pulsar)".
function pointTypeLabel(point) {
  var label = POINT_TYPE_LABELS[point.type] || capitalize(String(point.type).replace(/_/g, " "));
  if (!point.descriptor) {
    return label;
  }
  if (point.type === "neutron_star") {
    return label + " (" + (point.descriptor === "non-pulsing" ? "Not pulsing" : capitalize(point.descriptor) + " pulsar") + ")";
  }
  return label + " (" + capitalize(point.descriptor) + ")";
}

// A black hole, neutron star or quasar on the map (MAP.80,
// queryDb.galaxy_point_phenomena_in_box): what it is, and its page.
function pointFacts(point) {
  return {
    title: point.name || POINT_TYPE_LABELS[point.type] || "Phenomenon",
    fields: [
      ["Type", pointTypeLabel(point)],
      ["Luminosity", formatLuminositySol(point.luminosity_sol)],
    ].concat(point.mass_solar != null ? [["Mass", formatMassSolar(point.mass_solar)]] : [], [
      ["Distance from core", formatDistancePc(Math.hypot(point.x, point.y, point.z))],
    ]),
  };
}

function showPointInfo(point) {
  showInfo(withRegion(Object.assign(pointFacts(point), {
    more: [["Position x, y, z", [point.x, point.y, point.z].map(function (v) { return v.toFixed(1); }).join(", ") + " pc"]],
  }, phenomenonActions(point)), [point.x, point.y, point.z]));
}

// MAP.296: the nearest systems of a star system, under its position in the card: each name a link to that system's
// page, nearest first, with its distance. The list is fetched when the card opens (/galaxy/system/{id}/nearest) and
// the card fits itself again when it arrives.
function nearestSystemsNode(systemId) {
  if (!sceneData.systemNearestPath || systemId == null) {
    return null;
  }
  var box = document.createElement("div");
  box.className = "map-info-nearest";
  var heading = document.createElement("h4");
  heading.textContent = "Nearest systems";
  var note = document.createElement("p");
  note.className = "hint";
  note.textContent = "Looking up the nearest systems…";
  box.append(heading, note);
  function settle(message, list) {
    note.textContent = message || "";
    note.hidden = !message;
    if (list) {
      box.appendChild(list);
    }
    if (box.isConnected && infoPanel()) {
      infoPanel().refit();
    }
  }
  fetch(String(sceneData.systemNearestPath).replace("{id}", encodeURIComponent(systemId)), {
    headers: { Accept: "application/json" }, credentials: "same-origin",
  }).then(function (response) {
    return response.ok ? response.json() : Promise.reject(new Error("status " + response.status));
  }).then(function (body) {
    var near = body.nearest || [];
    if (!near.length) {
      settle("No other system is stored near it.");
      return;
    }
    var list = document.createElement("ol");
    list.className = "map-info-nearest-list";
    near.forEach(function (system) {
      var item = document.createElement("li");
      var link = document.createElement("a");
      link.href = system.url;
      link.textContent = system.name || "Unnamed system";
      var away = document.createElement("span");
      away.className = "map-info-nearest-distance";
      away.textContent = " " + threeFigures(system.distance_ly) + " ly";
      item.append(link, away);
      // MAP.349: each listed system has its Go to, which opens its sector in place and the system itself, as the find box does.
      var found = objectIdOf(system.id, parseSectorDesignation);
      if (found && regionOpener) {
        item.append(" ", goToButton(system.name || "the system", function () {
          regionOpener.place({ sector: found.sector, open: true, system: found.system, label: system.name || null });
        }));
      }
      list.appendChild(item);
    });
    settle("", list);
  }).catch(function () {
    settle("The nearest systems could not be loaded.");
  });
  return box;
}

// A star system worth a mark on the map (MAP.220, MAP.305: a habitable world, a gas giant): every reason it has, where it is,
// its nearest systems (MAP.296), and its page.
function interestFacts(system) {
  return {
    title: system.name || "Star system",
    kind: "Star system",
    fields: interestLines(system).map(function (line) { return [line.label, line.text]; }).concat([
      ["Distance from core", formatDistancePc(Math.hypot(system.x, system.y, system.z))],
    ]),
  };
}

function showInterestSystemInfo(system) {
  var spec = Object.assign(interestFacts(system), {
    more: [["Position x, y, z", [system.x, system.y, system.z].map(function (v) { return v.toFixed(1); }).join(", ") + " pc"]],
  });
  var nearest = nearestSystemsNode(system.id);
  if (nearest) {
    spec.after = [nearest];
  }
  if (sceneData.systemUrl) {
    spec.view = { href: String(sceneData.systemUrl).replace("{id}", encodeURIComponent(system.id)), label: "View system →" };
  }
  showInfo(withRegion(spec, [system.x, system.y, system.z]));
}

// MAP.243: a star beyond the galaxy's outline (the halo, a satellite galaxy), where no block or sector of the
// whole-galaxy view holds it: what is known of it and where it is.
// A star of the map (MAP.304): what the node lists of it, and its system's page when it has one.
function starFacts(star) {
  return {
    title: star.name || "Star",
    kind: "Star",
    fields: [
      ["Luminosity", formatLuminositySol(star.luminosity_sol)],
      ["Temperature", star.temperature_k != null ? formatNumber(star.temperature_k) + " K" : null],
      ["Class", star.star_type || null],
      ["Distance from core", formatDistancePc(Math.hypot(star.x, star.y, star.z))],
    ],
  };
}

function showStarInfo(star, nearby) {
  var spec = Object.assign(starFacts(star), {
    more: [["Position x, y, z", [star.x, star.y, star.z].map(function (v) { return v.toFixed(1); }).join(", ") + " pc"]],
  });
  if (star.system_id != null && sceneData.systemUrl) {
    var systemPage = String(sceneData.systemUrl).replace("{id}", encodeURIComponent(star.system_id));
    spec.view = { href: systemPage, label: "View system →" };
    // MAP.122: the ☆ saves the star's system (a star in unmade space has no system to keep).
    spec.bookmark = endpointBookmark("system:" + star.system_id, star.name || "Star system", systemPage);
  }
  if (nearby) {
    spec.after = [nearby];
  }
  showInfo(withRegion(spec, [star.x, star.y, star.z]));
}

function haloStarFacts(star) {
  return {
    title: "Star beyond the galaxy's outline",
    fields: [
      ["Luminosity", formatLuminositySol(star.luminosity_sol)],
      ["Temperature", star.temperature_k != null ? formatNumber(star.temperature_k) + " K" : null],
      ["Class", star.star_type || null],
      ["Distance from core", formatDistancePc(Math.hypot(star.x, star.y, star.z))],
    ],
  };
}

function showHaloStarInfo(star) {
  showInfo(Object.assign(haloStarFacts(star), {
    more: [["Position x, y, z", [star.x, star.y, star.z].map(function (v) { return v.toFixed(1); }).join(", ") + " pc"]],
  }));
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

// --- Filled-sector colors ---------------------------------------------
//
// Unlike planetgen/web/maps/starmap.py (size/color baked server-side, once), these run
// client-side for EVERY fetch (initial payload and every live re-fetch
// alike) -- see galaxymap3d.py's own module docstring for why: almost
// everything drawn here arrives through a live fetch that never passes
// through that Python module again, so a server-computed style would
// only ever apply to the first frame.

// A placed sector's own REAL stellar density (system_count / edge_ly^3),
// relative to physical_constants.LOCAL_STELLAR_DENSITY_LY3 (the real
// local-neighborhood average this whole generator already calibrates
// against -- see planetgen/web/maps/galaxymap3d.py's own referenceDensityPerLy3
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
  var BASE_NEAR = Math.max(data.minViewRadiusPc / 50, 0.001);
  var camera = new THREE.PerspectiveCamera(FOV_DEG, 1, BASE_NEAR, farPlane);
  camera.up.set(0, 0, 1);

  var MIN_RADIUS = data.minViewRadiusPc;
  var MAX_RADIUS = data.maxViewRadiusPc;
  var GALAXY_RADIUS = data.galaxyRadiusPc || data.maxViewRadiusPc;

  var target = new THREE.Vector3(data.initialCenter[0], data.initialCenter[1], data.initialCenter[2]);
  var initialRadius = data.initialRadiusPc;

  // theta: azimuth from +x in the xy-plane; phi: polar angle from +z --
  // the same (r, theta, phi) convention docs/design/galaxy/galaxy-coordinate-
  // system.md and planetgen.galaxy.geometry.sector_position_pc use,
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

  // MAP.338: the camera as a frame plus an offset. `camera.position` is a double, but at 8 kpc one step of a double is 1.8e-12 pc (55 km),
  // as big as the nearest view (1e-12 pc), so the position the camera draws from would jitter. `eyeBack` is the camera's own offset from
  // `target`, exact however far the target is; once the camera is nearer its target than DEEP_FRAME_PC and a system is open in place, the
  // system is drawn in a second pass from the eye's own frame (the eye at the origin, the system at its offset from it, all scaled so the
  // view is about a unit across: the same picture, with the log depth buffer's resolution where the bodies are).
  var eyeBack = new THREE.Vector3();
  var DEEP_FRAME_PC = 1e-6;
  var deepFrame = { active: false, scale: 1, passes: 0 };
  var DEEP_CLEAR = new THREE.Vector3();

  function renderDeepFrame(frame) {
    var dist = eyeBack.length();
    var scale = Math.pow(10, Math.round(-Math.log10(Math.max(dist, 1e-15))));
    var group = frame.group;
    var saved = { position: group.position.clone(), scale: group.scale.clone(), visible: group.visible,
      cameraPosition: camera.position.clone(), near: camera.near, far: camera.far };
    var hidden = [];
    scene.children.forEach(function (child) {
      if (child !== group && child.visible && !child.isLight) { child.visible = false; hidden.push(child); }
    });
    // The offset of the system from the eye, worked out as (centre - target) - back so no step of a double at the galaxy's size is in it.
    group.position.set((frame.centre[0] - target.x - eyeBack.x) * scale, (frame.centre[1] - target.y - eyeBack.y) * scale,
      (frame.centre[2] - target.z - eyeBack.z) * scale);
    group.scale.multiplyScalar(scale);
    camera.position.copy(DEEP_CLEAR);
    camera.near = Math.min(BASE_NEAR, dist * 1e-3) * scale;
    camera.far = Math.max(dist * 1e6, 1) * scale;
    camera.updateProjectionMatrix();
    camera.updateMatrixWorld();
    group.updateMatrixWorld(true);
    renderer.autoClear = false;
    renderer.clearDepth();
    renderer.render(scene, camera);
    renderer.autoClear = true;
    hidden.forEach(function (child) { child.visible = true; });
    group.position.copy(saved.position);
    group.scale.copy(saved.scale);
    group.visible = saved.visible;
    camera.position.copy(saved.cameraPosition);
    camera.near = saved.near;
    camera.far = saved.far;
    camera.updateProjectionMatrix();
    camera.updateMatrixWorld();
    group.updateMatrixWorld(true);
    deepFrame.scale = scale;
  }

  function renderFrame() {
    var frame = stageView.systemFrame();
    var deep = !!frame && frame.group.visible && eyeBack.length() > 0 && eyeBack.length() < DEEP_FRAME_PC;
    deepFrame.active = deep;
    if (!deep) {
      renderer.render(scene, camera);
      return;
    }
    // The map without the system first, then the system over it from its own frame.
    frame.group.visible = false;
    renderer.render(scene, camera);
    frame.group.visible = true;
    renderDeepFrame(frame);
    deepFrame.passes++;
  }

  // Nothing of the galaxy is clipped to the window the drill-down shows (MAP.264): the stars and clouds all stay
  // when you zoom, and the visibility law decides what is drawn. There are no wedge lines: the galaxy shows its stars
  // and spiral structure with no grid drawn over it (MAP.85).

  // A sector opened in place (galaxysector.js) draws its own stars and
  // phenomena, so the galaxy's are left out of its cell (bounds as a
  // block has them); null brings them back.
  var hiddenCell = null;

  function hideCell(bounds) {
    hiddenCell = bounds;
    markHiddenStars();
    updateClouds();
    applySectorView();
  }

  // The middle of a cell given as bounds (radius, colatitude, longitude), in the galaxy's frame.
  function cellCentre(b) {
    var r = (b.r0 + b.r1) / 2;
    var c = (b.c0 + b.c1) / 2;
    var t = (b.t0 + b.t1) / 2;
    return [r * Math.sin(c) * Math.cos(t), r * Math.sin(c) * Math.sin(t), r * Math.cos(c)];
  }

  // MAP.327: with a sector opened in place, what lies outside it is dimmed to SECTOR_OUTSIDE_SHARE (no glow, no names, no
  // comets or rogue planets), the sector is framed by a wireframe that stays where the cell is, and each wall that faces a
  // neighbouring sector carries an icon that flies there. With none opened the map is as it was.
  var SECTOR_OUTSIDE_SHARE = 0.3;
  var sectorFrame = null;
  var regionFrame = null;
  var framedRegion = "";
  var framedAt = "";
  var glowLayerOn = true;
  function applySectorView() {
    starMaterial.uniforms.outsideDim.value = hiddenCell ? SECTOR_OUTSIDE_SHARE : 1;
    if (glow) glow.setVisible(glowLayerOn && !hiddenCell);
    if (regionLayer) regionLayer.setVisible(regionsOn && !hiddenCell);
    if (nebulaLayer) nebulaLayer.setVisible(!hiddenCell);
    if (markerLayer) markerLayer.setOutsideShare(hiddenCell ? SECTOR_OUTSIDE_SHARE : 1);    // MAP.341: the marks dim as the stars do
    if (sectorFrame) sectorFrame.set(hiddenCell);
  }

  // How far round from bearing a0 bearing `angle` is, 0 to a full turn.
  function bearingFrom(angle, a0) {
    return (((angle - a0) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  }

  // Whether a point is drawn by the galaxy's own layers: not inside the cell of the sector opened in place, which
  // draws its own.
  function drawnByTheGalaxy(x, y, z) {
    if (!hiddenCell) {
      return true;
    }
    var plane = Math.hypot(x, y);
    var r = Math.hypot(plane, z);
    var colat = Math.atan2(plane, z);
    var angle = Math.atan2(y, x);
    var inside = function (b) {
      if (r < b.r0 || r > b.r1 || colat < b.c0 || colat > b.c1) {
        return false;
      }
      return isFullTurn(b) || plane === 0 || bearingFrom(angle, b.t0) <= b.t1 - b.t0;
    };
    return !inside(hiddenCell);
  }

  // --- The block solid -------------------------------------------------------
  //
  // A stage's blocks (galaxyprisms.js) fill their whole cells; a thin
  // brighter outline along each face's own edges (found from the face's
  // 0..1 uv, about a screen pixel wide) keeps neighboring blocks apart.
  // Each part is two meshes sharing one shader:
  // - solid: blocks opaque enough to write depth (none now: no fill is
  //   more than 50% opaque, galaxyblocks.MAX_FILL_OPACITY);
  // - glass: everything else, translucent (no depth writes), its blocks
  //   sorted back to front from the camera when built. Unfilled space is
  //   10% to 30% opaque by density, and a filled sector more opaque by what
  //   it holds, up to 50% (galaxyblocks.statsOpacity, MAP.128).
  // `fade` scales a mesh's opacity (the drill-down's fades and dimming);
  // `gridEdges` (0 or 1) turns the outlines off, as on the whole galaxy,
  // which shows no sector or block lines (MAP.85).
  // The logdepthbuf chunks match the renderer's logarithmic depth buffer.
  // `faceAlpha` 0 (a sector page's map, MAP.61) draws only the outlines:
  // sector space is not tinted at all.
  // Filled blocks' edges are a saturated amber (MAP.37: "a much higher
  // contrast"), a deeper shade on the light theme's pale background.
  // Their faces take the color of what their sectors hold
  // (galaxyblocks.statsColor, MAP.128), so no amber over them.
  var FILLED_TINT = new THREE.Color(isLightBackground() ? "#d06a00" : "#ffb02e");
  var FILLED_FACE_MIX = 0;
  // How far an unfilled block's edges brighten toward white: the ring,
  // wedge and layer boundaries of the grid, kept faint like the wedge
  // lines (a filled block's amber edges add up to 0.7 more).
  var GRID_EDGE_MIX = 0.1;

  // MAP.149, the near field (static/nearfield.js): the camera-to-focus
  // distance, the focus's radius (0: unknown) and the center of the block the
  // camera is inside (far away: none), shared by every block and star
  // material and set with the camera (setCamera).
  var NO_BLOCK_INSIDE = 1e12;
  // MAP.188: the share every block face and line is drawn at so the
  // blocks along a ray add up to no more than GRID_ACCUMULATION_CAP
  // (static/gridopacity.js), set by the stage view with the camera.
  var gridUniform = { value: 1 };
  var nearUniforms = {
    nearD: { value: 1 }, nearR: { value: 0 },
    insideCenter: { value: new THREE.Vector3(NO_BLOCK_INSIDE, NO_BLOCK_INSIDE, NO_BLOCK_INSIDE) },
  };

  // The near field's share (0..1) at a world point, as the shaders work it
  // out (the picker uses the same numbers, so picking agrees with drawing).
  function nearFieldAt(x, y, z) {
    return nearFieldAtWorld(camera.matrixWorldInverse.elements, x, y, z, nearUniforms.nearD.value, nearUniforms.nearR.value);
  }

  function makeBlockMaterial(translucent) {
    return new THREE.ShaderMaterial({
      uniforms: {
        filledTint: { value: FILLED_TINT }, fade: { value: 1 }, gridEdges: { value: 1 },
        faceAlpha: { value: data.pinned ? 0 : 1 },
        nearD: nearUniforms.nearD, nearR: nearUniforms.nearR, insideCenter: nearUniforms.insideCenter,
        gridScale: gridUniform,
      },
      transparent: translucent,
      depthWrite: !translucent,
      // Both sides drawn: the back faces only of the block the camera is
      // inside (the fragment shader drops the others), so a container is
      // seen from the inside (MAP.149). Picking ignores back faces.
      side: THREE.DoubleSide,
      vertexShader: [
        "#include <common>",
        "#include <logdepthbuf_pars_vertex>",
        "attribute vec3 prismColor;",
        "attribute float prismAlpha;",
        "attribute float prismFill;",
        "attribute vec2 faceUv;",
        "attribute vec3 prismCenter;",
        "varying vec3 vColor;",
        "varying float vAlpha;",
        "varying float vFill;",
        "varying vec2 vUv;",
        "varying vec3 vView;",
        "varying vec3 vCenter;",
        "void main() {",
        "  vColor = prismColor;",
        "  vAlpha = prismAlpha;",
        "  vFill = prismFill;",
        "  vUv = faceUv;",
        "  vCenter = prismCenter;",
        "  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);",
        "  vView = mvPosition.xyz;",
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
        "uniform float gridEdges;",
        "uniform float faceAlpha;",
        "uniform float nearD;",
        "uniform float nearR;",
        "uniform vec3 insideCenter;",
        "uniform float gridScale;",
        "varying float vAlpha;",
        "varying float vFill;",
        "varying vec2 vUv;",
        "varying vec3 vView;",
        "varying vec3 vCenter;",
        NEAR_FIELD_GLSL,
        "void main() {",
        "  #include <logdepthbuf_fragment>",
        // Back faces are only the walls of the block the camera is in.
        "  bool own = distance(vCenter, insideCenter) < 0.05;",
        "  if (!gl_FrontFacing && !own) discard;",
        "  vec2 toEdge = min(vUv, 1.0 - vUv) / max(fwidth(vUv), vec2(1e-6));",
        "  float edge = gridEdges * (1.0 - smoothstep(0.5, 1.5, min(toEdge.x, toEdge.y)));",
        "  vec3 face = mix(vColor, filledTint, step(0.001, vFill) * " + FILLED_FACE_MIX.toFixed(3) + ");",
        "  vec3 edgeColor = mix(vec3(1.0), filledTint, step(0.001, vFill));",
        "  float shown = mix(edge, 1.0, faceAlpha);",
        "  if (!gl_FrontFacing) shown = mix(edge, 1.0, " + INSIDE_FACE_ALPHA.toFixed(3) + ");",
        "  gl_FragColor = vec4(mix(face, edgeColor, (" + GRID_EDGE_MIX.toFixed(3) + " + 0.7 * vFill) * edge), vAlpha * fade * shown * gridScale * nearField(vView, nearD, nearR));",
        "  #include <colorspace_fragment>",
        "}",
      ].join("\n"),
    });
  }

  // The selection ring and, fainter, the hover ring (mappick.js), each
  // HIGHLIGHT_PX across on screen whatever its distance, drawn over the
  // blocks.
  var selectionRing = createRing(scene, camera, canvasEl, accentColor, { depthTest: false, renderOrder: 5 });
  var hoverRing = createRing(scene, camera, canvasEl, accentColor, { depthTest: false, renderOrder: 5, opacity: 0.45 });

  function highlightPosition(x, y, z) {
    selectionRing.at(x, y, z, { px: HIGHLIGHT_PX * 2 });
  }

  // Density shading's dim-to-accent-to-white ramp (galaxyblocks.js).
  var PRISM_DIM = new THREE.Color(0x1d2340);
  var PRISM_ACCENT = new THREE.Color(accentColor);
  var PRISM_HOT = new THREE.Color(0xeef0ff);
  var galaxyShape = data.densityShape || null;
  var edgePc = data.edgePc || 1;

  // Parsecs per screen pixel at the view's focus (the orbit target), for
  // the camera at orbit radius `radius` (default: where it is). Per CSS
  // pixel, on purpose, not per device pixel: the block size (a block at
  // least blockMinPx across) and the scale line are about what
  // a person can see and click, which is the same on a 2x phone screen as
  // on a 1x monitor. The renderer draws at up to 2 device pixels per CSS
  // pixel (setPixelRatio) for sharpness only, and the bright stars' sizes
  // are CSS pixels too (scaled by pixelRatio in their shader).
  function pcPerPixelAtTarget(radius) {
    return worldUnitsPerPixel(camera, radius || orbit.radius, canvasEl.clientHeight);
  }

  // --- Blocks ------------------------------------------------------------------
  //
  // The drill-down (galaxystageview.js) draws each stage's blocks with
  // galaxyblocks.js's scene and the block shader above.
  var reducedMotion = typeof window.matchMedia === "function"
    && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  var blockSceneConfig = {
    shape: galaxyShape,
    palette: {
      dim: PRISM_DIM.toArray(), accent: PRISM_ACCENT.toArray(), hot: PRISM_HOT.toArray(),
      placedLow: new THREE.Color(PLACED_LOW_DENSITY_COLOR).toArray(),
      placedHigh: new THREE.Color(PLACED_HIGH_DENSITY_COLOR).toArray(),
    },
  };
  var localBlockScene = createBlockScene(blockSceneConfig);

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
    // MAP.240: the part's vertices are offsets from its origin, which the mesh carries.
    if (part.origin) mesh.position.set(part.origin[0], part.origin[1], part.origin[2]);
    mesh.userData.relativeLayer = "block prisms";
    mesh.frustumCulled = false;
    mesh.visible = part.vertexCount > 0;
    mesh.userData.part = part;
    return mesh;
  }

  // --- Sector-block nodes --------------------------------------------------
  //
  // The map asks for the nodes of its cut rather than "everything within R of the
  // camera target": the blocks of the sector ladder (a block of 3 to 177,147 shells, or
  // a sector; planetgen.galaxy.drill) that the camera needs, big ones far and small
  // ones near, chosen by their size on screen (static/galaxynodecut.js, and
  // planetgen.galaxy.nodecut for the server's first frame). A node's key is
  // "m.shell.band.wedge". The cut opens blocks, biggest on screen first, until
  // NODE_BUDGET nodes are held; a block hands its place over to its children as it
  // grows on screen, so neither pops. Once the view's own nodes are in and the camera
  // has been still a moment, the nodes one click-zoom step in are fetched too
  // (prefetchNodes), so zooming in never waits on them; not on a metered or slow
  // connection.
  //
  // A node's contents depend only on its key and the database's content stamp, so
  // nodes are cached three ways: in memory here, in this browser's IndexedDB (so a
  // revisit or reload doesn't refetch them), and on the server's disk
  // (planetgen/web/lib/tilecache.py) -- only nodes in none of those reach the API and
  // database. Every response carries the current stamp; when it differs from ours,
  // the response also lists which nodes changed since (history), and only those are
  // dropped and refetched. A change the history can't account for (a deleted
  // sector, a new release, a stamp older than the history) drops every cached node.

  var NODE_EDGE_PC = data.nodeEdgePc;
  var MAX_NODES_PER_REQUEST = data.maxNodesPerRequest || 128;
  var NODE_CUT_CONFIG = { edgePc: NODE_EDGE_PC, bounds: data.nodeBoundsPc || null, maxNodes: data.nodeBudget || 1200 };

  // Where the camera would be at `radius` from the target along its present
  // line of sight. The stage view turns the camera freely (`setCamera`), so
  // the line comes from the camera itself, not from the orbit's angles.
  function cameraAtRadius(radius) {
    var back = camera.position.clone().sub(target);
    var length = back.length();
    var offset = length > 0 ? back.multiplyScalar(radius / length) : offsetFromOrbit({ radius: radius, theta: orbit.theta, phi: orbit.phi });
    return [target.x + offset.x, target.y + offset.y, target.z + offset.z];
  }

  // The camera's focal length in pixels, which turns a block's size into its size on screen.
  function focalPixels() {
    var height = canvasEl && canvasEl.clientHeight ? canvasEl.clientHeight : 800;
    return height / (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2));
  }

  // The cut of the camera at orbit radius `radius` (default: where it is): `keys` are its nodes, nearest first, and
  // `expanded` the blocks it opened (galaxynodecut.nodeCut). Fresh unless the camera has moved by only a few hundredths
  // of its distance from the target since the last cut at its own radius (and `fresh` is not asked for): the cut is the
  // costly part of a frame.
  var cutMemo = null;
  function neededNodes(radius, fresh) {
    var own = radius === undefined;
    var orbitRadius = Math.max(MIN_RADIUS, Math.min(MAX_RADIUS, own ? orbit.radius : radius));
    var cam = own ? camera.position.toArray() : cameraAtRadius(orbitRadius);
    var look = [target.x, target.y, target.z];
    var focal = focalPixels();
    var span = 0.03 * orbitRadius;
    if (own && !fresh && cutMemo && cutMemo.focal === focal
        && [0, 1, 2].every(function (a) { return Math.abs(cam[a] - cutMemo.cam[a]) <= span && Math.abs(look[a] - cutMemo.look[a]) <= span; })) {
      return cutMemo.need;
    }
    var cut = nodeCut(cam, look, Object.assign({ focalPx: focal }, NODE_CUT_CONFIG));
    if (own) {
      cutMemo = { cam: cam, look: look, focal: focal, need: cut };
    }
    return cut;
  }

  // --- Node caches ---------------------------------------------------------

  var NODE_MEMORY_MAX = 4000;
  // Keyed by database name.
  var STORAGE_PREFIX = "planetgen:node:" + data.storageKey + ":";
  // Where the stamp our stored nodes are at is kept, and the generation
  // (the label our stored nodes are filed under since they last all went
  // stale).
  var STAMP_RECORD = STORAGE_PREFIX + "@stamp";
  var remembered = readStampRecord();
  var currentStamp = remembered.stamp || "";
  var currentGeneration = remembered.generation || "";
  var nodeMemory = new Map();

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

  // The small record of where our stored nodes are at lives in localStorage,
  // which can be missing or throw on any access (private browsing, blocked
  // site data): every use is wrapped, and the map works the same without
  // it, just refetching from the server cache.
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

  // The nodes themselves live in IndexedDB (MAP.158), which holds far more
  // than localStorage's few megabytes (one sector link fills that). Every
  // use fails open to "nothing stored". Reads are asynchronous, so the
  // memory cache stays the one `getNode` reads; a node missing from it is
  // looked up here (loadStoredNodes) before the server is asked.
  var NODE_DB_NAME = "planetgen-nodes";
  var NODE_DB_STORE = "nodes";
  // The most packed node bytes kept, oldest use first out.
  var NODE_DB_MAX_BYTES = 64 * 1024 * 1024;
  var NODE_DB_TRIM_MS = 30000;
  var nodeDbPromise = null;
  var nodeDbTrimmedAt = 0;

  function openNodeDb() {
    if (!nodeDbPromise) {
      nodeDbPromise = new Promise(function (resolve) {
        try {
          if (typeof indexedDB === "undefined" || !indexedDB) {
            resolve(null);
            return;
          }
          var request = indexedDB.open(NODE_DB_NAME, 1);
          request.onupgradeneeded = function () {
            var db = request.result;
            var nodes = db.createObjectStore(NODE_DB_STORE, { keyPath: "id" });
            // gen: this database's generation ("db|generation"); at: when
            // the node was last stored or read.
            nodes.createIndex("gen", "gen");
            nodes.createIndex("scope", "scope");
            nodes.createIndex("at", "at");
          };
          request.onsuccess = function () { resolve(request.result); };
          request.onerror = function () { resolve(null); };
          request.onblocked = function () { resolve(null); };
        } catch (err) {
          resolve(null);
        }
      });
    }
    return nodeDbPromise;
  }

  function generationName(generation) {
    return data.storageKey + "|" + generation;
  }

  function nodeRecordId(generation, key) {
    return generationName(generation) + "|" + key;
  }

  // Runs `work(store, done)` in one transaction on the node store and
  // resolves with whatever `done` was given (undefined on any failure).
  function nodeDbRun(mode, work) {
    return openNodeDb().then(function (db) {
      if (!db) {
        return undefined;
      }
      return new Promise(function (resolve) {
        var result;
        try {
          var tx = db.transaction(NODE_DB_STORE, mode);
          tx.oncomplete = function () { resolve(result); };
          tx.onerror = function () { resolve(undefined); };
          tx.onabort = function () { resolve(undefined); };
          work(tx.objectStore(NODE_DB_STORE), function (value) { result = value; });
        } catch (err) {
          resolve(undefined);
        }
      });
    });
  }

  // The stored nodes among `keys` for `generation`, as a Map key -> packed
  // bytes (tilepack.js reads them).
  function readStoredNodes(generation, keys) {
    if (!generation || !keys.length) {
      return Promise.resolve(new Map());
    }
    return nodeDbRun("readwrite", function (store, done) {
      var found = new Map();
      var now = Date.now();
      keys.forEach(function (key) {
        var request = store.get(nodeRecordId(generation, key));
        request.onsuccess = function () {
          var record = request.result;
          if (record && record.bytes) {
            found.set(key, record.bytes);
            record.at = now;
            store.put(record);
          }
        };
      });
      done(found);
    }).then(function (found) { return found || new Map(); });
  }

  function writeStoredNodes(generation, nodes) {
    if (!generation || !nodes.length) {
      return;
    }
    var now = Date.now();
    nodeDbRun("readwrite", function (store) {
      nodes.forEach(function (entry) {
        store.put({
          id: nodeRecordId(generation, entry.key), scope: data.storageKey, gen: generationName(generation),
          key: entry.key, at: now, size: entry.bytes.length, bytes: entry.bytes,
        });
      });
    }).then(trimStoredNodes);
  }

  function removeStoredNodes(generation, keys) {
    if (!generation || !keys.length) {
      return;
    }
    nodeDbRun("readwrite", function (store) {
      keys.forEach(function (key) { store.delete(nodeRecordId(generation, key)); });
    });
  }

  // Drops this database's stored nodes of every other generation (and of
  // this one too when `keepCurrent` is false).
  function purgeStoredNodes(keepCurrent) {
    nodeDbRun("readwrite", function (store) {
      var request = store.index("scope").openCursor(IDBKeyRange.only(data.storageKey));
      request.onsuccess = function () {
        var cursor = request.result;
        if (!cursor) {
          return;
        }
        if (!keepCurrent || cursor.value.gen !== generationName(currentGeneration)) {
          cursor.delete();
        }
        cursor.continue();
      };
    });
  }

  // Keeps the store under NODE_DB_MAX_BYTES by dropping the nodes used
  // longest ago, at most every NODE_DB_TRIM_MS.
  function trimStoredNodes() {
    var now = Date.now();
    if (now - nodeDbTrimmedAt < NODE_DB_TRIM_MS) {
      return;
    }
    nodeDbTrimmedAt = now;
    nodeDbRun("readwrite", function (store) {
      var kept = 0;
      var request = store.index("at").openCursor(null, "prev");
      request.onsuccess = function () {
        var cursor = request.result;
        if (!cursor) {
          return;
        }
        kept += cursor.value.size || 0;
        if (kept > NODE_DB_MAX_BYTES) {
          cursor.delete();
        }
        cursor.continue();
      };
    });
  }

  // What the browser held before nodes replaced cube tiles: the old stamp record and any tiles localStorage held, and
  // the old IndexedDB database.
  function dropOldTileEntries() {
    try {
      if (typeof indexedDB !== "undefined" && indexedDB) {
        indexedDB.deleteDatabase("planetgen-tiles");
      }
    } catch (err) {
      // Nothing more to do.
    }
    var store = storage();
    if (!store) {
      return;
    }
    try {
      var doomed = [];
      for (var i = 0; i < store.length; i++) {
        var name = store.key(i);
        if (name && name.indexOf("planetgen:tile:") === 0 && name !== STAMP_RECORD) {
          doomed.push(name);
        }
      }
      doomed.forEach(function (name) { store.removeItem(name); });
    } catch (err) {
      // Nothing more to do.
    }
  }

  // Looks the nodes that are not in memory up in IndexedDB and puts what it
  // has into memory. Resolves with how many it found; a result for a
  // generation that has since changed is dropped.
  function loadStoredNodes(keys) {
    var generation = currentGeneration;
    return readStoredNodes(generation, keys).then(function (found) {
      if (generation !== currentGeneration) {
        return 0;
      }
      var loaded = 0;
      found.forEach(function (bytes, key) {
        if (memoryGet(nodeMemory, key) !== undefined) {
          return;
        }
        try {
          memorySet(nodeMemory, key, unpackNode(bytes), NODE_MEMORY_MAX);
          loaded += 1;
        } catch (err) {
          // A damaged record is fetched again, and replaced.
          removeStoredNodes(generation, [key]);
        }
      });
      return loaded;
    });
  }

  function getNode(key) {
    return memoryGet(nodeMemory, key);
  }

  // Remembers nodes from a response in memory and, once, in IndexedDB (as
  // the packed bytes the server sent: `packed` maps a key to them; the
  // page's own opening nodes have none and is only kept in memory).
  function putNodes(nodes, packed) {
    var entries = [];
    Object.keys(nodes).forEach(function (key) {
      var node = nodes[key];
      memorySet(nodeMemory, key, node, NODE_MEMORY_MAX);
      // A node the database was too busy to finish (PERF.76) is drawn as it is
      // but not kept for the next visit; it is asked for again below.
      if (!node.incomplete && packed && packed.has(key)) {
        entries.push({ key: key, bytes: packed.get(key) });
      }
    });
    writeStoredNodes(currentGeneration, entries);
  }

  // The node keys that changed between fromStamp and toStamp, from a
  // response's history (checks oldest first, each {from, to, nodes}), or
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
      if (!entry || entry.from !== at || !Array.isArray(entry.nodes)) {
        return null;
      }
      keys = keys.concat(entry.nodes);
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
    if (hadStamp) {
      loadRegions();
    }
    if (stale) {
      stale.forEach(function (key) {
        nodeMemory.delete(key);
      });
      removeStoredNodes(currentGeneration, stale);
    } else {
      currentGeneration = payload.generation || stamp;
      nodeMemory.clear();
      purgeStoredNodes(true);
    }
    rememberStamp();
    return true;
  }

  function absorb(payload, packed) {
    if (!payload) {
      return false;
    }
    var stampChanged = adoptStamp(payload);
    putNodes(payload.nodes || {}, packed);
    return stampChanged;
  }

  absorb(data.initial);
  purgeStoredNodes(true);
  dropOldTileEntries();

  // --- Clouds ----------------------------------------------------------------
  //
  // Each node lists the nebulae and supernova remnants that reach into it
  // (queryDb.galaxy_clouds_in_box); a cloud is a sprite the size of its
  // sphere, drawn over the blocks without a depth test (a flat sprite
  // through a cloud's middle would otherwise be cut in half by the blocks
  // around it, though the cloud fills that space). One that is only a pixel or two across
  // is hidden, and one the camera comes close to fades out so it never
  // fills the screen.
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
  // A nebula's shape mesh (MAP.103): fetched the first time the nebula is
  // big enough on screen, then drawn in place of its sprite.
  var cloudMeshes = new Map();
  var meshRequested = new Set();
  var loadNebulaShape = sceneData.nebulaShapePath ? createShapeLoader(sceneData.nebulaShapePath) : null;

  function requestNebulaMesh(sprite) {
    var cloud = sprite.userData.cloud;
    var key = cloudKey(cloud);
    if (!loadNebulaShape || cloud.type !== "nebula" || meshRequested.has(key)) {
      return;
    }
    meshRequested.add(key);
    loadNebulaShape(cloud.id, "low").then(function (shape) {
      if (cloudSprites.get(key) !== sprite) {
        return;
      }
      var look = NEBULA_LOOKS[cloud.descriptor] || DEFAULT_NEBULA_LOOK;
      var mesh = buildNebulaMesh(THREE, shape, cloud.radius_pc, look);
      mesh.position.copy(sprite.position);
      mesh.renderOrder = 4;
      mesh.visible = false;
      cloudGroup.add(mesh);
      cloudMeshes.set(key, mesh);
      updateClouds();
    }, function () {
      meshRequested.delete(key);
    });
  }

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
        var gone = cloudMeshes.get(key);
        if (gone) {
          cloudGroup.remove(gone);
          disposeNebulaMesh(gone);
        }
        cloudMeshes.delete(key);
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
    var drawn = 0;
    cloudSprites.forEach(function (sprite, key) {
      var cloud = sprite.userData.cloud;
      var distance = camera.position.distanceTo(sprite.position);
      var clipped = !drawnByTheGalaxy(sprite.position.x, sprite.position.y, sprite.position.z);
      var pixels = cloudPixels(cloud, distance);
      var mesh = cloudMeshes.get(key);
      if (cloud.type === "nebula" && !clipped && pixels >= NEBULA_MESH_MIN_PX) {
        requestNebulaMesh(sprite);
      }
      if (mesh && !clipped && pixels >= NEBULA_MESH_MIN_PX) {
        // The mesh takes over from the sprite, which stays (unseen) for aiming.
        mesh.visible = true;
        sprite.visible = true;
        sprite.material.opacity = 0;
        sprite.userData.share = 1;
        drawn++;
        return;
      }
      if (mesh) {
        mesh.visible = false;
      }
      var near = THREE.MathUtils.smoothstep(distance / cloud.radius_pc, CLOUD_NEAR_FADE[0], CLOUD_NEAR_FADE[1]);
      // MAP.155: opacity ramps in between 1 and 4 px of radius (a hard
      // edge at 2 px before); `share` is what can be aimed at.
      var share = near * cloudSizeShare(pixels);
      sprite.visible = share > 0 && !clipped;
      sprite.material.opacity = CLOUD_OPACITY * share;
      sprite.userData.share = clipped ? 0 : share;
    });
    // How many nebulae are drawn from their shape now (read by the browser tests).
    canvasEl.dataset.nebulaMeshes = String(drawn);
  }

  // The smallest visible cloud under a screen point that is still small
  // enough to aim at, as {cloud, core} (core: the point is within the
  // inner CLOUD_CORE of its radius), or null.
  var CLOUD_CORE = 0.5;

  function cloudAtRay(ray) {
    var best = null;
    var bestOffset = 0;
    var bestDistance = 0;
    cloudSprites.forEach(function (sprite) {
      var cloud = sprite.userData.cloud;
      if (!sprite.visible || !visibleEnough(sprite.userData.share) || (best && best.radius_pc <= cloud.radius_pc)) {
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
        bestDistance = distance;
      }
    });
    return best ? { cloud: best, core: bestOffset <= CLOUD_CORE * best.radius_pc, distance: bestDistance } : null;
  }

  // --- Stars ---------------------------------------------------------------
  //
  // Each node lists its most luminous pre-placed stars (bright_stars, v43:
  // every star of 500 L☉ or more, placed before any sector is filled), so
  // the arms show before anything is generated, and (MAP.51) the most
  // luminous stars of its generated systems down to a floor that drops
  // fourfold with each finer node level (queryDb.generated_star_floor_sol):
  // closing in on filled sectors brings out fainter and fainter stars,
  // down to the red dwarfs at sector depth. Boss: "no matter how far the
  // user zooms in they should be very small with a big glow", and "rough
  // sizes relative to the size of the star, brightness relative to
  // luminosity, and color relative to temperature". MAP.208: far from the
  // camera each is a point of light POINT_PX across whatever its size
  // (static/starpoint.js), only its brightness telling stars apart; from
  // about a sector away (NEAR_SECTOR_PC, the same constant the sector scene
  // uses) it grows smoothly into a core sized by the star's radius
  // (a little over one pixel to two and a half) in a soft halo of nearly the same size
  // for every star, its strength growing with luminosity (MAP.241: "glow brighter but not be bigger"), the core dimmer
  // for a fainter star, all in the star's blackbody color; the faint end is drawn brighter (MAP.87). One look for
  // both maps: static/starlight.js's `starLight`. Halos blend normally rather than adding up, so a crowded arm zoomed
  // out glows in its stars' colors instead of burning to white; stars are
  // drawn over the blocks, never hidden by them.
  var STAR_LOG_LUMINOSITY = [-4, 6];
  // A click within this many pixels of a black hole, neutron star or
  // quasar's center picks it (stars themselves are never picked, MAP.101).
  var POINT_PICK_PX = 7;
  // MAP.243: how near the pointer must be to a star beyond the outline to pick it (a star is a faint point).
  var HALO_PICK_PX = 10;
  // Stars new to the view fade in over this long (MAP.48: nothing pops
  // in); with reduced motion they just appear.
  var STAR_FADE_IN_MS = reducedMotion ? 0 : 300;

  // MAP.80: the color a black hole, neutron star or quasar is drawn in
  // among the stars -- off the blackbody line no star is drawn off (purple,
  // deep blue, pink). MAP.208: everything else about it is a star's: it is a
  // point of light and its brightness follows its mass, as bright as the
  // star of that mass (the server sends that luminosity, a black hole a
  // little dimmer), never brighter than the brightest star drawn.
  var POINT_COLORS = {
    black_hole: "#a43cff",
    neutron_star: "#3200ff",
    quasar: "#ff7fd0",
  };
  var DEFAULT_POINT_COLOR = "#ffffff";

  // MAP.111: with "Charted only" on, a star outside charted (generated)
  // sectors is drawn at this share of its opacity, as the blocks are
  // (galaxyblocks.buildCells). A generated star is charted, a bright
  // star once its sector is filled (`charted`); black holes, neutron
  // stars and quasars are landmarks and stay lit.

  var starMaterial = new THREE.ShaderMaterial({
    uniforms: {
      pixelRatio: { value: renderer.getPixelRatio() }, now: { value: 0 },
      fadeIn: { value: Math.max(STAR_FADE_IN_MS, 1) / 1000 },
      chartedOnly: { value: 0 }, unchartedFarKeep: { value: UNCHARTED_FAR_KEEP }, focusPc: { value: 0 },
      zoomFadeOn: { value: 1 },
      // MAP.327: the share of its light a star outside the opened sector keeps (1 with none opened).
      outsideDim: { value: 1 },
      magLimit: { value: 1000 }, magOn: { value: 1 },
      camHigh: { value: new THREE.Vector3() }, camLow: { value: new THREE.Vector3() },
      viewRotation: { value: new THREE.Matrix3() },
    },
    vertexShader: [
      "#include <common>",
      "#include <logdepthbuf_pars_vertex>",
      "attribute float starBorn;",
      // MAP.267: positions are relative-to-eye pairs: `position` the high float32, `starLow` the rest, and the camera
      // has its own pair; the sum of the two differences is the star's place in the camera's frame.
      "attribute vec3 starLow;",
      "uniform vec3 camHigh;",
      "uniform vec3 camLow;",
      "uniform mat3 viewRotation;",
      "uniform float now;",
      "uniform float fadeIn;",
      // MAP.152: the share of the node that lists the star (static/starfade.js),
      // worked out on the page as the camera moves (updateStarShares).
      "attribute float starNode;",
      "uniform float zoomFadeOn;",
      // MAP.148: the apparent-magnitude law (static/starmagnitude.js's twin).
      "attribute float starAbsMag;",
      // MAP.251: how many magnitudes dimmer a halo or satellite-galaxy star is judged (starmagnitude.js).
      "attribute float starTierMag;",
      "uniform float magLimit;",
      "uniform float magOn;",
      STAR_MAGNITUDE_GLSL,
      STAR_POINT_GLSL,
      "varying float vShown;",
      "attribute float starSize;",
      "attribute float starCore;",
      // x: the grown halo's strength; y, z: the far halo's strength and width beyond the point (one attribute: a
      // vertex shader has 16 at most).
      "attribute vec3 starGlow;",
      "attribute float starBright;",
      "attribute float starPointBright;",
      "attribute vec3 starColor;",
      "uniform float pixelRatio;",
      "attribute float starClipped;",
      // MAP.252: the share of its light an object keeps under "Charted only" (1 for a charted one).
      "attribute float starKeep;",
      "uniform float chartedOnly;",
      "uniform float outsideDim;",
      "uniform float unchartedFarKeep;",
      "uniform float focusPc;",
      UNCHARTED_FAR_GLSL,
      "varying vec3 vColor;",
      "varying float vCore;",
      "varying float vGlow;",
      "varying float vBright;",
      "varying float vWhiten;",
      "void main() {",
      "  vColor = starColor;",
      "  float zoomShown = starNode;",
      "  vShown = clamp((now - starBorn) / fadeIn, 0.0, 1.0) * mix(1.0, zoomShown, zoomFadeOn);",
      "  vec3 relative = (position - camHigh) + (starLow - camLow);",
      "  vec4 mvPosition = vec4(viewRotation * relative, 1.0);",
      // MAP.187: stars are not part of the near field -- it thins the blocks'
      // density and cell boxes only, so a star never vanishes in front of
      // the camera or on the line to the focus.
      // MAP.148: opacity and brightness follow the apparent magnitude from the camera.
      "  float dist = length(mvPosition.xyz);",
      // "Charted only" (MAP.111, MAP.252): an object outside charted space dimmed by its place in the cut order, and
      // the farther from the focus, the more.
      "  float keep = starKeep < 0.999 ? starKeep * mix(1.0, unchartedFarKeep, unchartedFar(dist, focusPc)) : 1.0;",
      "  vShown *= mix(1.0, keep, chartedOnly);",
      "  vShown *= outsideDim;",
      "  float appMag = apparentMagnitude(starAbsMag, dist) + starTierMag;",
      "  vShown *= mix(1.0, magnitudeShare(appMag, magLimit), magOn);",
      "  gl_Position = projectionMatrix * mvPosition;",
      // MAP.208: a point of light POINT_PX across until the camera is about
      // a sector away, then grown to its size and glow (static/starpoint.js's
      // pointSizePx and pointCoreShare); MAP.155: a star arriving with the
      // zoom holds back as a point until its zoom share rises.
      "  float arrive = mix(1.0, zoomShown, zoomFadeOn);",
      "  float grow = pointGrow(dist) * arrive;",
      "  float pointPx = " + POINT_PX.toFixed(2) + ";",
      // A far point keeps a halo of its own luminosity in its own color (starpoint.js's pointSizePx).
      "  float farPx = pointPx + starGlow.z;",
      "  float sizePx = farPx + max(starSize - farPx, 0.0) * grow;",
      "  float corePx = pointPx + max(starCore - pointPx, 0.0) * grow;",
      "  gl_PointSize = sizePx * pixelRatio;",
      "  vCore = corePx / sizePx;",
      "  vGlow = mix(starGlow.y, starGlow.x, grow);",
      "  vWhiten = mix(" + POINT_WHITEN.toFixed(3) + ", " + GROWN_WHITEN.toFixed(3) + ", grow);",
      "  vBright = mix(starPointBright, starBright, grow) * mix(1.0, fluxGain(appMag, magLimit), magOn);",
      // Inside the cell of the sector opened in place: dropped (it draws its own).
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
      "varying float vWhiten;",
      "void main() {",
      "  #include <logdepthbuf_fragment>",
      "  float r = length(gl_PointCoord * 2.0 - 1.0);",
      "  if (r > 1.0) discard;",
      "  float core = 1.0 - smoothstep(vCore * 0.5, vCore, r);",
      "  float halo = vGlow * exp(-r * r * 4.0) * (1.0 - r);",
      "  gl_FragColor = vec4(mix(vColor, vec3(1.0), core * vWhiten * vBright), clamp(core * vBright + halo, 0.0, 1.0) * vShown);",
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
  // MAP.266: the unresolved glow, under the stars (the Sector Map has none: its sector is all resolved).
  var glow = data.pinned ? null : createGlowLayer(scene, camera, canvasEl, { light: isLightBackground() });
  // MAP.316: the charted regions, soft spheres at their extents coloured by share charted, drawn by a layer of their own
  // (under the glow's) that fades each as its stars resolve; the list in the Menu and the names come from the same blobs.
  var regionLayer = data.pinned || !data.regionsPath ? null : createGlowLayer(scene, camera, canvasEl, { light: isLightBackground(), look: regionLook, renderOrder: 3 });
  // MAP.176: the groups of nebulae each drawn node lists (`nebula_regions`), soft sprites in the dominant family's colour with a thin
  // ring, a layer of their own under the charted regions' that holds itself to REGION_BUDGET and fades each as its node hands it down.
  var nebulaLayer = data.pinned ? null : createGlowLayer(scene, camera, canvasEl, {
    light: isLightBackground(), look: nebulaRegionLook, cut: cutNebulaRegions, renderOrder: 3,
  });
  var nebulaEntries = [];
  var nebulaVersion = 0;
  var nebulaHeld = false;
  var nebulaFamilyOn = {};
  NEBULA_FAMILIES.forEach(function (family) { nebulaFamilyOn[family] = true; });
  function setNebulaRegions() {
    if (!nebulaLayer) return;
    nebulaLayer.setBlobs(nebulaRegionBlobs(nebulaEntries, NEBULA_LOOKS, isLightBackground(), function (family) {
      return nebulaFamilyOn[family] !== false;
    }), [0, 0, 0]);
    nebulaVersion += 1;
    // The legend's line about rings shows only while some group is drawn.
    var held = nebulaLayer.count() > 0;
    if (held !== nebulaHeld) {
      nebulaHeld = held;
      if (stageView && stageView.layers) stageView.layers.refresh();
    }
  }
  var regionBlobList = [];
  var regionsVersion = 0;
  var regionsOn = true;
  var regionList = null;
  var regionByKey = new Map();
  // A click lands on a region's sphere within this share of its radius on the screen (at least this many pixels), and only
  // while it is drawn at least this strongly: the soft rim is not the region.
  var REGION_PICK_SHARE = 0.7;
  var REGION_PICK_MIN_PX = 9;
  var REGION_PICK_MIN_ALPHA = 0.05;
  var regionSpot = new THREE.Vector3();
  function setRegions(list) {
    regionBlobList = list;
    regionByKey = new Map(list.map(function (blob) { return [blob.key, blob]; }));
    regionsVersion += 1;
    markersDirty = true;    // the inhabited tokens of the regions (MAP.321 c)
    regionLayer.setBlobs(list, [0, 0, 0]);
    if (regionList) regionList.set(list);
    refreshRegionPois();
  }
  // MAP.357: from a galactic view every region with filled-in sectors stands as a point of interest (the region symbol) with a link that goes
  // to it, the camera framing the whole region. Regions closer together on the screen than REGION_POI_MERGE_PX are one point (the biggest
  // region's place, a count beside it) whose link frames all of them, so a point stands for what the zoom can tell apart. A region the camera
  // is in (within REGION_POI_REACH of its radii of its middle) gets none: its sectors and stars are the points there.
  var REGION_POI_REACH = 2.5;
  var REGION_POI_MERGE_PX = 48;
  var regionPoiLayer = null;
  var regionPoiEls = new Map();
  var regionPoiGroups = [];
  var lastPoiSignature = "";
  var poiSpot = new THREE.Vector3();
  function regionPoiGroupsNow() {
    var width = canvasEl.clientWidth || 1;
    var height = canvasEl.clientHeight || 1;
    var items = [];
    regionBlobList.forEach(function (blob) {
      if (!(blob.sectorCount > 0)) return;
      var away = Math.hypot(camera.position.x - blob.x, camera.position.y - blob.y, camera.position.z - blob.z) > REGION_POI_REACH * blob.radiusPc;
      if (!away) return;
      poiSpot.set(blob.x, blob.y, blob.z).project(camera);
      if (poiSpot.z < -1 || poiSpot.z > 1 || Math.abs(poiSpot.x) > 1 || Math.abs(poiSpot.y) > 1) return;
      items.push({ blob: blob, x: (poiSpot.x + 1) / 2 * width, y: (1 - poiSpot.y) / 2 * height });
    });
    items.sort(function (a, b) { return (b.blob.sectorCount - a.blob.sectorCount) || (a.blob.key < b.blob.key ? -1 : 1); });
    var groups = [];
    items.forEach(function (item) {
      var home = groups.find(function (group) { return Math.hypot(group.x - item.x, group.y - item.y) <= REGION_POI_MERGE_PX; });
      if (home) home.members.push(item.blob);
      else groups.push({ x: item.x, y: item.y, head: item.blob, members: [item.blob] });
    });
    return groups;
  }
  // The sphere that holds every region of a group: its middle is the members' mean, its radius reaches the farthest member's far side.
  function regionGroupExtent(group) {
    var n = group.members.length;
    var centre = [0, 0, 0];
    group.members.forEach(function (blob) { centre[0] += blob.x / n; centre[1] += blob.y / n; centre[2] += blob.z / n; });
    var radius = 0;
    group.members.forEach(function (blob) { radius = Math.max(radius, Math.hypot(blob.x - centre[0], blob.y - centre[1], blob.z - centre[2]) + blob.radiusPc); });
    return { centre: centre, radius: radius };
  }
  function openRegionGroup(group) {
    if (group.members.length === 1) {
      openRegionBlob(group.head);
      return;
    }
    if (!regionOpener) return;
    regionOpener.region({ block: group.head.block, name: group.head.name, extent: regionGroupExtent(group) });
  }
  function refreshRegionPois() {
    if (!viewport || data.pinned) return;
    var show = regionsOn && !hiddenCell && regionBlobList.length > 0;
    regionPoiGroups = show ? regionPoiGroupsNow() : [];
    if (!regionPoiLayer) {
      if (!regionPoiGroups.length) return;
      regionPoiLayer = document.createElement("div");
      regionPoiLayer.className = "map-region-pois";
      viewport.appendChild(regionPoiLayer);
    }
    var seen = new Set();
    regionPoiGroups.forEach(function (group) {
      var key = group.head.key;
      seen.add(key);
      var el = regionPoiEls.get(key);
      if (!el) {
        el = document.createElement("button");
        el.type = "button";
        el.className = "map-region-poi";
        el.dataset.region = key;
        var symbol = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        symbol.setAttribute("class", "poi-icon poi-region");
        symbol.setAttribute("aria-hidden", "true");
        symbol.setAttribute("focusable", "false");
        symbol.dataset.poi = "region";
        var use = document.createElementNS("http://www.w3.org/2000/svg", "use");
        use.setAttribute("href", new URL("poi/region.svg" + VERSION_QUERY + "#poi", import.meta.url).href);
        symbol.appendChild(use);
        var count = document.createElement("span");
        count.className = "map-region-poi-count";
        el.append(symbol, count);
        el.addEventListener("click", function () {
          var mine = regionPoiGroups.find(function (g) { return g.head.key === key; });
          if (mine) openRegionGroup(mine);
        });
        regionPoiLayer.appendChild(el);
        regionPoiEls.set(key, el);
      }
      var many = group.members.length > 1;
      var name = group.head.name || "Charted region";
      var text = many
        ? group.members.length + " charted regions, the biggest " + name + " (" + group.head.sectorCount + " sectors): go to all of them"
        : name + " (" + group.head.sectorCount + " sectors): go to the region";
      el.title = text;
      el.setAttribute("aria-label", text);
      el.querySelector(".map-region-poi-count").textContent = many ? String(group.members.length) : "";
      el.hidden = false;
      el.style.left = group.x.toFixed(1) + "px";
      el.style.top = group.y.toFixed(1) + "px";
    });
    regionPoiEls.forEach(function (el, key) {
      if (seen.has(key)) return;
      el.remove();
      regionPoiEls.delete(key);
    });
  }
  canvasEl.galaxyRegionPois = function () {
    return regionPoiGroups.map(function (group) {
      return { key: group.head.key, name: group.head.name, count: group.members.length, x: group.x, y: group.y, keys: group.members.map(function (blob) { return blob.key; }) };
    });
  };
  function loadRegions() {
    if (!regionLayer) return;
    fetch(data.regionsPath).then(function (response) {
      if (!response.ok) throw new Error(data.regionsPath + " returned " + response.status);
      return response.json();
    }).then(function (payload) {
      setRegions(regionBlobs(payload.regions));
    }).catch(function () {
      // The regions are an aid to finding known space; the map is whole without them, and the next change of the galaxy asks again.
    });
  }
  // MAP.267: the stars are one buffer filled node by node. Each drawn node's stars are a chunk at the end of the
  // buffer (`starSlots`, starbuffer.js); a node that leaves is cut out by moving the tail into its place, so the
  // buffer stays dense from 0 to `starList.length` and a node landing is a write at the end, not a rebuild. Positions
  // are relative-to-eye pairs (starbuffer.js): `position` is the high float32 of each coordinate and `starLow` what
  // that left over, and the shader subtracts the camera's own pair (`camHigh`, `camLow`), so a star keeps under a
  // quarter pixel down to about 0.2 AU from the camera at any distance from the galaxy's origin.
  var starList = [];
  var starSlots = new SlotStore();
  var STAR_ATTRIBUTES = [
    ["position", 3], ["starLow", 3], ["starColor", 3], ["starSize", 1], ["starCore", 1], ["starGlow", 3],
    ["starBright", 1], ["starPointBright", 1], ["starAbsMag", 1], ["starTierMag", 1], ["starBorn", 1],
    ["starClipped", 1], ["starKeep", 1], ["starNode", 1],
  ];
  var STAR_CAPACITY_MIN = 4096;
  var starCapacity = 0;
  var starArrays = {};
  // When each drawn star first showed (seconds, performance.now's clock) and how many chunks list it, by starKey, so
  // a star already on screen never fades in again.
  var starBornByKey = new Map();
  // The brightest star drawn so far: a phenomenon is never drawn brighter (MAP.208). It only grows until the buffer empties.
  var brightestSolSeen = 0;
  // Slots to upload: `[start, count]` runs, or true for all of them.
  var starDirty = [];

  function starClock() {
    return performance.now() / 1000;
  }

  // A star's key among the drawn ones: bright_stars, stars and each
  // phenomenon table's rows are numbered separately.
  function starKey(star) {
    return star.phenomenon ? star.type + star.id : (star.generated ? "s" : "b") + star.id;
  }

  // Where `value`'s log10 falls in `range`, 0..1.
  function logShare(value, range) {
    return THREE.MathUtils.clamp((Math.log10(Math.max(value, 1e-12)) - range[0]) / (range[1] - range[0]), 0, 1);
  }

  // The node fade (static/starfade.js) of each drawn star, as the node's
  // block (MAP.152), in the order of starList.
  // Counts the star sets drawn, so a new set is noticed by the magnitude limit (MAP.148).
  var starSetVersion = 0;
  var starFadeList = [];

  // Makes room for `n` stars: the buffers double (a new geometry, copied across) when they are full.
  function ensureStarCapacity(n) {
    if (n <= starCapacity) {
      return;
    }
    var capacity = Math.max(STAR_CAPACITY_MIN, starCapacity);
    while (capacity < n) {
      capacity *= 2;
    }
    var used = starList.length;
    var geometry = new THREE.BufferGeometry();
    var arrays = {};
    STAR_ATTRIBUTES.forEach(function (pair) {
      var array = new Float32Array(capacity * pair[1]);
      if (starArrays[pair[0]]) {
        array.set(starArrays[pair[0]].subarray(0, used * pair[1]));
      }
      arrays[pair[0]] = array;
      geometry.setAttribute(pair[0], new THREE.BufferAttribute(array, pair[1]).setUsage(THREE.DynamicDrawUsage));
    });
    starArrays = arrays;
    starCapacity = capacity;
    starPoints.geometry.dispose();
    starPoints.geometry = geometry;
    starDirty = true;
  }

  function markStarsDirty(start, count) {
    if (starDirty !== true && count > 0) {
      starDirty.push([start, count]);
    }
  }

  // Sets the geometry's draw range and uploads what changed since the last call.
  function flushStars() {
    var geometry = starPoints.geometry;
    var n = starList.length;
    geometry.setDrawRange(0, n);
    var runs = starDirty;
    starDirty = [];
    if (runs !== true && runs.length > 48) {
      runs = true;
    }
    STAR_ATTRIBUTES.forEach(function (pair) {
      var attribute = geometry.getAttribute(pair[0]);
      if (!attribute) {
        return;
      }
      if (runs === true) {
        attribute.clearUpdateRanges();
        attribute.needsUpdate = true;
      } else if (runs.length) {
        attribute.clearUpdateRanges();
        runs.forEach(function (run) {
          attribute.addUpdateRange(run[0] * pair[1], run[1] * pair[1]);
        });
        attribute.needsUpdate = true;
      }
    });
  }

  // Writes one star (an entry `{ star, fade }`) into slot `i`.
  function writeStar(i, star, now) {
    var a = starArrays;
    // A phenomenon is drawn as the star of its mass: `map_luminosity_sol` is the
    // server's (planetgen.physics.stellar_evolution.phenomenon_map_luminosity_sol).
    var luminosity = star.phenomenon ? phenomenonLuminositySol(star.map_luminosity_sol, brightestSolSeen) : star.luminosity_sol;
    var t = logShare(luminosity, STAR_LOG_LUMINOSITY);
    a.starKeep[i] = unchartedKeep(star, star.phenomenon ? massRampShare(star.mass_solar) : t);
    a.starAbsMag[i] = absoluteMagnitude(luminosity);
    a.starTierMag[i] = tierMagnitudePenalty(star.tier);
    var x = splitPair(star.x);
    var y = splitPair(star.y);
    var z = splitPair(star.z);
    a.position[3 * i] = x[0];
    a.position[3 * i + 1] = y[0];
    a.position[3 * i + 2] = z[0];
    a.starLow[3 * i] = x[1];
    a.starLow[3 * i + 1] = y[1];
    a.starLow[3 * i + 2] = z[1];
    var color = star.phenomenon ? new THREE.Color(POINT_COLORS[star.type] || DEFAULT_POINT_COLOR).toArray()
      : starColor(star.temperature_k);
    a.starColor[3 * i] = color[0];
    a.starColor[3 * i + 1] = color[1];
    a.starColor[3 * i + 2] = color[2];
    // MAP.241: the star's own look, one for both maps (static/starlight.js): a size in a narrow band, the luminosity in the
    // glow, the radius in the core; MAP.87's boost of the faint end is in it.
    var light = starLight(luminosity, star.radius_sol);
    a.starCore[i] = light.corePx;
    a.starSize[i] = light.sizePx;
    a.starGlow[3 * i] = light.glow;
    a.starGlow[3 * i + 1] = pointHaloGlow(t);
    a.starGlow[3 * i + 2] = pointHaloPx(t);
    a.starBright[i] = light.bright;
    a.starPointBright[i] = pointBrightness(t);
    a.starClipped[i] = clippedByTheSector(star);
    a.starNode[i] = 0;
    var key = starKey(star);
    var born = starBornByKey.get(key);
    if (born) {
      born.count += 1;
    } else {
      born = { at: now, count: 1 };
      starBornByKey.set(key, born);
    }
    a.starBorn[i] = born.at;
  }

  // Adds a node's stars (`entries`: `{ star, fade }` for each) as a chunk of the buffer, tagged with the node's revision.
  function addStarChunk(key, tag, entries) {
    var n = entries.length;
    var start = starSlots.add(key, n, tag);
    ensureStarCapacity(start + n);
    var now = starClock();
    entries.forEach(function (entry) {
      if (!entry.star.phenomenon && entry.star.luminosity_sol > brightestSolSeen) {
        brightestSolSeen = entry.star.luminosity_sol;
      }
    });
    entries.forEach(function (entry, offset) {
      writeStar(start + offset, entry.star, now);
      starList[start + offset] = entry.star;
      starFadeList[start + offset] = entry.fade;
    });
    markStarsDirty(start, n);
    starSync.written += n;
  }

  // Takes a node's chunk out, moving the tail of the buffer into the holes.
  function removeStarChunk(key) {
    var slots = starSlots.slotsOf(key);
    slots.forEach(function (slot) {
      var id = starKey(starList[slot]);
      var born = starBornByKey.get(id);
      if (born && --born.count <= 0) {
        starBornByKey.delete(id);
      }
    });
    var moves = starSlots.remove(key);
    var destinations = [];
    moves.forEach(function (move) {
      STAR_ATTRIBUTES.forEach(function (pair) {
        var array = starArrays[pair[0]];
        array.copyWithin(move.to * pair[1], move.from * pair[1], (move.from + 1) * pair[1]);
      });
      starList[move.to] = starList[move.from];
      starFadeList[move.to] = starFadeList[move.from];
      destinations.push(move.to);
    });
    starSync.moved += moves.length;
    starList.length = starSlots.n;
    starFadeList.length = starSlots.n;
    destinations.sort(function (p, q) { return p - q; });
    runsOf(destinations).forEach(function (run) { markStarsDirty(run[0], run[1]); });
    if (!starSlots.n) {
      brightestSolSeen = 0;
    }
  }

  function clearStarChunks() {
    starSlots.keys().forEach(removeStarChunk);
  }

  // --- The limiting magnitude (MAP.148, static/starmagnitude.js) ----------------
  //
  // Each time the view changes (at most every MAG_RECOMPUTE_MS), the stars on
  // screen are counted by apparent magnitude, weighted by the zoom opacity
  // MAP.153 gives them, and the limit is set where the count reaches the
  // target; the shader eases toward it, so the picture never jumps.
  // The share of the node each drawn star is listed in (MAP.152), by star in
  // starList's order; redone when the camera moves or the stars change.
  var starShares = new Float32Array(0);
  var starShareSignature = "";
  function updateStarShares() {
    var n = starList.length;
    var signature = starSetVersion + "," + n + "," + camera.position.x + "," + camera.position.y + "," + camera.position.z;
    if (signature === starShareSignature || !n) {
      return;
    }
    starShareSignature = signature;
    if (starShares.length !== n) {
      starShares = new Float32Array(n);
    }
    var weigh = drawnSettled ? weightCalculator(camera.position.toArray(), NODE_EDGE_PC, focalPixels(), drawnSettled) : null;
    var byNode = new Map();
    for (var i = 0; i < n; i++) {
      var fade = starFadeList[i];
      var share = 1;
      if (fade && weigh) {
        share = byNode.get(fade);
        if (share === undefined) {
          share = starShare(fade, weigh);
          byNode.set(fade, share);
        }
      }
      starShares[i] = share;
    }
    var attribute = starPoints.geometry.getAttribute("starNode");
    if (attribute) {
      attribute.array.set(starShares);
      attribute.clearUpdateRanges();
      attribute.addUpdateRange(0, n);
      attribute.needsUpdate = true;
    }
  }

  // MAP.267: the camera as a relative-to-eye pair for the star shader, and the rotation of its view.
  function updateStarEye() {
    camera.updateMatrixWorld();
    var uniforms = starMaterial.uniforms;
    var p = camera.position;
    var x = splitPair(p.x);
    var y = splitPair(p.y);
    var z = splitPair(p.z);
    uniforms.camHigh.value.set(x[0], y[0], z[0]);
    uniforms.camLow.value.set(x[1], y[1], z[1]);
    uniforms.viewRotation.value.setFromMatrix4(camera.matrixWorldInverse);
  }

  var MAG_RECOMPUTE_MS = 150;
  var MAG_TAU_S = 0.3;
  var magTarget = window.matchMedia("(max-width: 700px), (pointer: coarse)").matches ? TARGET_STARS_PHONE : TARGET_STARS;
  var magGoal = 1000;
  var magSignature = "";
  var magComputedAt = 0;
  var magLastFrame = 0;
  var magOnScreen = 0;

  function updateMagnitudeLimit(now) {
    updateStarShares();
    var uniforms = starMaterial.uniforms;
    uniforms.focusPc.value = nearUniforms.nearD.value;
    var seconds = magLastFrame ? (now - magLastFrame) / 1000 : 0;
    magLastFrame = now;
    var n = starList.length;
    var e = camera.matrixWorldInverse.elements;
    var signature = [starSetVersion, n, e[0], e[1], e[2], e[4], e[5], e[6], e[8], e[9], e[10], e[12], e[13], e[14], canvasEl.clientHeight].join(",");
    if (n && signature !== magSignature && now - magComputedAt >= MAG_RECOMPUTE_MS) {
      magSignature = signature;
      magComputedAt = now;
      var tanY = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2) * 1.1;
      var tanX = tanY * camera.aspect;
      var clipped = starPoints.geometry.getAttribute("starClipped");
      var absolute = starPoints.geometry.getAttribute("starAbsMag");
      var tierMag = starPoints.geometry.getAttribute("starTierMag");
      var magnitudes = new Float64Array(n);
      var onScreen = new Float64Array(n);
      var zoomed = new Float64Array(n);
      var count = 0;
      for (var i = 0; i < n; i++) {
        var star = starList[i];
        var w = starShares[i];
        zoomed[i] = w;
        if (!(w > 0) || (clipped && clipped.array[i] > 0.5)) continue;
        var vx = e[0] * star.x + e[4] * star.y + e[8] * star.z + e[12];
        var vy = e[1] * star.x + e[5] * star.y + e[9] * star.z + e[13];
        var vz = e[2] * star.x + e[6] * star.y + e[10] * star.z + e[14];
        var depth = -vz;
        if (depth <= 0 || Math.abs(vx) > depth * tanX || Math.abs(vy) > depth * tanY) continue;
        magnitudes[i] = apparentMagnitude(absolute.array[i], Math.sqrt(vx * vx + vy * vy + vz * vz)) + tierMag.array[i];
        onScreen[i] = w;
        count += 1;
      }
      magOnScreen = count;
      var focus = nearUniforms.nearD.value;
      magGoal = limitingMagnitude(
        histogramLimit(magnitudes, onScreen, magTarget),
        calibrationLimit(dimmestShown(absolute.array.subarray(0, n), zoomed), focus));
      if (!Number.isFinite(magGoal)) magGoal = 1000;
    }
    var current = uniforms.magLimit.value;
    uniforms.magLimit.value = reducedMotion || current >= 999 ? magGoal : approachLimit(current, magGoal, seconds, MAG_TAU_S);
  }

  // For tests (MAP.148): the limiting magnitude in use and the one the view
  // asks for, how many stars were counted on screen, the target, and the
  // apparent magnitude and opacity of every drawn star that is on screen.
  canvasEl.galaxyStarMagnitude = function () {
    var absolute = starPoints.geometry.getAttribute("starAbsMag");
    var tierMag = starPoints.geometry.getAttribute("starTierMag");
    var limit = starMaterial.uniforms.magLimit.value;
    var e = camera.matrixWorldInverse.elements;
    var shares = [];
    starList.forEach(function (star, i) {
      var d = Math.sqrt(Math.pow(star.x - camera.position.x, 2) + Math.pow(star.y - camera.position.y, 2) + Math.pow(star.z - camera.position.z, 2));
      var m = apparentMagnitude(absolute.array[i], d) + tierMag.array[i];
      shares.push([m, Math.min(1, Math.max(0, (limit - m) / RAMP))]);
    });
    return { limit: limit, goal: magGoal, target: magTarget, onScreen: magOnScreen, count: starList.length, stars: shares, ramp: RAMP, focus: nearUniforms.nearD.value, e: e[0] };
  };

  // For tests (MAP.149): the near field the shaders were given (the
  // camera-to-focus distance, the focus's radius, the center of the block
  // the camera is inside or null), the camera and its target, and the
  // near-field share at a world point if one is given.
  canvasEl.galaxyNearField = function (x, y, z) {
    var inside = nearUniforms.insideCenter.value;
    return {
      distance: nearUniforms.nearD.value, radius: nearUniforms.nearR.value,
      inside: inside.x === NO_BLOCK_INSIDE ? null : inside.toArray(),
      camera: camera.position.toArray(), target: target.toArray(),
      share: x === undefined ? null : nearFieldAt(x, y, z),
    };
  };

  // For tests (MAP.187): the star material's vertex shader and uniform
  // names, so a test can show a star is not in the near field.
  canvasEl.galaxyStarShader = function () {
    return { vertex: starMaterial.vertexShader, uniforms: Object.keys(starMaterial.uniforms) };
  };

  // For tests (MAP.153, MAP.152): the sum of the drawn stars' node shares if
  // the camera orbited at `radius` pc (the same arithmetic as the shader's),
  // and how many of them show more than half, so a test bounds the change
  // in the picture over a small zoom step without reading pixels.
  canvasEl.galaxyStarOpacity = function (radius) {
    var weigh = weightCalculator(cameraAtRadius(radius), NODE_EDGE_PC, focalPixels(), drawnSettled || undefined);
    var fade = starFadeList;
    var shown = 0;
    fade.forEach(function (f) {
      if (starShare(f, weigh) > 0.5) shown++;
    });
    return { sum: opacitySum(fade, weigh), shown: shown, count: fade.length, camera: orbit.radius };
  };

  // For tests (MAP.252): the share of its light each drawn object keeps under "Charted only" (1 for a charted one)
  // with the pieces it was worked out from, and whether "Charted only" is on.
  canvasEl.galaxyUnchartedKeeps = function () {
    var keeps = starPoints.geometry.getAttribute("starKeep");
    return {
      on: starMaterial.uniforms.chartedOnly.value > 0,
      far: starMaterial.uniforms.unchartedFarKeep.value,
      entries: starList.map(function (star, i) {
        return { id: star.id, phenomenon: Boolean(star.phenomenon), type: star.type || null, tier: star.tier || 0,
          mass: star.mass_solar == null ? null : star.mass_solar, keep: keeps.array[i] };
      }),
    };
  };

  // For tests: switches the zoom fade off (every star at full opacity, as
  // before MAP.153) so a pixel test of picking can find stars at any zoom.
  canvasEl.galaxyStarZoomFade = function (on) {
    starMaterial.uniforms.zoomFadeOn.value = on ? 1 : 0;
    starMaterial.uniforms.magOn.value = on ? 1 : 0;
  };

  // For tests: the frame the stars are drawn from and the farthest star
  // from it, in the points' own (offset) coordinates.
  // For tests (MAP.267): how the star buffer was filled: syncs so far, stars written and moved in all, and the last
  // sync's own figures (chunks added and removed, stars kept untouched, written, moved, in all, and its milliseconds).
  canvasEl.galaxyStarFill = function () {
    return { syncs: starSync.syncs, written: starSync.written, moved: starSync.moved, last: starSync.last,
      capacity: starCapacity, count: starList.length };
  };

  // For tests (MAP.241): the look each drawn star was given (the halo's width and strength, the core's width and opacity, the
  // blackbody colour), with the star's own luminosity and temperature.
  canvasEl.galaxyStarLooks = function () {
    var a = starArrays;
    return starList.map(function (star, i) {
      return {
        id: star.id, luminosity: star.luminosity_sol, temperature: star.temperature_k, phenomenon: !!star.phenomenon,
        sizePx: a.starSize[i], corePx: a.starCore[i], glow: a.starGlow[3 * i], bright: a.starBright[i],
        color: [a.starColor[3 * i], a.starColor[3 * i + 1], a.starColor[3 * i + 2]],
      };
    });
  };

  // For tests (MAP.338): whether the open system is being drawn from the eye's own frame, the scale it was drawn at, the camera's exact offset from
  // its target (not the difference of two positions), and how many frames were drawn so.
  canvasEl.galaxyDeepFrame = function () {
    return { active: deepFrame.active, scale: deepFrame.scale, passes: deepFrame.passes, back: eyeBack.toArray(), dist: eyeBack.length() };
  };

  canvasEl.galaxyStarFrame = function () {
    // MAP.267: the stars are relative-to-eye pairs. `far` is the farthest star from the camera, measured as the shader
    // does (each difference in float32), and `error` the largest gap between that and the star's own double position.
    var high = starArrays.position;
    var low = starArrays.starLow;
    var cam = camera.position.toArray().map(splitPair);
    var far = 0;
    var error = 0;
    for (var i = 0; i < starList.length; i++) {
      var star = starList[i];
      var truth = [star.x, star.y, star.z];
      var length = 0;
      for (var a = 0; a < 3; a++) {
        var offset = relativeOffset(high[3 * i + a], low[3 * i + a], cam[a][0], cam[a][1]);
        length += offset * offset;
        error = Math.max(error, Math.abs(offset - (truth[a] - camera.position.toArray()[a])));
      }
      far = Math.max(far, Math.sqrt(length));
    }
    return { far: far, error: error, count: starList.length, camera: camera.position.toArray(), capacity: starCapacity };
  };

  // For tests (MAP.264): how many stars the galaxy has, and how many it leaves out of its own drawing (those in the
  // opened sector's cell, which that sector draws).
  canvasEl.galaxyHiddenStars = function () {
    var attribute = starPoints.geometry.getAttribute("starClipped");
    var hidden = 0;
    if (attribute) {
      for (var i = 0; i < starList.length; i++) {
        if (attribute.array[i] > 0.5) hidden++;
      }
    }
    return { total: starList.length, hidden: hidden };
  };

  // 1 for a star the shader drops: one inside the opened sector's cell (it draws its own) and, with a sector opened, a
  // phenomenon outside it (no comets, rogue planets or black holes of the space around, MAP.327).
  function clippedByTheSector(star) {
    if (!drawnByTheGalaxy(star.x, star.y, star.z)) return 1;
    return hiddenCell && star.phenomenon ? 1 : 0;
  }

  // Marks the stars inside the opened sector's cell, which the shader drops.
  function markHiddenStars() {
    var attribute = starPoints.geometry.getAttribute("starClipped");
    if (!attribute) {
      return;
    }
    starList.forEach(function (star, i) {
      attribute.array[i] = clippedByTheSector(star);
    });
    attribute.clearUpdateRanges();
    attribute.addUpdateRange(0, starList.length);
    attribute.needsUpdate = true;
  }

  // "Charted only" (MAP.111): dims the stars outside charted sectors and
  // the blocks holding none, and outlines the charted blocks.
  function setChartedOnly(on) {
    starMaterial.uniforms.chartedOnly.value = on ? 1 : 0;
    stageView.setChartedOnly(on);
  }

  // --- Drawing from nodes --------------------------------------------------

  var starSignature = "";

  // A node's stars, bright and generated (each generated one marked),
  // and its black holes, neutron stars and quasars (MAP.80, marked
  // `phenomenon`, drawn as stars in their own colors), or none for a
  // missing node.
  function nodeStars(node) {
    if (!node) {
      return [];
    }
    if (!node.allStars) {
      Object.defineProperty(node, "allStars", {
        value: (node.stars || []).concat((node.generated || []).map(function (star) {
          return Object.assign({ generated: true }, star);
        }), (node.points || []).map(function (point) {
          return Object.assign({ phenomenon: true }, point);
        })),
      });
    }
    return node.allStars;
  }

  // MAP.123: whether the star filters (classes hidden, dimmest star shown)
  // let `star` through. Phenomena are never filtered here: they have their
  // own kind toggles.
  var STAR_FILTER_CLASSES = ["O", "B", "A", "F", "G", "K", "M", "other"];
  var lastNeed = null;
  function filterStamp() {
    if (typeof stageView === "undefined" || !stageView) return "";
    return stageView.hiddenClasses().sort().join("") + "/" + stageView.minLuminosity();
  }
  function starShown(star) {
    if (star.phenomenon || typeof stageView === "undefined" || !stageView) return true;
    var letter = String(star.star_type || "").charAt(0).toUpperCase();
    if (stageView.starClassHidden(STAR_FILTER_CLASSES.indexOf(letter) >= 0 ? letter : "other")) return false;
    var floor = stageView.minLuminosity();
    return !(floor > 0) || (star.luminosity_sol || 0) >= floor;
  }
  // The dimmest and brightest luminosity of the stars at this zoom (the
  // drawn nodes' stars, whatever the filters), the luminosity slider's
  // ends. Null with none loaded.
  var nodeLumRange = null;
  function noteLuminosityRange(keys) {
    var lo = Infinity, hi = 0;
    keys.forEach(function (key) {
      var node = getNode(key);
      if (node === undefined) return;
      nodeStars(node).forEach(function (star) {
        if (star.phenomenon) return;
        var lum = star.luminosity_sol;
        if (lum > 0) { lo = Math.min(lo, lum); hi = Math.max(hi, lum); }
      });
    });
    // Whatever the nodes of this view carry: far out their per-node
    // star caps leave only the brightest stars, and the scale follows that.
    nodeLumRange = hi > 0 ? [lo, hi] : null;
    if (typeof updateLumSlider === "function") updateLumSlider();
  }
  // Draws the stars again after a filter changed.
  function refreshStarFilters() {
    starSignature = "";
    if (lastNeed) renderFromCache(lastNeed);
  }

  // MAP.267: brings the star buffer to the stars of the nodes `drawn`, node by node: a node no longer drawn (or read
  // again, a new object) is cut out, a node new to the view is written at the end, and the nodes held already are not
  // touched. A change of the star filters starts the buffer afresh (the filters decide which stars a chunk holds).
  var starFilterStamp = "";
  var starSync = { syncs: 0, written: 0, moved: 0, last: null };
  function syncStars(drawn) {
    var began = performance.now();
    var filters = filterStamp();
    var kept = 0;
    var written0 = starSync.written;
    var moved0 = starSync.moved;
    var removed = 0;
    var added = 0;
    if (filters !== starFilterStamp) {
      starFilterStamp = filters;
      clearStarChunks();
    }
    var wanted = new Set(drawn);
    starSlots.keys().forEach(function (key) {
      var node = getNode(key);
      if (!wanted.has(key) || node === undefined || starSlots.tagOf(key) !== nodeRevision(node)) {
        removeStarChunk(key);
        removed += 1;
      } else {
        kept += starSlots.slotsOf(key).length;
      }
    });
    drawn.forEach(function (key) {
      if (starSlots.has(key)) {
        return;
      }
      var node = getNode(key);
      var fade = { block: parseNodeKey(key) };
      var entries = [];
      nodeStars(node).forEach(function (star) {
        if (starShown(star)) {
          entries.push({ star: star, fade: fade });
        }
      });
      addStarChunk(key, nodeRevision(node), entries);
      added += 1;
    });
    starSetVersion += 1;
    flushStars();
    starSync.syncs += 1;
    starSync.last = {
      added: added, removed: removed, kept: kept, written: starSync.written - written0, moved: starSync.moved - moved0,
      count: starList.length, ms: performance.now() - began,
    };
  }

  // Draws what is in hand of the cut `need` (neededNodes) and returns the keys of its nodes still to fetch. A node is
  // drawn for the share the cut gives it (galaxynodecut.weightCalculator): a block hands over to its children only
  // once they are all in hand (settledExpanded), so there is never a hole where the children are on their way; the
  // children of a node that is being drawn are in the buffer too, at no share yet, so they are there as the hand-over
  // reaches them.
  var drawnSettled = null;
  function renderFromCache(need) {
    lastNeed = need;
    markersDirty = true;
    // MAP.249: a view whose nodes are the ones already drawn (a fetch that brought nothing new, a settled camera)
    // is known from the nodes alone, so nothing is built or sorted.
    var view = drawnSignature(need.keys, getNode);
    var settled = settledExpanded(need, function (key) { return getNode(key) !== undefined; });
    var weigh = weightCalculator(camera.position.toArray(), NODE_EDGE_PC, focalPixels(), settled);
    var held = need.keys.filter(function (key) { return getNode(key) !== undefined; });
    var drawn = held.filter(function (key) {
      var block = parseNodeKey(key);
      if (weigh(block) > 0) return true;
      var parent = drillParent(block);
      return parent !== null && weigh(parent) > 0;
    });
    var nodesDrawn = currentStamp + "|" + filterStamp() + "|" + drawn.map(function (key) {
      return key + "#" + nodeRevision(getNode(key));
    }).sort().join(";");
    drawnSettled = settled;
    if (nodesDrawn === starSignature) {
      return view.missing;
    }
    var clouds = new Map();
    drawn.forEach(function (key) {
      (getNode(key).clouds || []).forEach(function (cloud) {
        clouds.set(cloudKey(cloud), cloud);
      });
    });
    noteLuminosityRange(drawn);
    starSignature = nodesDrawn;
    syncStars(drawn);
    var drawnNodes = drawn.map(function (key) {
      return { key: key, block: parseNodeKey(key), node: getNode(key) };
    });
    if (glow) {
      glow.setBlobs(blobsOf(drawnNodes, NODE_EDGE_PC), [target.x, target.y, target.z]);
    }
    nebulaEntries = drawnNodes;
    setNebulaRegions();
    var cloudKeys = Array.from(clouds.keys()).sort().join(",");
    if (cloudKeys !== cloudSignature) {
      cloudSignature = cloudKeys;
      setClouds(Array.from(clouds.values()));
    }
    return view.missing;
  }

  // --- Live node fetching --------------------------------------------------

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
  // fetch shortly after it stops. (With nodes, an aborted request's work
  // isn't wasted either: the server caches whatever it computed.)
  var MAX_CLICKS_PER_SECOND = 4;
  var MIN_MS_BETWEEN_IMMEDIATE_FETCHES = 1000 / MAX_CLICKS_PER_SECOND;
  var lastImmediateFetchAt = 0;

  function scheduleFetch(immediate) {
    if (fetchTimer) {
      clearTimeout(fetchTimer);
    }
    // The camera moved: not still, so no prefetch until it has been again.
    cancelPrefetch();
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

  // Bumped by every doFetch: a lookup in IndexedDB that finishes after a
  // newer fetch began lets that one do the work.
  var fetchRound = 0;

  function doFetch() {
    fetchTimer = null;
    var missing = renderFromCache(neededNodes(undefined, true));
    if (!missing.length) {
      schedulePrefetch();
      return;
    }
    var round = ++fetchRound;
    // What this browser kept from an earlier visit first (MAP.158).
    loadStoredNodes(missing).then(function (found) {
      if (round !== fetchRound) {
        return;
      }
      if (found) {
        missing = renderFromCache(neededNodes());
      }
      if (!missing.length) {
        schedulePrefetch();
        return;
      }
      fetchFromServer(missing);
    });
  }

  function fetchFromServer(missing) {
    if (activeAbort) {
      activeAbort.abort();
    }
    var controller = typeof AbortController !== "undefined" ? new AbortController() : null;
    activeAbort = controller;
    var params = new URLSearchParams({
      keys: missing.slice(0, MAX_NODES_PER_REQUEST).join(","),
    });
    if (currentStamp) {
      params.set("stamp", currentStamp);
    }
    fetch(data.fetchPath + "?" + params.toString(), controller ? { signal: controller.signal } : undefined)
      .then(function (response) {
        if (!response.ok) {
          throw new Error(data.fetchPath + " returned " + response.status);
        }
        return response.arrayBuffer();
      })
      .then(function (body) {
        var unpacked = unpackResponse(body);
        var payload = Object.assign({ nodes: unpacked.nodes }, unpacked.header);
        if (activeAbort === controller) {
          activeAbort = null;
        }
        fetchFailures = 0;
        var stampChanged = absorb(payload, unpacked.packed);
        askAgainForIncomplete(payload.nodes || {});
        var stillMissing = renderFromCache(neededNodes(undefined, true));
        // A new stamp dropped changed (or all) cached nodes, and a view
        // needing more than one request's worth of nodes has more to
        // fetch -- either way, go again (nodes already fetched are
        // cached by now).
        if (stampChanged || stillMissing.length) {
          scheduleFetch(false);
        } else {
          schedulePrefetch();
        }
      })
      .catch(function (err) {
        if (err && err.name === "AbortError") {
          return;
        }
        // A failed fetch leaves the currently-drawn content in place. A
        // database busy with a fill answers 5xx for a while (PERF.64), so
        // ask again after a growing pause (2 s, 4 s ... 30 s) instead of
        // waiting for the next camera move; a move still fetches at once.
        if (activeAbort === controller) {
          activeAbort = null;
        }
        fetchFailures += 1;
        var pause = Math.min(30000, 2000 * Math.pow(2, Math.min(fetchFailures, 5) - 1));
        if (fetchTimer) {
          clearTimeout(fetchTimer);
        }
        fetchTimer = setTimeout(doFetch, pause);
      });
  }
  var fetchFailures = 0;

  // MAP.315: the node of one block, from memory, the stored nodes or the server (one key, like the prefetch), for the
  // region an infobox names: the opening view and the view's cut rarely hold a 9-block while the camera is deep in it.
  function nodeOf(key) {
    var held = getNode(key);
    if (held !== undefined) {
      return Promise.resolve(held);
    }
    return loadStoredNodes([key]).then(function () {
      var stored = getNode(key);
      if (stored !== undefined) {
        return stored;
      }
      var params = new URLSearchParams({ keys: key });
      if (currentStamp) {
        params.set("stamp", currentStamp);
      }
      return fetch(data.fetchPath + "?" + params.toString()).then(function (response) {
        if (!response.ok) {
          throw new Error(data.fetchPath + " returned " + response.status);
        }
        return response.arrayBuffer();
      }).then(function (body) {
        var unpacked = unpackResponse(body);
        if (absorb(Object.assign({ nodes: unpacked.nodes }, unpacked.header), unpacked.packed)) {
          scheduleFetch(false);
        }
        return getNode(key);
      });
    });
  }

  // The charted region holding the point (x, y, z) pc: its 9-block with the name and the family the node carries, and
  // the 81-block of that family; null when the block has no name (nothing generated in it, API.35).
  function regionAt(x, y, z) {
    var block = cellAt(9, x, y, z, edgePc);
    return nodeOf(block.m + "." + block.shell + "." + block.band + "." + block.wedge).then(function (node) {
      if (!node || !node.name) {
        return null;
      }
      return {
        name: node.name, parent: node.parent || null, block: block, extent: node.extent || null,
        parentBlock: cellAt(81, x, y, z, edgePc),
      };
    });
  }

  // PERF.76: a node whose slowest piece timed out arrives marked `incomplete`
  // (the map shows what it has). Its pieces are fetched again after a growing
  // pause (5 s, 10 s ... 30 s), up to INCOMPLETE_ROUNDS times, while the
  // database catches up.
  var INCOMPLETE_ROUNDS = 6;
  var incompleteRounds = 0;
  var incompleteTimer = null;
  function askAgainForIncomplete(nodes) {
    var partial = Object.keys(nodes).filter(function (key) {
      return nodes[key] && nodes[key].incomplete;
    });
    if (incompleteTimer) {
      clearTimeout(incompleteTimer);
      incompleteTimer = null;
    }
    if (!partial.length) {
      incompleteRounds = 0;
      return;
    }
    if (incompleteRounds >= INCOMPLETE_ROUNDS) {
      return;
    }
    incompleteRounds += 1;
    incompleteTimer = setTimeout(function () {
      incompleteTimer = null;
      fetchFromServer(partial);
    }, Math.min(30000, 5000 * incompleteRounds));
  }

  // The nodes one pick in (about PREFETCH_ZOOM times closer), fetched at low
  // priority once the view's own are in and the camera has been still for
  // PREFETCH_IDLE_MS (MAP.158). Zooming out needs a coarser, smaller view
  // the page mostly holds already, so only the way in is fetched ahead.
  // Not on a metered or slow connection, and no more than
  // PREFETCH_MAX_BYTES a page visit. One prefetch runs at a time and a real
  // fetch never waits on or cancels it (or the other way round); its nodes
  // just land in the cache.
  var prefetching = false;
  var prefetchTimer = null;
  var prefetchedBytes = 0;
  var PREFETCH_ZOOM = 3;
  var PREFETCH_IDLE_MS = 1500;
  var PREFETCH_MAX_BYTES = 1536 * 1024;
  var PREFETCH_MAX_NODES = 12;

  function prefetchAllowed() {
    var connection = typeof navigator !== "undefined" ? navigator.connection : null;
    if (connection && (connection.saveData || /^(slow-2g|2g|3g)$/.test(connection.effectiveType || ""))) {
      return false;
    }
    return prefetchedBytes < PREFETCH_MAX_BYTES;
  }

  function cancelPrefetch() {
    if (prefetchTimer) {
      clearTimeout(prefetchTimer);
      prefetchTimer = null;
    }
  }

  function schedulePrefetch() {
    cancelPrefetch();
    if (prefetching || activeAbort || !prefetchAllowed()) {
      return;
    }
    prefetchTimer = setTimeout(prefetchNodes, PREFETCH_IDLE_MS);
  }

  function prefetchNodes() {
    prefetchTimer = null;
    if (prefetching || activeAbort || !prefetchAllowed()) {
      return;
    }
    var missing = neededNodes(Math.max(MIN_RADIUS, Math.min(MAX_RADIUS, orbit.radius / PREFETCH_ZOOM)), true).keys
      .filter(function (key) { return getNode(key) === undefined; });
    if (!missing.length) {
      return;
    }
    prefetching = true;
    loadStoredNodes(missing).then(function () {
      missing = missing.filter(function (key) { return getNode(key) === undefined; });
      if (!missing.length || fetchTimer || activeAbort) {
        prefetching = false;
        return undefined;
      }
      var params = new URLSearchParams({ keys: missing.slice(0, PREFETCH_MAX_NODES).join(",") });
      if (currentStamp) {
        params.set("stamp", currentStamp);
      }
      return fetch(data.fetchPath + "?" + params.toString(), { priority: "low", headers: { "X-Map-Prefetch": "1" } })
        .then(function (response) {
          if (!response.ok) {
            throw new Error(data.fetchPath + " returned " + response.status);
          }
          return response.arrayBuffer();
        })
        .then(function (body) {
          var unpacked = unpackResponse(body);
          prefetching = false;
          prefetchedBytes += body.byteLength;
          // A new stamp means the view's own nodes may be stale too.
          if (absorb(Object.assign({ nodes: unpacked.nodes }, unpacked.header), unpacked.packed)) {
            scheduleFetch(false);
          } else {
            schedulePrefetch();
          }
        });
    }).catch(function () {
      // Only a head start: the real fetch gets these when they're needed.
      prefetching = false;
    });
  }

  // Normally every node of the first frame is already cached (the page embeds them), so this
  // only fetches if the page's own node fetch came up short.
  // First frame's nodes, then (normally nothing is missing) the prefetch.
  renderFromCache(neededNodes(undefined, true));
  scheduleFetch(false);

  // --- Territories -----------------------------------------------------------
  //
  // The Territories button draws who holds what (docs/design/habitability/population-
  // and-politics.md): each polity's reach as a translucent ball around
  // its capital, and the systems it owns as small dots in its own color
  // (/galaxy/territories, at most 20,000 nearest their capitals). Both
  // are cheap at galaxy scale -- a handful of balls, one points object --
  // and the dots keep a fixed size on screen, so zooming out never turns
  // them into a smear. Fetched once, the first time the button is
  // pressed, and listed under the map with each polity's system count.
  // MAP.240: the cube, in parsecs, a layer of many scattered points is split into so each part has an origin near it
  // (static/camerarelative.js); half its diagonal is under SAFE_EXTENT_PC.
  var RELATIVE_CELL_PC = 1000;
  var TERRITORY_DOT_PX = 3.5;
  var TERRITORY_BALL_OPACITY = 0.1;
  var territoryGroup = null;
  var territoryLoading = false;
  var territoryWanted = false;

  function lyToPc(ly) {
    return convert(ly, "ly", "pc");
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
      var flat = [];
      var tints = [];
      points.forEach(function (point) {
        flat.push(point.x, point.y, point.z);
        var color = territoryColor(point.color);
        tints.push(color.r, color.g, color.b);
      });
      // MAP.240: the dots in cubes of their own, each an object at its cube's middle with offsets from it.
      groupByCell(flat, RELATIVE_CELL_PC).forEach(function (chunk) {
        var colors = new Float32Array(chunk.indices.length * 3);
        chunk.indices.forEach(function (at, n) {
          colors[3 * n] = tints[3 * at];
          colors[3 * n + 1] = tints[3 * at + 1];
          colors[3 * n + 2] = tints[3 * at + 2];
        });
        var geometry = new THREE.BufferGeometry();
        geometry.setAttribute("position", new THREE.BufferAttribute(chunk.positions, 3));
        geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
        var dots = new THREE.Points(geometry, new THREE.PointsMaterial({
          size: TERRITORY_DOT_PX * renderer.getPixelRatio(), sizeAttenuation: false, vertexColors: true,
          transparent: true, depthWrite: false,
        }));
        dots.position.set(chunk.origin[0], chunk.origin[1], chunk.origin[2]);
        dots.userData.relativeLayer = "territory dots";
        dots.frustumCulled = false;
        group.add(dots);
      });
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
  // A line through `points` (Vector3s of doubles) as a group of short THREE.Lines (MAP.240): each piece is at most
  // LINE_PIECE_PC long, an object at its own middle with offsets from it, so a course across the galaxy stays exact too.
  // `dashed` has each piece's line distances worked out (a dash pattern starts again on each piece).
  var LINE_PIECE_PC = 4000;
  function lineFrom(points, material, dashed) {
    var group = new THREE.Group();
    linePieces(points.map(function (point) { return [point.x, point.y, point.z]; }), LINE_PIECE_PC).forEach(function (built) {
      var geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(built.positions, 3));
      var line = new THREE.Line(geometry, material);
      line.position.set(built.origin[0], built.origin[1], built.origin[2]);
      line.userData.relativeLayer = "course line";
      if (dashed) line.computeLineDistances();
      group.add(line);
    });
    return group;
  }
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
    // MAP.240: a line is an object at its first stop with offsets from it.
    courseGroup.add(lineFrom(linePoints, new THREE.LineBasicMaterial({ color: new THREE.Color(accentColor), depthTest: false, transparent: true })));
    // NAV.20: with stops on the way, the straight line from the start to the
    // end is drawn apart from the route, dashed and muted.
    if (course.direct && course.direct.length === 2) {
      var directPoints = course.direct.map(function (point) { return new THREE.Vector3(point.x, point.y, point.z); });
      courseGroup.add(lineFrom(directPoints, new THREE.LineDashedMaterial({
        color: new THREE.Color(0x9aa3b8), depthTest: false, transparent: true,
        dashSize: directPoints[0].distanceTo(directPoints[1]) / 60,
        gapSize: directPoints[0].distanceTo(directPoints[1]) / 120,
      }), true));
    }
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

  // --- Superbubbles (GEN.254) --------------------------------------------------------
  //
  // The cavities the supernovae of OB associations blow in the gas (the Local Bubble
  // is one), from the galaxy's seed: the page embeds each one's centre and radius in
  // galaxy-frame parsecs. Each is drawn as its sphere (MAP.236, superbubble.js): a faint fill and a
  // crisp rim where the sphere meets the line of sight, under everything else and never picked. No
  // molecular cloud lies inside one.
  var bubbles = data.superbubbles || [];
  var bubbleLook = superbubble.bubbleLook(isLightBackground());
  var bubbleMaterial = null;
  var bubbleMeshes = [];
  if (bubbles.length) {
    // MAP.240: one object for each bubble, at its centre, so no vertex is more than the radius from its origin.
    bubbleMaterial = new THREE.ShaderMaterial({
      uniforms: {
        color: { value: new THREE.Vector3().fromArray(superbubble.colorComponents(bubbleLook.color)) },
        fillAlpha: { value: bubbleLook.fillAlpha },
        rimAlpha: { value: bubbleLook.rimAlpha },
        rimPx: { value: superbubble.RIM_PX },
      },
      vertexShader: superbubble.VERTEX_SHADER,
      fragmentShader: superbubble.FRAGMENT_SHADER,
      transparent: true, depthWrite: false, depthTest: false, side: THREE.BackSide,
    });
    bubbles.forEach(function (bubble) {
      var bubbleMesh = new THREE.Mesh(new THREE.SphereGeometry(bubble[3], 48, 24), bubbleMaterial);
      bubbleMesh.position.set(bubble[0], bubble[1], bubble[2]);
      bubbleMesh.renderOrder = 0;
      bubbleMesh.frustumCulled = false;
      bubbleMesh.name = "superbubbles";
      bubbleMesh.userData.relativeLayer = "superbubbles";
      bubbleMesh.visible = false;    // off until the Menu's Superbubbles button turns them on (registered with the other layers)
      bubbleMeshes.push(bubbleMesh);
      scene.add(bubbleMesh);
    });
    renderer.domElement.dataset.superbubbles = String(bubbles.length);
  }

  // For the browser tests (MAP.236): the look the bubbles are drawn with.
  canvasEl.galaxyBubbleLook = function () {
    return bubbleMaterial ? {
      count: bubbles.length, shown: bubbleMeshes.some(function (mesh) { return mesh.visible; }), color: bubbleLook.color, fillAlpha: bubbleMaterial.uniforms.fillAlpha.value,
      rimAlpha: bubbleMaterial.uniforms.rimAlpha.value, rimPx: bubbleMaterial.uniforms.rimPx.value,
    } : null;
  };

  // --- Sectors a generate job made (ADM.31) ---------------------------------------
  //
  // With ?made=<since>,<until> (the job page's "Show on Galaxy Map") the page
  // embeds the centres and addresses of the sectors that run made, up to a
  // cap. They are ringed in one points object, over everything else and never
  // picked, and the map opens on the smallest stage holding them (the same
  // fit a NAV course gets), unless the URL names a stage itself.
  var made = data.made && (data.made.points || []).length ? data.made : null;
  var MADE_RING_PX = 16;
  var MADE_FIT_LIMIT = 300;
  if (made) {
    var madeFlat = made.points.reduce(function (all, point) { all.push(point.x, point.y, point.z); return all; }, []);
    // MAP.240: the rings in cubes of their own, each an object at its cube's middle.
    groupByCell(madeFlat, RELATIVE_CELL_PC).forEach(function (chunk) {
      var madeGeometry = new THREE.BufferGeometry();
      madeGeometry.setAttribute("position", new THREE.BufferAttribute(chunk.positions, 3));
      var madePoints = new THREE.Points(madeGeometry, new THREE.PointsMaterial({
        map: highlightTexture, size: MADE_RING_PX, sizeAttenuation: false, transparent: true,
        depthWrite: false, depthTest: false, alphaTest: 0.05,
      }));
      madePoints.position.set(chunk.origin[0], chunk.origin[1], chunk.origin[2]);
      madePoints.userData.relativeLayer = "made rings";
      madePoints.renderOrder = 4;
      madePoints.frustumCulled = false;
      scene.add(madePoints);
    });
    // A stage that holds the extremes and an even sample holds them all.
    var madeFit = made.sectors;
    if (madeFit.length > MADE_FIT_LIMIT) {
      var picked = {};
      [0, 1, 2].forEach(function (axis) {
        var low = madeFit[0], high = madeFit[0];
        madeFit.forEach(function (sector) {
          if (sector[axis] < low[axis]) low = sector;
          if (sector[axis] > high[axis]) high = sector;
        });
        picked[low.join(".")] = low;
        picked[high.join(".")] = high;
      });
      var step = madeFit.length / (MADE_FIT_LIMIT - 6);
      for (var at = 0; at < madeFit.length; at += step) {
        var sample = madeFit[Math.floor(at)];
        picked[sample.join(".")] = sample;
      }
      madeFit = Object.keys(picked).map(function (key) { return picked[key]; });
    }
    if (!courseSectors.length) {
      courseSectors = madeFit.map(function (sector) {
        return { shell: sector[0], band: sector[1], wedge: sector[2] };
      });
    }
  }

  // NAV.20: a click on a ringed stop opens its page. Registered ahead of the
  // stage view's own pointer handlers (which come later) and ending the event
  // for them, so a click on a ring never drills into the sector beneath it.
  if (course) {
    var ringDown = null;
    var ringVector = new THREE.Vector3();
    var ringAt = function (event) {
      var rect = canvasEl.getBoundingClientRect();
      var found = null;
      var best = COURSE_RING_PX * 0.9;
      course.points.forEach(function (point) {
        ringVector.set(point.x, point.y, point.z).project(camera);
        if (ringVector.z > 1) return;
        var dx = (ringVector.x + 1) / 2 * rect.width + rect.left - event.clientX;
        var dy = (1 - ringVector.y) / 2 * rect.height + rect.top - event.clientY;
        var distance = Math.sqrt(dx * dx + dy * dy);
        if (distance < best) { best = distance; found = point; }
      });
      return found;
    };
    canvasEl.addEventListener("pointerdown", function (event) {
      ringDown = event.button === 0 && ringAt(event) ? { x: event.clientX, y: event.clientY } : null;
      if (ringDown) event.stopImmediatePropagation();
    });
    canvasEl.addEventListener("pointerup", function (event) {
      if (!ringDown) return;
      var moved = Math.abs(event.clientX - ringDown.x) + Math.abs(event.clientY - ringDown.y);
      ringDown = null;
      event.stopImmediatePropagation();
      var point = moved < 6 ? ringAt(event) : null;
      if (point && point.url) window.location.href = point.url;
    });
    canvasEl.addEventListener("pointermove", function (event) {
      canvasEl.style.cursor = ringAt(event) ? "pointer" : "";
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
  // The map opens on the drill-down (galaxystageview.js): the whole
  // galaxy in 3D, a window of level-27 blocks around the block clicked, and
  // so on down to a window of sectors. Pointer and key input goes to the
  // stage view; nodes (stars, clouds) follow the camera.
  var picker = createPicker(camera, canvasEl);

  // The sector page's Contents "Show on map" buttons (MAP.46) name a
  // cloud by its key; selecting it shows its details and its ring, and
  // brings the map into view. A button shows once its entry is on the map.
  var mapTargetsOn = false;
  function showMapTargets(on) {
    mapTargetsOn = on;
    document.querySelectorAll("[data-map-target]").forEach(function (button) {
      button.hidden = !(on && stageView.entryByKey(button.dataset.mapTarget));
    });
  }
  // The Menu's "Show on the map" (MAP.79): one toggle per kind of object in
  // the open sector (all on to begin with), highlighted while shown.
  var kindsEl = document.getElementById("galaxymap3d-kinds");
  // MAP.217: the buttons are made from the layer registry the maps share (maplayers.js); it leaves out the
  // kinds the open sector has none of.
  var kindsLayersEl = document.createElement("span");
  kindsLayersEl.className = "galaxy-layers";
  var drawKindLayers = null;
  function rebuildKinds() {
    if (!kindsEl) return;
    var kinds = stageView.kinds();
    kindsEl.textContent = "";
    refreshStarFilters();
    // With no sector open the star filters still apply to the stars drawn
    // at galaxy scale, so they stay.
    kindsEl.hidden = false;
    kindsEl.appendChild(kindsLayersEl);
    if (drawKindLayers) drawKindLayers();
    else drawKindLayers = stageView.layers.renderMenu(kindsLayersEl);
    if (!kinds.length) {
      rebuildStarFilters();
      return;
    }
    rebuildHighlights(kinds);
    rebuildStarFilters();
  }

  // MAP.123: a toggle per kind of phenomenon that draws it larger and brighter.
  function rebuildHighlights(kinds) {
    const phenomena = kinds.filter(function (kind) { return kind.kind !== "star"; });
    if (!phenomena.length) return;
    const heading = document.createElement("span");
    heading.className = "galaxy-kinds-heading";
    heading.textContent = "Highlight";
    kindsEl.appendChild(heading);
    phenomena.forEach(function (kind) {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "starmap-btn starmap-toggle";
      button.dataset.highlightKind = kind.kind;
      button.textContent = kind.label;
      button.title = "Draw the " + kind.label.toLowerCase() + " larger and brighter";
      button.setAttribute("aria-pressed", String(stageView.kindMarked(kind.kind)));
      button.addEventListener("click", function () {
        stageView.setKindMarked(kind.kind, button.getAttribute("aria-pressed") !== "true");
        button.setAttribute("aria-pressed", String(stageView.kindMarked(kind.kind)));
      });
      kindsEl.appendChild(button);
    });
  }

  // MAP.123: the star classes shown (a toggle each, all on to begin with) and
  // the dimmest star shown (a slider on a log scale: 0 is every star).
  // The slider spans the dimmest to the brightest star at the current zoom
  // (an open sector's stars, else the loaded nodes'), on a log scale:
  // step 1 is the dimmest star there is, the last step the brightest.
  var LUM_STEPS = 100;
  var lumSliderEl = null;
  var lumTextEl = null;
  function lumRange() {
    var range = (typeof stageView !== "undefined" && stageView && stageView.sectorLuminosityRange()) || nodeLumRange;
    return range || [1e-5, 1e6];
  }
  function lumFromStep(step) {
    if (step <= 0) return 0;
    var range = lumRange();
    if (step >= LUM_STEPS || !(range[1] > range[0])) return step >= LUM_STEPS ? range[1] : range[0];
    var share = (step - 1) / (LUM_STEPS - 1);
    return range[0] * Math.pow(range[1] / range[0], share);
  }
  function stepFromLum(lum) {
    if (!(lum > 0)) return 0;
    var range = lumRange();
    if (!(range[1] > range[0])) return 1;
    var share = Math.log(lum / range[0]) / Math.log(range[1] / range[0]);
    return Math.max(1, Math.min(LUM_STEPS, 1 + Math.round(share * (LUM_STEPS - 1))));
  }
  // Moves the slider to where the chosen floor sits on the new scale.
  function updateLumSlider() {
    if (!lumSliderEl || !lumSliderEl.isConnected) return;
    lumSliderEl.value = String(stepFromLum(stageView.minLuminosity()));
    if (lumTextEl) lumTextEl.textContent = "Dimmest star shown: " + floorText(stageView.minLuminosity());
  }
  // MAP.166: "every star" only when the view holds every star there is, which
  // only an open sector's scene does. A galaxy view's nodes are capped and
  // budgeted (and sectors may be unmade), so with the floor at 0 it names
  // the dimmest star it actually carries.
  function viewIsComplete() {
    return !!(typeof stageView !== "undefined" && stageView && stageView.sectorLuminosityRange());
  }
  function floorText(lum) {
    if (lum > 0) return formatLum(lum);
    if (viewIsComplete()) return "every star";
    return nodeLumRange ? formatLum(nodeLumRange[0]) : "no stars in view";
  }
  function formatLum(lum) {
    if (!(lum > 0)) return "every star";
    return formatLuminositySol(lum) + " and up";
  }
  function rebuildStarFilters() {
    var classes = stageView.starClasses();
    if (!classes.length) {
      classes = STAR_FILTER_CLASSES.map(function (c) { return { starClass: c, count: null }; });
    }
    var heading = document.createElement("span");
    heading.className = "galaxy-kinds-heading";
    heading.textContent = "Star classes";
    kindsEl.appendChild(heading);
    classes.forEach(function (item) {
      var button = document.createElement("button");
      button.type = "button";
      button.className = "starmap-btn starmap-toggle";
      button.dataset.starClass = item.starClass;
      button.textContent = item.starClass === "other" ? "Other" : item.starClass;
      button.title = "Show or hide the class " + button.textContent + " stars" + (item.count == null ? "" : " (" + item.count + ")");
      button.setAttribute("aria-pressed", String(!stageView.starClassHidden(item.starClass)));
      button.addEventListener("click", function () {
        stageView.setStarClassHidden(item.starClass, button.getAttribute("aria-pressed") === "true");
        button.setAttribute("aria-pressed", String(!stageView.starClassHidden(item.starClass)));
        refreshStarFilters();
      });
      kindsEl.appendChild(button);
    });
    var label = document.createElement("label");
    label.className = "galaxy-lum";
    var text = document.createElement("span");
    var slider = document.createElement("input");
    slider.type = "range";
    slider.min = "0";
    slider.max = String(LUM_STEPS);
    slider.step = "1";
    slider.id = "galaxymap3d-lum";
    slider.value = String(stepFromLum(stageView.minLuminosity()));
    lumSliderEl = slider;
    lumTextEl = text;
    var show = function () { text.textContent = "Dimmest star shown: " + floorText(lumFromStep(Number(slider.value))); };
    slider.addEventListener("input", function () {
      show();
      stageView.setMinLuminosity(lumFromStep(Number(slider.value)));
      refreshStarFilters();
    });
    show();
    label.appendChild(text);
    label.appendChild(slider);
    kindsEl.appendChild(label);
  }

  // The Contents table draws its rows as they scroll (datatable.js), so the
  // buttons are found at the click, and shown again after each redraw.
  document.addEventListener("click", function (event) {
    var button = event.target.closest && event.target.closest("[data-map-target]");
    if (!button) return;
    var entry = stageView.entryByKey(button.dataset.mapTarget);
    if (!entry) return;
    stageView.selectEntry(entry);
    if (viewport) viewport.scrollIntoView({ block: "center" });
    canvasEl.focus({ preventScroll: true });
  });
  document.addEventListener("datatable:rows", function () { showMapTargets(mapTargetsOn); });

  // The legend of a single-statistic color mode (galaxyblocks' legendFor), or
  // none for the default.
  function showLegend(legend) {
    var el = document.getElementById("galaxymap3d-legend");
    if (!el) return;
    el.hidden = !legend;
    if (!legend) return;
    el.querySelector(".galaxy-legend-title").textContent = legend.label + " (" + legend.unit + ")";
    el.querySelector(".galaxy-legend-bar").style.background = "linear-gradient(to right, " + legend.stops.join(", ") + ")";
    el.querySelector(".galaxy-legend-low").textContent = legend.lowText;
    el.querySelector(".galaxy-legend-high").textContent = legend.highText;
    el.setAttribute("aria-label", legend.label + ": " + legend.lowText + " (dark blue) to " + legend.highText + " (yellow), " + legend.unit);
  }

  var navPad = null;
  var stageView = createStageView({
    THREE: THREE, scene: scene, camera: camera, canvasEl: canvasEl,
    edgePc: edgePc, shape: galaxyShape, galaxyRadius: GALAXY_RADIUS, reducedMotion: reducedMotion,
    accentColor: accentColor, canGenerate: !!data.generate, courseSectors: courseSectors,
    blockScene: localBlockScene,
    makeBlockMesh: makeBlockMesh,
    // MAP.131: the blocks' color modes, the legend of the chosen one and the
    // select that follows the address.
    colorModes: Object.keys(COLOR_MODES),
    showLegend: showLegend,
    colorChanged: function (mode) {
      var select = document.getElementById("galaxymap3d-color-by");
      if (select) select.value = mode;
    },
    // While picking a NAV start or destination, the map's own URLs keep
    // the pick, so Back, a reload or a bookmark of the page keep it too.
    pickQuery: function () { return navPick.query(); },
    navPick: navPick,
    // The stage view turns the camera freely (MAP.96): `quaternion` is
    // its turn, the camera sitting `dist` back along it from `target`.
    // `focusRadius` (optional) is how wide the focus is and `inside` (the
    // center of the block the camera is inside, or null) the block drawn
    // from the inside (MAP.149).
    nearFieldAt: nearFieldAt,
    // MAP.301: the stars the map draws now (their x, y, z in parsecs), which the wheel reads to head for a crowd of them.
    drawnStars: function () { return starList; },
    // MAP.188: the blocks' alpha share (see gridUniform).
    setGridScale: function (scale) { gridUniform.value = scale; },
    gridScale: function () { return gridUniform.value; },
    setCamera: function (v) {
      target.set(v.target[0], v.target[1], v.target[2]);
      orbit.radius = v.dist;
      nearUniforms.nearD.value = v.dist;
      nearUniforms.nearR.value = v.focusRadius || 0;
      if (v.inside) nearUniforms.insideCenter.value.set(v.inside[0], v.inside[1], v.inside[2]);
      else nearUniforms.insideCenter.value.set(NO_BLOCK_INSIDE, NO_BLOCK_INSIDE, NO_BLOCK_INSIDE);
      camera.quaternion.fromArray(v.quaternion);
      var back = new THREE.Vector3(0, 0, v.dist).applyQuaternion(camera.quaternion);
      eyeBack.copy(back);
      camera.position.copy(target).add(back);
      // A system opened in place is zoomed to a moon, far inside the usual
      // near plane (the depth buffer is logarithmic, so a near plane close
      // to the camera costs no precision).
      var near = Math.min(BASE_NEAR, v.dist * 1e-3);
      if (near !== camera.near) {
        camera.near = near;
        camera.updateProjectionMatrix();
      }
      camera.updateMatrixWorld();
      updateScaleBar();
      if (navPad) navPad.sync();
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
    fetchWindow: function (query) {
      return fetch((data.windowPath || "/galaxy/window") + query, { headers: { Accept: "application/json" } })
        .then(function (response) {
          if (!response.ok) {
            throw new Error("HTTP " + response.status);
          }
          return response.json();
        });
    },
    // MAP.300: the stars the map draws now, from the nodes' lists (the brightest last): what the wheel can zoom toward when no sector is open.
    drawnStars: function () {
      return pickableStars();
    },
    // MAP.335: the generated sectors near a point (the eye), asked for as the camera comes near them.
    fetchNear: function (query) {
      return fetch((data.nearPath || "/galaxy/near") + query, { headers: { Accept: "application/json" } })
        .then(function (response) {
          if (!response.ok) {
            throw new Error("HTTP " + response.status);
          }
          return response.json();
        });
    },
    showBlockInfo: function (info) { showBlockInfo(info, edgePc); },
    showPlacedInfo: showPlacedInfo,
    showCellInfo: showCellInfo,
    showHint: function (text) { showHint(text, false); },
    // MAP.298: the names over the map are clicked through here, not as elements of their own.
    labelAt: function (clientX, clientY) {
      if (!galaxyLabels) return null;
      var box = canvasEl.getBoundingClientRect();
      return galaxyLabels.at(clientX - box.left, clientY - box.top);
    },
    labelHot: function (id) { if (galaxyLabels) galaxyLabels.setHot(id); },
    labelAct: function (label) { actOnLabel(label); },
    picker: picker,
    clearSelection: function () { selectionRing.hide(); },
    restingHint: restingHint,
    sectorUrl: function (id) { return sceneData.sectorUrl ? sectorUrl(id) : null; },
    // A sector page's map is locked to its sector (MAP.68): {shell, band,
    // wedge}; another sector's click goes to that sector's page.
    pinned: data.pinned || null,
    openSectorPage: function (id) { if (sceneData.sectorUrl) location.assign(sectorUrl(id)); },
    // The sector page's scene JSON (planetgen/web/sector_page.py).
    // A system of the open sector opened in place (MAP.125): its page's /scene.
    fetchSystemScene: sceneData.sectorUrl ? function (href) {
      return fetch(href.split("?")[0] + "/scene", { headers: { Accept: "application/json" } })
        .then(function (response) {
          if (!response.ok) {
            throw new Error("HTTP " + response.status);
          }
          return response.json();
        });
    } : null,
    fetchSectorScene: sceneData.sectorUrl ? function (id) {
      var path = sectorUrl(id).split("?")[0] + "/scene";
      return fetch(path, { headers: { Accept: "application/json" } })
        .then(function (response) {
          if (!response.ok) {
            throw new Error("HTTP " + response.status);
          }
          return response.json();
        });
    } : null,
    // The scene of a cell nothing was generated in (MAP.162), from its address.
    fetchUnchartedScene: sceneData.unchartedSceneUrl ? function (at) {
      var path = String(sceneData.unchartedSceneUrl).replace("{shell}", encodeURIComponent(at.shell))
        .replace("{band}", encodeURIComponent(at.band)).replace("{wedge}", encodeURIComponent(at.wedge));
      return fetch(path, { headers: { Accept: "application/json" } })
        .then(function (response) {
          if (!response.ok) {
            throw new Error("HTTP " + response.status);
          }
          return response.json();
        });
    } : null,
    lightBackground: isLightBackground(),
    pixelRatio: function () { return renderer.getPixelRatio(); },
    showInfo: showInfo,
    hideCell: hideCell,
    // MAP.348: what the opened sector's things say first when hovered.
    preview: function (facts, entry) { previewHover(entry, facts, entry); },
    viewport: viewport,
    // The sector page's "Show on map" buttons (MAP.46) are there once the
    // sector is.
    opened: function () { showMapTargets(true); rebuildKinds(); },
    closed: function () { showMapTargets(false); rebuildKinds(); },
    kindsChanged: function () { rebuildKinds(); showMapTargets(true); },
    // MAP.336: the charted regions (the find-a-place box lists them by name) and where a Go to a region lands (the reader's setting).
    regions: function () { return regionBlobList || []; },
    regionLanding: function () {
      var storage = null;
      try { storage = window.localStorage; } catch (error) { storage = null; }
      return regionLanding(storage);
    },
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
      where: document.getElementById("galaxymap3d-where"),
      pageCrumb: document.querySelector("nav.breadcrumbs"),
      tooltip: document.getElementById("galaxymap3d-tooltip"),
      notice: document.getElementById("galaxymap3d-notice"),
      address: document.getElementById("galaxymap3d-address"),
      matches: document.getElementById("galaxymap3d-matches"),
      controls: document.getElementById("galaxymap3d-controls"),
      steps: document.getElementById("galaxymap3d-steps"),
    },
  });
  rebuildKinds();
  // MAP.217: the layers' legend lines, over the map; the browser tests register a layer on the canvas's registry.
  var layerLegendEl = document.getElementById("galaxymap3d-layer-legend");
  if (layerLegendEl) stageView.layers.renderLegend(layerLegendEl, { collapsible: true });
  // The wheel over the Key's fold zooms the map as it does over the other controls (MAP.234).
  if (layerLegendEl) layerLegendEl.addEventListener("wheel", function (event) { if (stageView.onWheel(event)) event.preventDefault(); }, { passive: false });
  canvasEl.__mapLayers = stageView.layers;
  // MAP.219: the marker layer (galaxymarkers.js): glyphs over the notable things, a toggle and a legend line per kind.
  var markerKinds = createMarkerKinds();
  var POI_GROUP = "Points of interest";
  var CLASS_GROUP = "World classes";   // MAP.220: one optional marker for each planet class, all off until chosen
  var LIFE_CLASS_GROUP = "Life classes";   // MAP.351: one optional marker for each life class, all off until chosen
  var TECH_GROUP = "Civilization technology";   // MAP.344: one optional marker for each technology band of a civilization, all off until chosen
  var INHABITED_MARKER_KIND = INTEREST_KINDS.find(function (kind) { return kind.interest === "inhabited"; });
  var markerLayer = null;
  var markersDirty = true;
  if (!data.pinned) {
    markerKinds.register(BLACK_HOLE_KIND);
    MAP_INTEREST_KINDS.forEach(function (kind) { markerKinds.register(kind); });
    markerLayer = createMarkerLayer(scene, camera, canvasEl, isLightBackground() ? "#ffffff" : "#0b0d14", {
      listIn: viewport, onChoose: function (entry) { chooseMarker(entry); },
    });
    markerLayer.setOutsideShare(hiddenCell ? SECTOR_OUTSIDE_SHARE : 1);    // MAP.341: a sector may be open before the layer exists
    stageView.layers.bulkGroup(POI_GROUP);
    stageView.layers.bulkGroup(CLASS_GROUP);
    stageView.layers.bulkGroup(LIFE_CLASS_GROUP);
    stageView.layers.bulkGroup(TECH_GROUP);
    // The Menu lists a group's layers together, so the technology bands (registered ahead of the Civilizations kind, which they pre-empt) come
    // after the points of interest.
    var familyOrder = { tech: 1, "class": 2, lifeclass: 3 };
    markerKinds.kinds().map(function (kind, index) { return { kind: kind, index: index }; }).sort(function (p, q) {
      return (familyOrder[p.kind.family] || 0) - (familyOrder[q.kind.family] || 0) || p.index - q.index;
    }).map(function (entry) { return entry.kind; }).forEach(function (kind) {
      stageView.layers.register({
        id: kind.id, label: kind.label, modes: ["galaxy"], forward: !!kind.interest && !kind.family,
        group: kind.family === "class" ? CLASS_GROUP : kind.family === "lifeclass" ? LIFE_CLASS_GROUP : kind.family === "tech" ? TECH_GROUP : kind.interest ? POI_GROUP : "Markers",
        data: { marker: kind.glyph || kind.icon },
        defaultOn: kind.defaultOn, urlKey: "hide", urlKeyOn: "show", fadedWithSectors: false,
        title: "Show or hide the " + kind.label.toLowerCase() + " markers", legend: markerLegendOf(kind),
      });
    });
    // The registry forgets ids nobody registered when it read the address, so read it again for these (without the
    // stage writing the address back before it has taken its marks from it).
    stageView.readLayersFromLocation();
    stageView.layers.refresh();   // the Menu and legend were drawn before these kinds existed
    // MAP.306: the point-of-interest rows open as the visitor left them unless the address chose.
    var poiIds = markerKinds.kinds().filter(function (kind) { return kind.interest; }).map(function (kind) { return kind.id; });
    var poiDefaults = {};
    poiIds.forEach(function (id) { poiDefaults[id] = markerKinds.get(id).defaultOn ? "shown" : "hidden"; });
    var poiStorage = null;
    try { poiStorage = window.localStorage; } catch (error) { poiStorage = null; }
    if (poiStorage) applyPoiPrefs(poiStorage, stageView.layers, poiIds, poiDefaults);
    stageView.layers.onChange(function (id) {
      if (markerKinds.get(id)) markersDirty = true;
      if (poiStorage && poiIds.indexOf(id) >= 0) savePoiPrefs(poiStorage, stageView.layers, poiIds);
    });
    canvasEl.__mapMarkers = markerLayer;    // the browser tests read what is drawn
  }
  // MAP.327: the frame of the sector opened in place: a locked wireframe and an icon on each wall that faces a neighbour.
  if (viewport && !data.pinned) {
    sectorFrame = createSectorFrame({
      THREE: THREE, scene: scene, camera: camera, canvasEl: canvasEl, parent: viewport, accent: accentColor, edgePc: edgePc,
      // The sector across a wall is the one the shared step helper finds (galaxyprisms.neighbourBlock, as the arrow keys use).
      neighbour: function (wall) {
        if (!hiddenCell) return null;
        var middle = cellCentre(hiddenCell);
        var here = sectorAddressAt(middle[0], middle[1], middle[2], edgePc);
        var found = neighbourBlock({ m: 1, shell: here.shell, band: here.band, wedge: here.wedge }, wall.dir, edgePc);
        return found && stageOfSector(found.shell, found.band, found.wedge) ? { shell: found.shell, band: found.band, wedge: found.wedge } : null;
      },
      label: function (sector) { return "Sector " + sectorDesignation(sector.shell, sector.band, sector.wedge); },
      cellPoints: cellEdgePoints,
      onMove: function (sector) { stageView.goTo({ sector: sector, open: true }); },
    });
    canvasEl.__mapSectorFrame = sectorFrame;    // the browser tests read it
    // For tests (MAP.327): what the opened sector does to the rest of the map.
    canvasEl.galaxySectorView = function () {
      return { outsideDim: starMaterial.uniforms.outsideDim.value, glowVisible: glow ? glow.isVisible() : null, opened: !!hiddenCell };
    };
    sectorFrame.set(hiddenCell);
    // MAP.327: in the galactic view the same locked wireframe round the charted region the view is in (the 9-block holding the
    // point the camera looks at), with no dimming and no wall icons; it is the sector's frame that takes over once a sector is opened.
    if (regionLayer) {
      regionFrame = createSectorFrame({
        THREE: THREE, scene: scene, camera: camera, canvasEl: canvasEl, parent: viewport, accent: accentColor, edgePc: edgePc,
        neighbour: function () { return null; }, label: function () { return ""; }, cellPoints: cellEdgePoints, onMove: function () {},
      });
      canvasEl.__mapRegionFrame = regionFrame;    // the browser tests read it
    }
  }
  // MAP.266: the glow is a layer of the Menu like the markers, on until the visitor says otherwise.
  if (glow) {
    stageView.layers.register({
      id: "glow", label: "Unresolved glow", modes: ["galaxy"], group: "Show on the map", data: { glow: "unresolved" },
      defaultOn: true, urlKey: "hide", urlKeyOn: "show", fadedWithSectors: false,
      title: "Show or hide the soft light of the stars too far away to tell apart",
      legend: null,
    });
    stageView.readLayersFromLocation();
    stageView.layers.refresh();
    glowLayerOn = stageView.layers.isOn("glow");
    glow.setVisible(glowLayerOn && !hiddenCell);
    stageView.layers.onChange(function (id, on) {
      if (id === "glow") {
        glowLayerOn = on;
        glow.setVisible(glowLayerOn && !hiddenCell);
      }
    });
    canvasEl.__mapGlow = glow;    // the browser tests read what is drawn
  }
  // MAP.176: the nebula regions are layers of the Menu, one for each family, all on until the visitor says otherwise, with a
  // legend line each and one line that says what a ring is.
  if (nebulaLayer) {
    var NEBULA_GROUP = "Nebulae";
    var nebulaLooksNow = NEBULA_LOOKS;
    var nebulaLightNow = isLightBackground();
    var nebulaSwatch = function (family) {
      var colour = nebulaFamilyColour(family, nebulaLooksNow, nebulaLightNow);
      return "rgb(" + colour.map(function (v) { return Math.round(255 * v); }).join(", ") + ")";
    };
    stageView.layers.bulkGroup(NEBULA_GROUP);
    NEBULA_FAMILIES.forEach(function (family) {
      stageView.layers.register({
        id: "nebulae-" + family, label: NEBULA_FAMILY_LABELS[family], modes: ["galaxy"], group: NEBULA_GROUP,
        data: { nebulae: family }, defaultOn: true, urlKey: "hide", urlKeyOn: "show", fadedWithSectors: false,
        title: "Show or hide the groups of " + NEBULA_FAMILY_LABELS[family].toLowerCase() + " on the map",
        legend: { text: NEBULA_FAMILY_LABELS[family], swatch: nebulaSwatch(family) },
      });
    });
    stageView.layers.register({
      id: "nebulae-groups", label: "Nebula groups", modes: ["galaxy"], group: NEBULA_GROUP, defaultOn: true, fadedWithSectors: false,
      available: function () { return false; },
      legend: function () {
        var any = NEBULA_FAMILIES.some(function (family) { return nebulaFamilyOn[family] !== false; });
        return any && nebulaLayer.count() > 0 ? { text: "Ring: nebulae grouped when close together; zoom in to separate", swatchClass: "map-nebula-group-swatch" } : null;
      },
    });
    stageView.readLayersFromLocation();
    stageView.layers.refresh();
    NEBULA_FAMILIES.forEach(function (family) { nebulaFamilyOn[family] = stageView.layers.isOn("nebulae-" + family); });
    setNebulaRegions();
    stageView.layers.onChange(function (id, on) {
      if (id.indexOf("nebulae-") !== 0) return;
      var family = id.slice("nebulae-".length);
      if (nebulaFamilyOn[family] === undefined) return;
      nebulaFamilyOn[family] = on;
      setNebulaRegions();
      stageView.layers.refresh();
    });
    canvasEl.__mapNebulae = nebulaLayer;    // the browser tests read what is drawn
  }
  // MAP.328: the stars and bodies of a system opened in the sector turn slowly; this switch (a standing preference, as the System
  // Map's) holds them still. The galactic views turn nothing, so it matters only inside a system.
  stageView.layers.register({
    id: "spin", label: "Spin", modes: ["galaxy"], group: "Show on the map", data: { setting: "spin" },
    defaultOn: stageView.spinSwitch.on(), fadedWithSectors: false,
    title: "Let the stars and bodies of an opened system turn slowly on their axes, or hold them still",
    legend: null, keyless: "a setting: it draws nothing of its own",
  });
  stageView.layers.refresh();
  stageView.layers.onChange(function (id, on) {
    if (id === "spin") stageView.spinSwitch.set(on);
  });
  // The superbubbles are a layer of the Menu, off unless the visitor shows them.
  if (bubbleMeshes.length) {
    stageView.layers.register({
      id: "superbubbles", label: "Superbubbles", modes: ["galaxy"], group: "Show on the map", data: { superbubbles: "all" },
      defaultOn: false, urlKey: "hide", urlKeyOn: "show", fadedWithSectors: false,
      title: "Show or hide the faint spheres that mark superbubbles, the cavities the supernovae of young star groups blew in the gas",
    });
    stageView.readLayersFromLocation();
    stageView.layers.refresh();
    bubbleMeshes.forEach(function (mesh) { mesh.visible = stageView.layers.isOn("superbubbles"); });
    stageView.layers.onChange(function (id, on) {
      if (id === "superbubbles") bubbleMeshes.forEach(function (mesh) { mesh.visible = on; });
    });
  }
  // MAP.316: the charted regions are a layer of the Menu too, on unless the visitor hides them.
  if (regionLayer) {
    stageView.layers.register({
      id: "regions", label: "Charted regions", modes: ["galaxy"], group: "Show on the map", data: { regions: "charted" },
      defaultOn: true, urlKey: "hide", urlKeyOn: "show", fadedWithSectors: false,
      title: "Show or hide the soft spheres that mark the regions with generated sectors, colored by how much is charted",
      // MAP.331: while the wireframe cube of the region the view is in is drawn, the legend says what it outlines.
      legend: function () {
        if (!framedRegion || !regionByKey.has(framedRegion)) return null;
        var blob = regionByKey.get(framedRegion);
        return { text: "Outline: " + (blob.name || "a charted region") + " (" + blob.block.m + " sectors a side)", swatchClass: "map-frame-swatch" };
      },
    });
    stageView.readLayersFromLocation();
    stageView.layers.refresh();
    regionsOn = stageView.layers.isOn("regions");
    regionLayer.setVisible(regionsOn && !hiddenCell);
    stageView.layers.onChange(function (id, on) {
      if (id === "regions") {
        regionsOn = on;
        markersDirty = true;
        regionLayer.setVisible(regionsOn && !hiddenCell);
      }
    });
    var regionListRoot = document.getElementById("galaxymap3d-regions");
    if (regionListRoot) regionList = createRegionList(regionListRoot, openRegionBlob);
    canvasEl.__mapRegions = regionLayer;    // the browser tests read what is drawn
    canvasEl.__mapRegionList = regionList;
    // For tests (MAP.316): each region as drawn, where its middle is on the screen (CSS pixels from the canvas's top left),
    // its radius there and how strongly it shows.
    canvasEl.galaxyRegionSpots = function () {
      var width = canvasEl.clientWidth;
      var height = canvasEl.clientHeight;
      var focal = focalPixels();
      var here = camera.position.toArray();
      var point = new THREE.Vector3();
      return regionBlobList.map(function (blob) {
        var look = regionLook(blob, here, focal);
        point.set(blob.x, blob.y, blob.z).project(camera);
        var x = (point.x + 1) / 2 * width;
        var y = (1 - point.y) / 2 * height;
        return { key: blob.key, name: blob.name, share: blob.share, color: blob.color, alpha: look.alpha, radiusPx: look.radiusPx, x: x, y: y,
          visible: point.z > -1 && point.z < 1 && x >= 0 && x <= width && y >= 0 && y <= height };
      });
    };
    // A click on a region's sphere shows its card, a double click flies to it; any star, marker or object under the
    // pointer comes first (the layer is added last).
    picker.addLayer({
      name: "regions",
      priority: 100,
      pick: function (ctx) {
        if (!regionsOn || hiddenCell) return null;
        var focal = focalPixels();
        var here = camera.position.toArray();
        var best = null;
        regionBlobList.forEach(function (blob) {
          var look = regionLook(blob, here, focal);
          if (!(look.alpha > REGION_PICK_MIN_ALPHA)) return;
          regionSpot.set(blob.x, blob.y, blob.z).project(camera);
          if (!(regionSpot.z > -1 && regionSpot.z < 1)) return;
          var x = ctx.rect.left + (regionSpot.x + 1) / 2 * ctx.rect.width;
          var y = ctx.rect.top + (1 - regionSpot.y) / 2 * ctx.rect.height;
          var reach = Math.max(REGION_PICK_MIN_PX, look.radiusPx * REGION_PICK_SHARE);
          var near = Math.hypot(ctx.clientX - x, ctx.clientY - y) / reach;
          if (near <= 1 && (!best || near < best.near)) best = { near: near, blob: blob, distance: look.distance };
        });
        return best ? { entry: best.blob, distance: best.distance } : null;
      },
      tooltip: function (blob) { return (blob.name || "Charted region") + ", " + chartedText(blob.share); },
      hover: function (blob) { previewHover(blob, blob ? function () { return regionFacts(blob); } : null); },
      select: function (blob) {
        highlightPosition(blob.x, blob.y, blob.z);
        showRegionInfo(blob);
      },
    });
    loadRegions();
  }
  regionLookup = regionAt;
  regionOpener = {
    region: function (region) {
      if (region.extent) {
        stageView.goTo({ region: { block: region.block, center: region.extent.centre, radius: region.extent.radius }, label: region.name || null });
      } else {
        stageView.goTo({ block: region.block });
      }
    },
    block: function (block) { stageView.goTo({ block: block }); },
    place: function (target) { stageView.goTo(target); },
  };
  canvasEl.__mapInfo = infoPanel();    // the browser tests put long cards in it (MAP.295)
  // A marker picked (clicked, or its button in the hidden list activated): ring it and show what it is.
  // MAP.305: the card of a crowd's mark: how many systems it stands for, each one a button that shows that system (and flies to it), and one
  // that zooms in on the crowd so it splits.
  function showCrowdInfo(crowd) {
    var list = document.createElement("ul");
    list.className = "map-info-crowd";
    crowd.members.forEach(function (system) {
      var line = document.createElement("li");
      var button = document.createElement("button");
      button.type = "button";
      button.className = "map-info-crowd-system";
      if (crowd.holes) {
        // MAP.340: a crowd of black holes, each a button that shows that hole (and flies to it).
        var mass = system.mass_solar;
        button.textContent = (system.name || "Black hole") + (mass ? ", " + formatMassSolar(mass) : "");
        button.addEventListener("click", function () {
          highlightPosition(system.x, system.y, system.z);
          showPointInfo(system);
          stageView.goTo({ entry: system });
        });
      } else {
        button.textContent = (system.name || "Star system") + ": " + interestLines(system).map(function (reason) { return reason.text; }).join("; ");
        button.addEventListener("click", function () {
          highlightPosition(system.x, system.y, system.z);
          showInterestSystemInfo(system);
          stageView.goTo({ entry: system });
        });
      }
      line.appendChild(button);
      // MAP.349: and a Go to of its own, which flies to it and leaves the card as it is.
      line.appendChild(document.createTextNode(" "));
      line.appendChild(goToButton(system.name || (crowd.holes ? "the black hole" : "the system"), function () { stageView.goTo({ entry: system }); }));
      list.appendChild(line);
    });
    var reach = 0;
    crowd.members.forEach(function (system) {
      reach = Math.max(reach, Math.hypot(system.x - crowd.x, system.y - crowd.y, system.z - crowd.z));
    });
    showInfo({
      title: crowd.count + (crowd.holes ? " black holes" : " points of interest"),
      fields: [["Distance from core", formatDistancePc(Math.hypot(crowd.x, crowd.y, crowd.z))]],
      after: [list],
      buttons: [{ label: "Go to them", onClick: function () { stageView.goTo({ entry: { x: crowd.x, y: crowd.y, z: crowd.z, r: Math.max(reach * 1.5, 0.5) } }); } }],
    });
  }
  function chooseMarker(entry) {
    if (entry.item.merged) {
      highlightPosition(entry.item.x, entry.item.y, entry.item.z);
      showCrowdInfo(entry.item);
      return;
    }
    if (entry.item.regionKey) {
      var region = regionByKey.get(entry.item.regionKey);    // a region's inhabited token picks the region (MAP.321 c)
      if (region) {
        highlightPosition(region.x, region.y, region.z);
        showRegionInfo(region);
      }
      return;
    }
    if (entry.item.bodyRef) {
      if (entry.item.follow && entry.item.follow.kind === "belt") {
        var belt = entry.item.follow.body;     // a belt is not a body to pick: its card says where it is
        showInfo({
          title: entry.item.name || "Asteroid belt",
          kind: "Asteroid belt",
          fields: [["Inner edge", formatDistanceKm(belt.inner_km)], ["Outer edge", formatDistanceKm(belt.outer_km)]],
        });
      } else {
        stageView.selectBodyRef(entry.item.bodyRef);   // a body's own symbol picks the body (MAP.321)
      }
      return;
    }
    highlightPosition(entry.x, entry.y, entry.z);
    if (entry.item.type) showPointInfo(entry.item);
    else showInterestSystemInfo(entry.item);
  }
  // MAP.234: the key. Each thing the map draws has a line in the on-map Key (and in the help dialog's Key section) while it is drawn:
  // a layer's own legend line (a marker kind, a nebula family, the framed region), or a line registered here. The sector scene's kinds
  // have one each here; a kind with no line is a layer the key says nothing about, which `registry.unkeyed()` lists and the browser test
  // holds empty, so a new layer cannot be added without one.
  var KIND_KEYS = {
    star: null,    // the stars' line below stands for it
    nebula: ["cloud", "Clouds: nebulae"],
    supernovaRemnant: ["cloud", "Shells: supernova remnants"],
    asteroidField: ["cloud", "Clouds: asteroid fields"],
    blackHole: ["point", "Bright points and dark spheres: black holes"],
    neutronStar: ["point", "Bright points: neutron stars"],
    quasar: ["point", "Bright points: quasars"],
    roguePlanet: ["point", "Faint points: rogue planets"],
    interstellarComet: ["point", "Small spheres: interstellar comets"],
    neighbor: ["outline", "Outlines: the neighboring sectors"],
  };
  (function registerKeyLines() {
    var layers = stageView.layers;
    layers.layers().filter(function (layer) { return layer.data && layer.data.kind; }).forEach(function (layer) {
      var entry = KIND_KEYS[layer.id];
      if (!entry) return;    // the stars', or left for `unkeyed()` to name
      layers.registerKey({
        id: "key-" + layer.id, layer: layer.id, text: entry[1], swatchClass: "map-key-" + entry[0],
        modes: ["galaxy", "sector"], drawn: function () { return layers.isOn(layer.id) && (!layer.available || layer.available()); },
      });
    });
    var starLayer = layers.get("star");
    layers.registerKey({
      id: "key-stars", layer: "star", text: "Points of light: stars, colored by temperature, brighter the more luminous", swatchClass: "map-key-point",
      modes: ["galaxy", "sector"],
      drawn: function () { return data.pinned ? !!starLayer && layers.isOn("star") && (!starLayer.available || starLayer.available()) : starList.length > 0; },
    });
    // The legend follows what is drawn as the view changes.
    window.setInterval(function () { layers.refreshKey(); }, 400);
    if (data.pinned) return;
    var lines = function () { return canvasEl.galaxyLines ? canvasEl.galaxyLines() : null; };
    layers.registerKey({
      id: "key-blocks", text: "Blocks: predicted density, brighter is denser; more solid the more of their sectors are generated",
      swatchClass: "map-key-blocks", drawn: function () { var l = lines(); return !!l && l.blocksDrawn > 0; },
    });
    layers.registerKey({
      id: "key-charted", text: "Outlines: the blocks that hold generated sectors", swatchClass: "map-key-outline",
      drawn: function () { var l = lines(); return !!l && l.chartedLines > 0; },
    });
    layers.registerKey({
      id: "key-glow", layer: "glow", text: "Haze: stars too faint and close together to tell apart", swatchClass: "map-key-cloud",
      drawn: function () { return layers.isOn("glow") && starList.length > 0; },
    });
    layers.registerKey({
      id: "key-bubbles", text: "Faint blue spheres with a rim: superbubbles, cavities the supernovae of young star groups blew in the gas",
      swatchClass: "map-key-sphere", drawn: function () { return bubbles.length > 0 && stageView.layers.isOn("superbubbles"); },
    });
    if (regionLayer) {
      layers.registerKey({
        id: "key-regions", layer: "regions", text: "Soft spheres: regions with generated sectors, colored by how much is charted",
        swatchClass: "map-key-cloud", drawn: function () { return regionsOn && !hiddenCell && regionLayer.drawn() > 0; },
      });
    }
    if (course) {
      layers.registerKey({ id: "key-course", text: "Line: the course, through every stop", swatchClass: "map-key-course", drawn: function () { return true; } });
      if (course.direct && course.direct.length === 2) {
        layers.registerKey({ id: "key-direct", text: "Dashed line: the direct line from the start to the end", swatchClass: "map-key-direct", drawn: function () { return true; } });
      }
    }
  })();
  // MAP.268: the names (galaxylabels.js, through the label engine maplabels.js): sectors, the open sector's systems and
  // nebulae, a layer each in the Menu's Labels group.
  var galaxyLabels = null;
  if (viewport) {
    LABEL_KINDS.forEach(function (kind) {
      stageView.layers.register({
        id: kind.id, label: kind.label, modes: ["galaxy", "sector"], group: "Labels", data: { label: kind.kind },
        defaultOn: true, urlKey: "hide", urlKeyOn: "show", fadedWithSectors: false, title: kind.title,
        keyless: "names laid over the map, not a thing it draws",
      });
    });
    stageView.readLayersFromLocation();
    stageView.layers.refresh();
    galaxyLabels = createGalaxyLabels({
      THREE: THREE, camera: camera, canvasEl: canvasEl, parent: viewport, layers: stageView.layers, light: isLightBackground(),
      // MAP.299: the clear margin between two names, a preference read at each placement.
      padding: function () {
        try { return labelBufferPx(window.localStorage); } catch (error) { return labelBufferPx(null); }
      },
      sources: {
        sectors: stageView.sectorLabels, systems: stageView.systemLabels, bodies: bodyNames, nebulae: cloudLabels,
        stars: starNames, objects: objectNames, worlds: worldNames, regions: regionNames,
      },
    });
    canvasEl.__mapLabels = galaxyLabels;    // the browser tests read what is placed
    canvasEl.__mapSectorNames = stageView.sectorLabels;    // ... and the sector names with the total mass each is ranked by (MAP.322)
    canvasEl.__mapBodyNames = bodyNames;    // ... and the slice each body's name is offered in (MAP.321)
  }
  // MAP.314: the names of the charted regions the nodes of the view offer, by the size of the nodes (regionlabels.js): none
  // in a view with no charted region, and none while the nodes are not in. Only the nodes the view is drawn with are asked,
  // not a block that has handed over to its children, so a family is named once, at the size the view shows it.
  function regionNames() {
    var need = neededNodes();
    var settled = settledExpanded(need, function (key) { return getNode(key) !== undefined; });
    var drawn = need.keys.filter(function (key) { return !settled.has(key); });
    return regionLabels(drawn.map(function (key) { return { key: key, node: getNode(key) }; })).filter(function (label) {
      return drawnByTheGalaxy(label.x, label.y, label.z);
    });
  }
  // The nebulae and remnants in view that are large enough on screen to read: their name, else their type.
  function cloudLabels() {
    var found = [];
    cloudSprites.forEach(function (sprite) {
      var cloud = sprite.userData.cloud;
      var distance = camera.position.distanceTo(sprite.position);
      if (!drawnByTheGalaxy(sprite.position.x, sprite.position.y, sprite.position.z) || distance < 1.2 * cloud.radius_pc) {
        return;
      }
      found.push({
        id: "cloud:" + cloudKey(cloud), text: cloud.name || cloudTypeLabel(cloud),
        x: sprite.position.x, y: sprite.position.y, z: sprite.position.z, pixels: 2 * cloudPixels(cloud, distance),
        ref: { kind: "nebula", cloud: cloud },
        slice: { domain: "nebula", volume: 4 / 3 * Math.PI * Math.pow(cloud.radius_pc || 0, 3) },
      });
    });
    return found;
  }
  // MAP.327: the wireframe round the charted region the view is in; worked out when the point the camera looks at, the regions or
  // the opened sector change, not every frame.
  function updateRegionFrame() {
    if (!regionFrame) return;
    var at = [target.x, target.y, target.z, regionsVersion, regionsOn, !!hiddenCell].join(",");
    if (at === framedAt) return;
    framedAt = at;
    var key = "";
    if (regionsOn && !hiddenCell && regionBlobList.length) {
      var cell = cellAt(9, target.x, target.y, target.z, edgePc);
      var cellKey = cell.m + "." + cell.shell + "." + cell.band + "." + cell.wedge;
      if (regionByKey.has(cellKey)) key = cellKey;
    }
    if (key === framedRegion) return;
    framedRegion = key;
    regionFrame.set(key ? drillBlockBounds(regionByKey.get(key).block, edgePc) : null);
    stageView.layers.refresh();    // the legend names the outline while it is drawn (MAP.331)
  }
  var SECTOR_MARKER_BUDGET = 1500;
  var SECTOR_MERGE_PX = 12;
  var markerFrustum = new THREE.Frustum();
  var markerMatrix = new THREE.Matrix4();
  var markerPoint = new THREE.Vector3();
  // Chooses what the marker layer draws from the nodes of this view (galaxymarkers.selectMarkers): the kinds that are
  // on, the items in front of the camera and outside any clipped wedge, no more than the budget.
  function refreshMarkers() {
    if (!markerLayer) return;
    markersDirty = false;
    camera.updateMatrixWorld();
    markerMatrix.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);
    markerFrustum.setFromProjectionMatrix(markerMatrix);
    var nodes = neededNodes().keys.map(function (key) { return { key: key, node: getNode(key) }; });
    // MAP.339: the black holes are the ones near what the camera looks at, as many as the density preference lets the screen hold.
    var densityStorage = null;
    try { densityStorage = window.localStorage; } catch (error) { densityStorage = null; }
    // In an opened sector there is room to show every system that fits: the budget is raised and the marks are merged only when they touch.
    var sectorOpen = !!(canvasEl.galaxyOpenSector && canvasEl.galaxyOpenSector());
    var chosen = selectMarkers({
      budget: sectorOpen ? SECTOR_MARKER_BUDGET : undefined,
      near: nearFieldRule({
        centre: [target.x, target.y, target.z], distance: camera.position.distanceTo(target),
        width: canvasEl.clientWidth || 1, height: canvasEl.clientHeight || 1, densityPx: blackHoleDensityPx(densityStorage),
      }),
      kinds: markerKinds.kinds(),
      isOn: function (id) { return stageView.layers.isOn(id); },
      isPulled: function (id) { return stageView.layers.isPulled(id); },
      nodes: nodes,
      inView: function (item) {
        return markerFrustum.containsPoint(markerPoint.set(item.x, item.y, item.z)) && drawnByTheGalaxy(item.x, item.y, item.z);
      },
    });
    // MAP.321: inside the open system the star does not repeat its symbol, its bodies carry their own.
    var inside = bodyHandoffSource();
    if (inside) {
      chosen = chosen.filter(function (entry) { return entry.item.id !== inside.id; })
        .concat(bodyMarkers(INTEREST_KINDS, function (id) { return stageView.layers.isOn(id); }, inside.entries));
    }
    // MAP.321 c: out at galactic scale a region with an inhabited world carries a token, while its sphere shows and no system's own mark
    // stands for it.
    var inhabitedKind = INHABITED_MARKER_KIND;
    if (regionLayer && regionsOn && !hiddenCell && inhabitedKind && stageView.layers.isOn(inhabitedKind.id)) {
      var here = camera.position.toArray();
      var focal = focalPixels();
      var marked = chosen.filter(function (entry) { return entry.kind.interest === "inhabited" && !entry.item.bodyRef; });
      chosen = chosen.concat(regionMarkers(inhabitedKind, regionBlobList,
        function (blob) { return regionLook(blob, here, focal); },
        function (blob) { return markerFrustum.containsPoint(markerPoint.set(blob.x, blob.y, blob.z)); },
        function (blob) {
          return marked.some(function (entry) {
            return Math.hypot(entry.item.x - blob.x, entry.item.y - blob.y, entry.item.z - blob.z) <= blob.radiusPc;
          });
        }));
    }
    // MAP.306: a symbol of a reason pulled forward is drawn over the others and named first, whatever level it is drawn at.
    chosen.forEach(function (entry) {
      if (!entry.pulled && entry.kind.interest && stageView.layers.isPulled(entry.kind.id)) entry.pulled = true;
    });
    // MAP.305: marks too close on the screen to tell apart are one mark with a count.
    var viewW = canvasEl.clientWidth || 1;
    var viewH = canvasEl.clientHeight || 1;
    chosen = mergeMarks(chosen, function (item) {
      markerPoint.set(item.x, item.y, item.z).project(camera);
      return { x: (markerPoint.x + 1) / 2 * viewW, y: (1 - markerPoint.y) / 2 * viewH,
        depth: Math.hypot(item.x - camera.position.x, item.y - camera.position.y, item.z - camera.position.z) };
    }, sectorOpen ? SECTOR_MERGE_PX : undefined);
    markerLayer.set(chosen);
  }
  // MAP.321: the open system while the camera is close enough for the bodies to carry the symbols (bodymarkers.createHandoff's
  // gap keeps a camera resting at the boundary from flickering between the levels), or null.
  var bodyHandoff = createHandoff();
  var lastBodyEntries = null;
  function bodyHandoffSource() {
    var source = stageView.systemMarkerSource();
    if (!source || !source.centre) {
      bodyHandoff.update(0, 0, false);
      return null;
    }
    var distance = Math.hypot(camera.position.x - source.centre[0], camera.position.y - source.centre[1], camera.position.z - source.centre[2]);
    return bodyHandoff.update(distance, source.radiusPc, true) ? source : null;
  }
  canvasEl.__bodyMarkersOn = function () { return bodyHandoff.isOn(); };
  // MAP.200: the on-map zoom slider and direction star.
  if (viewport) navPad = createNavPad(viewport, stageView, document.getElementById("galaxymap3d-centre"));

  // --- Pointer and keys ------------------------------------------------------
  //
  // Everything goes to the drill-down, which turns and zooms the view
  // where it may.
  // A click on a black hole, neutron star or quasar, or a cloud small
  // enough to aim at, shows it instead of picking what's under it (never a
  // star: MAP.101).
  canvasEl.addEventListener("pointerdown", function (event) { stageView.onPointerDown(event); });
  canvasEl.addEventListener("pointermove", function (event) { stageView.onPointerMove(event); });
  canvasEl.addEventListener("pointerup", function (event) { stageView.onPointerUp(event); });
  canvasEl.addEventListener("pointercancel", function (event) { stageView.onPointerUp(event); });
  canvasEl.addEventListener("pointerleave", function (event) { stageView.onPointerLeave(event); });
  canvasEl.addEventListener("keydown", function (event) { stageView.onKey(event); });
  // The wheel zooms within the stage's zoom policy.
  canvasEl.addEventListener("wheel", function (event) {
    if (stageView.onWheel(event)) {
      event.preventDefault();
    }
  }, { passive: false });

  // --- Picking (mappick.js) --------------------------------------------------
  //
  // Inside a container a black hole, neutron star or quasar, or a cloud
  // small enough to aim at, takes a click before the block under it (a
  // point inside a cloud, before the cloud)
  // (the drill-down's own layer, added last); stars never do (MAP.101):
  // in a dense sector they would cover it, so a click on one picks the
  // block, slab or sector under it, and a star's details are on its
  // sector's page.
  function phenomenonHover(entry) {
    if (entry) hoverRing.at(entry.x, entry.y, entry.z, { px: HIGHLIGHT_PX * 2 });
    else hoverRing.hide();
  }

  // MAP.348: a thing hovered shows a small preview (its card's first view) in a corner of the map that is clear of the pointer, of the thing
  // and of the open card. It takes no pointer event, and a touch screen has no hover to show it.
  var hoverPreview = createHoverPreview(canvasEl.closest(".galaxy-map-main") || viewport);
  var lastPointer = null;
  canvasEl.addEventListener("pointermove", function (event) { lastPointer = { x: event.clientX, y: event.clientY }; });
  ["pointerleave", "pointerdown", "wheel", "keydown"].forEach(function (type) { canvasEl.addEventListener(type, function () { hoverPreview.hide(); }); });
  canvasEl.__mapHoverPreview = hoverPreview;    // the browser tests read what it shows
  var previewSpot = new THREE.Vector3();
  function previewPoints(entry) {
    var points = lastPointer ? [lastPointer] : [];
    previewSpot.set(entry.x, entry.y, entry.z).project(camera);
    if (previewSpot.z > -1 && previewSpot.z < 1) {
      var rect = canvasEl.getBoundingClientRect();
      points.push({ x: rect.left + (previewSpot.x + 1) / 2 * rect.width, y: rect.top + (1 - previewSpot.y) / 2 * rect.height });
    }
    return points;
  }
  // `thing` is what is under the pointer (null when nothing), `facts()` its card's spec; `at` is the entry with the thing's x, y, z.
  function previewHover(thing, facts, at) {
    if (!thing) hoverPreview.hide();
    else hoverPreview.show(thing, facts, function () { return previewPoints(at || thing); });
  }
  // The hover of a layer whose things have a card: the ring (phenomenonHover) and the preview.
  function hoverWith(facts) {
    return function (entry) {
      phenomenonHover(entry);
      previewHover(entry, entry ? function () { return facts(entry); } : null);
    };
  }
  // A marker is a star system, a black hole, neutron star or quasar, or a region's token; a crowd or a body has no first view of its own.
  function markerFacts(marker) {
    var item = marker.item;
    if (item.merged || item.bodyRef) return null;
    if (item.regionKey) {
      var region = regionByKey.get(item.regionKey);
      return region ? regionFacts(region) : null;
    }
    return item.type ? pointFacts(item) : interestFacts(item);
  }

  picker.addLayer({
    name: "clouds",
    enabled: function () { return stageView.inContainer(); },
    pick: function (ctx) {
      var found = cloudAtRay(ctx.ray);
      return found && found.core ? { entry: found.cloud, distance: found.distance } : null;
    },
    tooltip: function (cloud) { return (cloud.name ? cloud.name + ", " : "") + cloudTypeLabel(cloud); },
    hover: hoverWith(cloudFacts),
    select: function (cloud) {
      highlightPosition(cloud.x, cloud.y, cloud.z);
      showCloudInfo(cloud);
    },
  });

  picker.addLayer({
    name: "phenomena",
    enabled: function () { return stageView.inContainer(); },
    points: function () { return starList; },
    reach: function () { return POINT_PICK_PX; },
    accept: function (star) { return star.phenomenon && drawnByTheGalaxy(star.x, star.y, star.z); },
    lastWins: true,
    tooltip: function (point) { return (point.name ? point.name + ", " : "") + pointTypeLabel(point); },
    hover: hoverWith(pointFacts),
    select: function (point) {
      highlightPosition(point.x, point.y, point.z);
      showPointInfo(point);
    },
  });

  // MAP.304: a star of the map takes a click (selected, its panel shown) and a double click (the camera flies to it) from
  // any distance: the nearest drawn star under the pointer, the brightest of those level with it. A star is pickable
  // while it is drawn, that is, its node's share and its apparent brightness both leave it visible, as the shader works
  // them out (starmagnitude.js's twin). The list is worked out once per view and star set, not per pointer event.
  // MAP.122: the pointer's reach on a star, 10 px for a mouse and 22 px for a finger (a thumb covers more than a cursor),
  // and when several stars lie inside it the panel offers the nearest few to choose from.
  var STAR_MOUSE_PICK_PX = 10;
  var STAR_TOUCH_PICK_PX = 22;
  var STAR_CHOICES = 6;
  var lastPointer = { type: "mouse", x: 0, y: 0 };
  ["pointerdown", "pointerup", "pointermove"].forEach(function (name) {
    canvasEl.addEventListener(name, function (event) {
      lastPointer = { type: event.pointerType === "touch" ? "touch" : "mouse", x: event.clientX, y: event.clientY };
    }, true);
  });
  function starReach() {
    return lastPointer.type === "touch" ? STAR_TOUCH_PICK_PX : STAR_MOUSE_PICK_PX;
  }
  // The stars within the reach of the last pointer, nearest first (the brightest of those level with it first).
  function starsNearPointer() {
    return entriesNearScreenPoint(pickableStars(), camera, canvasEl.getBoundingClientRect(), lastPointer.x, lastPointer.y,
      { reach: starReach, lastWins: true, limit: STAR_CHOICES }).map(function (found) { return found.entry; });
  }
  // The choice of a crowded click: a button for each star within reach, the shown one marked.
  function starChoiceList(stars, chosen) {
    var list = document.createElement("ul");
    list.className = "map-star-choices";
    list.setAttribute("aria-label", "Stars within reach of the click");
    stars.forEach(function (star) {
      var line = document.createElement("li");
      var button = document.createElement("button");
      button.type = "button";
      button.className = "starmap-btn";
      button.dataset.starChoice = String(star.id);
      button.textContent = (star.name || "Star") + ", " + formatLuminositySol(star.luminosity_sol);
      button.setAttribute("aria-pressed", String(star === chosen));
      button.addEventListener("click", function () {
        highlightPosition(star.x, star.y, star.z);
        showStarInfo(star, starChoiceList(stars, star));
      });
      line.appendChild(button);
      list.appendChild(line);
    });
    return list;
  }
  var STAR_PICK_MIN_SHARE = 0.25;
  var pickableCache = { signature: "", list: [] };

  function pickableStars() {
    var n = starList.length;
    var p = camera.position;
    var uniforms = starMaterial.uniforms;
    var container = stageView.inContainer();
    var signature = [starSetVersion, n, container, p.x, p.y, p.z, uniforms.magLimit.value, uniforms.magOn.value].join(",");
    if (signature === pickableCache.signature) return pickableCache.list;
    updateStarShares();
    var absolute = starPoints.geometry.getAttribute("starAbsMag");
    var tierMag = starPoints.geometry.getAttribute("starTierMag");
    var clipped = starPoints.geometry.getAttribute("starClipped");
    var list = [];
    var objects = [];
    if (absolute && tierMag) {
      for (var i = 0; i < n; i++) {
        var star = starList[i];
        if (!(starShares[i] >= STAR_PICK_MIN_SHARE) || (clipped && clipped.array[i] > 0.5)) continue;
        // The other objects that are named (MAP.298): kept for the labels, picked by their own layer.
        if (star.phenomenon) {
          if (star.name && drawnByTheGalaxy(star.x, star.y, star.z)) objects.push(star);
          continue;
        }
        if (!container && stageView.outsideOutline(star)) continue;  // the halo layer below names those
        if (uniforms.magOn.value > 0.5) {
          var dx = star.x - p.x;
          var dy = star.y - p.y;
          var dz = star.z - p.z;
          var appMag = apparentMagnitude(absolute.array[i], Math.sqrt(dx * dx + dy * dy + dz * dz)) + tierMag.array[i];
          if (magnitudeShare(appMag, uniforms.magLimit.value) < STAR_PICK_MIN_SHARE) continue;
        }
        if (!drawnByTheGalaxy(star.x, star.y, star.z)) continue;
        list.push(star);
      }
    }
    // The brightest last, so among stars level on the screen the brightest wins (the layer's lastWins).
    list.sort(function (a, b) { return a.luminosity_sol - b.luminosity_sol; });
    pickableCache = { signature: signature, list: list, objects: objects };
    return list;
  }

  // MAP.298: what the label engine may name besides the open sector's systems: the named stars the map draws and the
  // named objects of the other kinds, the brightest few hundred of each (the engine keeps the ones with room).
  var NAMED_LIMIT = 400;
  function starNames() {
    var found = [];
    if (hiddenCell) return found;
    var list = pickableStars();
    for (var i = list.length - 1; i >= 0 && found.length < NAMED_LIMIT; i--) {
      var star = list[i];
      if (star.name) {
        found.push({ id: "star:" + star.id, text: star.name, x: star.x, y: star.y, z: star.z,
          luminositySol: star.luminosity_sol || 0, ref: { kind: "star", star: star },
          slice: { domain: "star", radiusSol: star.radius_sol || 0, luminositySol: star.luminosity_sol || 0 } });
      }
    }
    return found;
  }
  function objectNames() {
    if (hiddenCell) return [];
    pickableStars();
    return pickableCache.objects.slice(0, NAMED_LIMIT).map(function (point) {
      return { id: "object:" + point.type + ":" + point.id, text: point.name, x: point.x, y: point.y, z: point.z,
        luminositySol: point.map_luminosity_sol || point.luminosity_sol || 0, ref: { kind: "object", point: point },
        slice: { domain: "uninhabited", category: point.type === "rogue_planet" ? "rogue planet" : "non-luminous star-like object", mass: point.mass_solar || 0 } };
    });
  }

  // The world markers the map draws, offered as names of the inhabited slice (MAP.317, MAP.321): a name is only offered for a
  // marker that is drawn, and an inhabited world (a star's inhabited token) is always a candidate, named before a habitable one.
  // The node lists a system's worlds as a count and the best score, which rank them until the stored population does (MAP.318).
  var NAMED_WORLD_REASONS = { inhabited: "inhabited world", habitable: "habitable world" };
  function worldNames() {
    if (!markerLayer) return [];
    return markerLayer.entries().filter(function (entry) {
      return (NAMED_WORLD_REASONS[entry.kind.interest] || entry.pulled) && entry.item.name && !entry.item.bodyRef && !entry.item.regionKey;
    }).map(function (entry) {
      var world = NAMED_WORLD_REASONS[entry.kind.interest];
      // A reason pulled forward (MAP.306) is named first and kept when names are cut; one with no inhabitants of its own is in the uninhabited slice.
      var slice = world
        ? { domain: "inhabited", category: world, inhabited: entry.kind.interest === "inhabited", habitability: entry.item.score || 0, population: entry.item.population || 0 }
        : { domain: "uninhabited", category: entry.kind.bodyName || entry.kind.interest, size: 0 };
      if (entry.pulled) slice.first = true;
      return { id: "world:" + entry.item.id, text: entry.item.name, x: entry.item.x, y: entry.item.y, z: entry.item.z, ref: { kind: "marker", entry: entry }, slice: slice };
    });
  }
  // The bodies of the open system, named as the stage offers them, except that a planet or moon whose inhabited or habitable symbol is
  // drawn takes its place in the inhabited slice with the worlds (MAP.321), by its habitability score, and one whose reason is pulled forward
  // is named first (MAP.306).
  function bodyNames() {
    var labels = stageView.bodyLabels();
    if (!markerLayer) return labels;
    var worlds = new Map();
    markerLayer.entries().forEach(function (entry) {
      if (entry.item.bodyRef && (NAMED_WORLD_REASONS[entry.kind.interest] || entry.pulled)) worlds.set(entry.item.bodyRef, entry);
    });
    if (!worlds.size) return labels;
    return labels.map(function (label) {
      var entry = worlds.get(label.id);
      if (!entry) return label;
      var body = label.ref.entry.body || {};
      var slice = NAMED_WORLD_REASONS[entry.kind.interest]
        ? { domain: "inhabited", category: label.ref.entry.kind, inhabited: entry.kind.interest === "inhabited", habitability: body.phi4 || 0, population: 0 }
        : Object.assign({}, label.slice);
      if (entry.pulled) slice.first = true;
      return Object.assign({}, label, { slice: slice });
    });
  }

  // A click on a name (MAP.298): the item is selected as a click on it would and the camera flies to it.
  function actOnLabel(label) {
    var ref = label.ref;
    if (!ref) return;
    if (ref.kind === "star") {
      highlightPosition(ref.star.x, ref.star.y, ref.star.z);
      showStarInfo(ref.star);
      stageView.goTo({ entry: ref.star });
    } else if (ref.kind === "object") {
      highlightPosition(ref.point.x, ref.point.y, ref.point.z);
      showPointInfo(ref.point);
      stageView.goTo({ entry: ref.point });
    } else if (ref.kind === "marker") {
      chooseMarker(ref.entry);
      stageView.goTo({ entry: ref.entry });
    } else if (ref.kind === "system") {
      stageView.selectEntry(ref.entry);
      stageView.goTo({ entry: ref.entry });
    } else if (ref.kind === "body") {
      stageView.selectBody(ref.entry);
    } else if (ref.kind === "nebula") {
      highlightPosition(ref.cloud.x, ref.cloud.y, ref.cloud.z);
      showCloudInfo(ref.cloud);
      stageView.goTo({ entry: { x: ref.cloud.x, y: ref.cloud.y, z: ref.cloud.z, r: ref.cloud.radius_pc } });
    } else if (ref.kind === "sector") {
      stageView.goTo({ sector: ref.sector, label: label.text || null });
    } else if (ref.kind === "region") {
      // MAP.314: a charted region is framed by its sphere, a family or a region listed without one by its block.
      var block = parseNodeKey(ref.key);
      if (ref.extent) regionOpener.region({ block: block, extent: ref.extent });
      else regionOpener.block(block);
    }
  }

  // For tests (MAP.314): the keys of the nodes the view is cut into, and whether each is in yet.
  canvasEl.galaxyNodeKeys = function () {
    return neededNodes().keys.map(function (key) { return { key: key, loaded: !!getNode(key) }; });
  };

  // For tests (MAP.304): the stars a click can pick now, brightest first, with where each is on the screen (CSS pixels
  // from the canvas's top left; `visible` false off the screen).
  // For tests (MAP.240): every layer built relative to an origin of its own -- its name, how far its farthest buffer value is
  // from the object's origin in parsecs (the value times the object's scale), and where the origin is.
  canvasEl.galaxyBufferExtents = function () {
    var found = [];
    scene.updateMatrixWorld(true);
    scene.traverse(function (object) {
      var layer = object.userData && object.userData.relativeLayer;
      var position = layer && object.geometry && object.geometry.getAttribute("position");
      if (!position) return;
      found.push({
        layer: layer, extentPc: extentOf(position.array) * object.matrixWorld.getMaxScaleOnAxis(),
        origin: new THREE.Vector3().setFromMatrixPosition(object.matrixWorld).toArray(), count: position.count,
      });
    });
    return found;
  };

  canvasEl.galaxyPickableStars = function () {
    var width = canvasEl.clientWidth;
    var height = canvasEl.clientHeight;
    var point = new THREE.Vector3();
    return pickableStars().slice().reverse().map(function (star) {
      point.set(star.x, star.y, star.z).project(camera);
      var x = (point.x + 1) / 2 * width;
      var y = (1 - point.y) / 2 * height;
      return { id: star.id, system_id: star.system_id, name: star.name || null, luminosity: star.luminosity_sol, x: x, y: y,
        visible: point.z > -1 && point.z < 1 && x >= 0 && x <= width && y >= 0 && y <= height };
    });
  };

  picker.addLayer({
    name: "stars",
    points: pickableStars,
    reach: starReach,
    lastWins: true,
    tooltip: function (star) { return (star.name ? star.name + ", " : "") + "star" + (star.star_type && star.star_type !== "fixture" ? " " + star.star_type : ""); },
    hover: hoverWith(starFacts),
    select: function (star) {
      highlightPosition(star.x, star.y, star.z);
      var crowd = starsNearPointer();
      showStarInfo(star, crowd.length > 1 && crowd.indexOf(star) >= 0 ? starChoiceList(crowd, star) : null);
    },
  });

  // MAP.243: on the whole-galaxy view the halo and the satellite galaxies lie beyond the outline, where no block
  // is, so a click would find nothing. A star or a phenomenon there is picked on the screen like a phenomenon
  // inside a window (a double click flies to it; the wheel then opens the window round its cell).
  picker.addLayer({
    name: "halo",
    enabled: function () { return !stageView.inContainer(); },
    points: function () { return starList; },
    reach: function () { return HALO_PICK_PX; },
    accept: function (star) { return stageView.outsideOutline(star); },
    lastWins: true,
    tooltip: function (point) {
      return point.phenomenon ? (point.name ? point.name + ", " : "") + pointTypeLabel(point) : "Star beyond the galaxy's outline";
    },
    hover: hoverWith(function (point) { return point.phenomenon ? pointFacts(point) : haloStarFacts(point); }),
    select: function (point) {
      highlightPosition(point.x, point.y, point.z);
      if (point.phenomenon) showPointInfo(point);
      else showHaloStarInfo(point);
    },
  });

  // MAP.219: a click or hover on a marker picks it before the point of light under it (a black hole shows its own page
  // card, as the point does).
  if (markerLayer) {
    picker.addLayer({
      name: "markers",
      priority: -1,
      enabled: function () { return markerLayer.entries().length > 0; },
      points: function () { return markerLayer.entries(); },
      reach: function () { return MARKER_PX / 2 + 3; },
      tooltip: function (entry) { return entry.kind.describe(entry.item); },
      hover: function (entry) {
        phenomenonHover(entry);
        var facts = entry && markerFacts(entry);
        previewHover(facts ? entry.item : null, function () { return facts; }, entry);
      },
      select: function (entry) { chooseMarker(entry); },
    });
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
      setChartedOnly: setChartedOnly,
    }));
  }
  var menuEl = document.getElementById("galaxymap3d-menu");
  if (menuEl) {
    wireMapMenu(menuEl);
  }
  var colorSelect = document.getElementById("galaxymap3d-color-by");
  if (colorSelect) {
    colorSelect.addEventListener("change", function () { stageView.setColorBy(colorSelect.value); });
    colorSelect.value = stageView.colorBy();
  }

  // --- Scale bar -----------------------------------------------------------

  var scaleEl = document.getElementById("galaxymap3d-scale");

  // Up to 3 significant figures, grouped: 0.0512, 3.4, 1,280.
  function formatCount(value) {
    if (!(value > 0)) return "0";
    if (value >= 100) return formatNumber(value);
    return threeFigures(value);
  }

  function plural(count, word) {
    return count === 1 ? word : word + "s";
  }

  // The readout's label: "<sectors> · <distance>", the distance on the
  // shared ladder (static/unitladder.js), which adds ly in parentheses to
  // parsec values.
  function scaleLabel(pc) {
    var sectors = pc / edgePc;
    var label = document.createElement("span");
    label.className = "starmap-scale-label";
    label.textContent = formatCount(sectors) + " " + plural(sectors, "sector") + " · " + formatDistancePc(pc);
    return label;
  }

  // One scale line (MAP.60): a bar of about SCALE_BAR_PX rounded to a nice
  // number of parsecs, and its length. Per CSS pixel: it is about what a
  // person can see, not device pixels.
  var SCALE_BAR_PX = 70;

  function updateScaleBar() {
    if (!scaleEl) {
      return;
    }
    var pcPerScreenPx = pcPerPixelAtTarget();
    var nicePc = niceScaleValue(SCALE_BAR_PX * pcPerScreenPx);
    if (!nicePc) {
      return;
    }
    var bar = document.createElement("span");
    bar.className = "starmap-scale-bar";
    bar.style.width = (nicePc / pcPerScreenPx).toFixed(1) + "px";
    scaleEl.replaceChildren(bar, scaleLabel(nicePc));
  }

  // --- Resize/render loop ----------------------------------------------

  // MAP.291: the window is as wide as its column and as tall as the browser window leaves room for.
  // MAP.350: the Galaxy Map is as tall as the browser window less 20 pixels above and below, and the page rests with it there.
  // MAP.359: a sector page's map (pinned) is the same size, so the sector's map matches the Galaxy Map for the window.
  watchWindowFit(viewport, null, { fill: true });

  watchResize(viewport, function () {
    fitRendererToCanvas(renderer, camera, canvasEl);
    updateScaleBar();
  });

  // While a NAV course is picked only a choice holding something generated
  // can be taken (an end lives in a generated sector); every control,
  // "Charted only" included, works as when browsing (NAV.32).
  function takeOnlyGenerated(on) {
    stageView.setNeedGenerated(on);
  }

  // The pick moved on or ended (navpick.js): the banner, the Bookmarks menu,
  // the URL and the panel follow, and the view stays where it is.
  function pickChanged() {
    takeOnlyGenerated(navPick.active());
    renderPickBanner();
    keepPickInBookmarks();
    stageView.pickChanged();
    var panel = infoPanel();
    if (panel) {
      panel.hint(navPick.active()
        ? navPick.bannerParts("object")[0] + ": zoom out or in to find it, then click it."
        : "NAV pick cancelled.", false);
    }
  }

  onPickChanged = pickChanged;

  // The banner under the map's title: the page's own when the NAV page began
  // the pick, else made here.
  function renderPickBanner() {
    var banner = document.querySelector(".pick-banner");
    var parts = navPick.bannerParts("object");
    if (!parts) {
      if (banner) banner.remove();
      return;
    }
    if (!banner) {
      banner = document.createElement("p");
      banner.className = "pick-banner";
      banner.setAttribute("role", "status");
      var layout = document.querySelector(".galaxymap3d-panel .starmap-layout");
      if (!layout || !layout.parentNode) return;
      layout.parentNode.insertBefore(banner, layout);
    }
    banner.textContent = "";
    var lead = document.createElement("strong");
    lead.textContent = parts[0] + ":";
    banner.appendChild(lead);
    banner.appendChild(document.createTextNode(parts[1]));
    var cancelUrl = navPick.cancelUrl();
    var cancel = document.createElement(cancelUrl ? "a" : "button");
    cancel.textContent = "Cancel";
    if (cancelUrl) {
      cancel.href = cancelUrl;
    } else {
      cancel.type = "button";
      cancel.className = "btn btn-small";
      cancel.addEventListener("click", function () { navPick.clear(); });
    }
    banner.appendChild(cancel);
  }

  // The Bookmarks menus carry the pick (NAV.40, bookmarks.js reads these).
  function keepPickInBookmarks() {
    Array.prototype.forEach.call(document.querySelectorAll("[data-bookmarks-menu]"), function (menu) {
      if (navPick.active()) {
        menu.setAttribute("data-pick", navPick.pick());
        menu.setAttribute("data-keep-name", navPick.keepName());
        menu.setAttribute("data-keep-value", navPick.other() || "");
        menu.setAttribute("data-nav-url", String(sceneData.navUrl || "/nav"));
      } else {
        ["data-pick", "data-keep-name", "data-keep-value", "data-nav-url"].forEach(function (name) {
          menu.removeAttribute(name);
        });
      }
    });
  }

  // NAV's "Pick on Galaxy Map": endpoints live only in generated sectors.
  if (navPick.active()) {
    takeOnlyGenerated(true);
  }

  // The drill-down, at the stage the URL names.
  stageView.openFromLocation();
  stageView.setActive(true);

  // Every frame: the drill-down's camera move, then (when the view moved)
  // the stars and clouds from cached nodes and a fetch for missing ones.
  var lastViewSignature = "";
  var lastRenderAt = 0;
  // A camera that has not left the cut it has still moves the nodes' shares, and so which nodes are drawn, slowly.
  var RENDER_EVERY_MS = 100;

  (function animate(now) {
    requestAnimationFrame(animate);
    now = now || performance.now();
    stageView.step(now);
    var viewSignature = [camera.position.x, camera.position.y, camera.position.z, target.x, target.y, target.z,
      canvasEl.clientHeight].join(",");
    if (viewSignature !== lastViewSignature) {
      var need = neededNodes();
      if (need !== lastNeed || now - lastRenderAt >= RENDER_EVERY_MS) {
        lastViewSignature = viewSignature;
        markersDirty = true;
        lastRenderAt = now;
        var missing = renderFromCache(need);
        if (missing.length && !fetchTimer && !activeAbort) {
          scheduleFetch(false);
        }
      }
    }
    starMaterial.uniforms.now.value = starClock();
    updateStarEye();
    updateStarShares();
    if (regionLayer) {
      regionLayer.update(null, focalPixels(), regionsVersion);
      updateRegionFrame();
      // MAP.357: the region points follow the camera every frame it moves (a few projections), and the layers' switches.
      var poiSignature = viewSignature + "|" + regionsVersion + "|" + regionsOn + "|" + hiddenCell;
      if (poiSignature !== lastPoiSignature) {
        lastPoiSignature = poiSignature;
        refreshRegionPois();
      }
    }
    if (glow) {
      glow.update(drawnSettled ? function () {
        return weightCalculator(camera.position.toArray(), NODE_EDGE_PC, focalPixels(), drawnSettled);
      } : null, focalPixels(), starSetVersion);
    }
    if (nebulaLayer) {
      nebulaLayer.update(drawnSettled ? function () {
        return weightCalculator(camera.position.toArray(), NODE_EDGE_PC, focalPixels(), drawnSettled);
      } : null, focalPixels(), nebulaVersion);
    }
    updateMagnitudeLimit(now);
    updateClouds();
    selectionRing.update();
    hoverRing.update();
    updateCourseMarkers();
    if (markerLayer) {
      // MAP.321: the symbols change level when the hand-off turns (or the open system's bodies change).
      var bodiesOn = bodyHandoff.isOn();
      var handoffSource = stageView.systemMarkerSource();
      if (handoffSource && handoffSource.entries !== lastBodyEntries) {
        lastBodyEntries = handoffSource.entries;    // a new scale rebuilds the bodies, a new system brings its own
        if (bodiesOn) markersDirty = true;
      }
      if (handoffSource && handoffSource.centre) {
        var handoffDistance = Math.hypot(camera.position.x - handoffSource.centre[0], camera.position.y - handoffSource.centre[1], camera.position.z - handoffSource.centre[2]);
        if (bodyHandoff.update(handoffDistance, handoffSource.radiusPc, true) !== bodiesOn) markersDirty = true;
      } else if (bodiesOn) {
        bodyHandoff.update(0, 0, false);
        markersDirty = true;
      }
      if (markersDirty) refreshMarkers();
      markerLayer.update();
    }
    if (sectorFrame) sectorFrame.update();
    if (galaxyLabels) galaxyLabels.update(now, false);
    renderFrame();
  })();
}

if (canvas && sceneData) {
  initGalaxyMap3d(canvas, sceneData);
}
