// html/static/generatebuttons.js
//
// The admin Generate buttons for a sector that isn't generated yet, shared
// by the Sector Map (sectormap.js, a neighbor) and the Galaxy Map
// (galaxymap3d.js, a sector cell). `target` is the server's
// {url, csrfField, csrfToken} (web/helpers.generate_target), given to the
// page only for a logged-in admin who can use the Generate page.

const VERSION_QUERY = new URL(import.meta.url).search;

const { LIGHTYEAR_M, PARSEC_M } = await import(`./distance.js${VERSION_QUERY}`);

// The neighborhood dialog's radius in light years: the default, and the
// bounds generationLimits allows (a sector's own edge up to
// MAX_GENERATE_RADIUS_PC, about 652 ly).
export var NEIGHBORHOOD_DEFAULT_LY = 100;
export var NEIGHBORHOOD_MIN_LY = 13;
export var NEIGHBORHOOD_MAX_LY = 652;
// More sectors than this asks before it starts.
export var NEIGHBORHOOD_CONFIRM_SECTORS = 5000;

export function lyToPc(ly) {
  return (ly * LIGHTYEAR_M) / PARSEC_M;
}

// About how many sector slots a sphere of `radiusLy` covers, at sectors
// `edgeLy` light years a side: the sphere's volume over a sector's. An
// estimate for the dialog, not a count -- the real job skips slots the
// galaxy's density leaves out, so it is an upper bound.
export function neighborhoodSectors(radiusLy, edgeLy) {
  if (!(edgeLy > 0) || !(radiusLy > 0)) {
    return null;
  }
  return Math.max(1, Math.round(((4 / 3) * Math.PI * Math.pow(radiusLy, 3)) / Math.pow(edgeLy, 3)));
}

// One Generate button: a plain POST form to the admin Generate page
// (web/generate_page.py), which starts the job and redirects there, so it
// works like the page's own forms (CSRF token included) with no fetch.
export function generateForm(target, label, fields, confirmText) {
  var form = document.createElement("form");
  form.method = "post";
  form.action = target.url;
  form.className = "starmap-generate-form";
  var all = [[target.csrfField, target.csrfToken], ["action", "galaxy"]].concat(fields);
  all.forEach(function (pair) {
    var input = document.createElement("input");
    input.type = "hidden";
    input.name = pair[0];
    input.value = String(pair[1]);
    form.appendChild(input);
  });
  var button = document.createElement("button");
  button.type = "submit";
  button.className = "btn";
  button.textContent = label;
  form.appendChild(button);
  if (confirmText) {
    form.addEventListener("submit", function (event) {
      if (!window.confirm(confirmText)) {
        event.preventDefault();
      }
    });
  }
  return form;
}

// "Generate neighborhood", which asks for the radius in light years
// first (the drill-down design's section 6): a number field between
// NEIGHBORHOOD_MIN_LY and NEIGHBORHOOD_MAX_LY, with about how many
// sectors that reaches when the sector size is known (`edgeLy`), and a
// confirmation past NEIGHBORHOOD_CONFIRM_SECTORS.
export function neighborhoodForm(target, ring, layer, slot, edgeLy) {
  var form = generateForm(target, "Generate neighborhood",
    [["mode", "slot"], ["slot_ring", ring], ["slot_layer", layer], ["slot", slot]]);
  var radiusPc = document.createElement("input");
  radiusPc.type = "hidden";
  radiusPc.name = "slot_radius_pc";
  form.insertBefore(radiusPc, form.firstChild);

  var field = document.createElement("label");
  field.className = "starmap-generate-radius";
  var caption = document.createElement("span");
  caption.textContent = "Radius (light years)";
  var input = document.createElement("input");
  input.type = "number";
  input.min = String(NEIGHBORHOOD_MIN_LY);
  input.max = String(NEIGHBORHOOD_MAX_LY);
  input.step = "1";
  input.value = String(NEIGHBORHOOD_DEFAULT_LY);
  field.appendChild(caption);
  field.appendChild(input);
  var reach = document.createElement("span");
  reach.className = "hint";
  form.insertBefore(field, form.lastChild);
  form.insertBefore(reach, form.lastChild);

  function radiusLy() {
    return Math.min(NEIGHBORHOOD_MAX_LY, Math.max(NEIGHBORHOOD_MIN_LY, Number(input.value) || 0));
  }

  function update() {
    var ly = radiusLy();
    var sectors = neighborhoodSectors(ly, edgeLy);
    radiusPc.value = lyToPc(ly).toFixed(3);
    reach.textContent = sectors == null
      ? "Reaches about " + Math.round(lyToPc(ly)) + " pc"
      : "Up to about " + sectors.toLocaleString() + " sector" + (sectors === 1 ? "" : "s");
  }

  input.addEventListener("input", update);
  input.addEventListener("change", update);
  update();

  form.addEventListener("submit", function (event) {
    update();
    var sectors = neighborhoodSectors(radiusLy(), edgeLy);
    if (sectors != null && sectors > NEIGHBORHOOD_CONFIRM_SECTORS
        && !window.confirm("A " + radiusLy() + " light year neighborhood covers up to about "
          + sectors.toLocaleString() + " sectors, which can run for hours. Start it anyway?")) {
      event.preventDefault();
    }
  });
  return form;
}

// The four buttons for sector (ring, layer, slot): it alone, its
// neighborhood (with a radius to set), its column, and its whole shell
// (after a confirm). `edgeLy` (the sector edge, light years) is optional:
// without it the neighborhood shows its reach in parsecs rather than
// about how many sectors it covers.
export function generateButtons(target, ring, layer, slot, edgeLy) {
  var address = [["slot_ring", ring], ["slot_layer", layer], ["slot", slot]];
  var box = document.createElement("div");
  box.className = "starmap-generate";
  box.appendChild(generateForm(target, "Generate this sector", [["mode", "slot"]].concat(address)));
  box.appendChild(neighborhoodForm(target, ring, layer, slot, edgeLy));
  box.appendChild(generateForm(target, "Generate column",
    [["mode", "column"], ["column_ring", ring], ["column_slot", slot]]));
  box.appendChild(generateForm(target, "Generate the entire shell (not recommended)",
    [["mode", "shell"], ["shell_ring", ring], ["whole_shell", "1"]],
    "This generates every sector of ring " + ring + " through every layer, often thousands of " +
    "sectors, and can run for hours. Start it anyway?"));
  var hint = document.createElement("p");
  hint.className = "hint";
  hint.textContent = "Starts a background job on the Generate page.";
  box.appendChild(hint);
  return box;
}
