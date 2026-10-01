// html/static/generatebuttons.js
//
// The admin Generate buttons for a sector that isn't generated yet, shared
// by the Sector Map (sectormap.js, a neighbor) and the Galaxy Map
// (galaxymap3d.js, a sector cell). `target` is the server's
// {url, csrfField, csrfToken} (web/helpers.generate_target), given to the
// page only for a logged-in admin who can use the Generate page.

// One Generate button: a plain POST form to the admin Generate page
// (web/generate_page.py), which starts the job and redirects there, so it
// works like the page's own forms (CSRF token included) with no fetch.
export var NEIGHBORHOOD_RADIUS_PC = "30.7";

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

// The four buttons for sector (ring, layer, slot): it alone, its
// neighborhood, its column, and its whole shell (after a confirm).
export function generateButtons(target, ring, layer, slot) {
  var address = [["slot_ring", ring], ["slot_layer", layer], ["slot", slot]];
  var box = document.createElement("div");
  box.className = "starmap-generate";
  box.appendChild(generateForm(target, "Generate this sector", [["mode", "slot"]].concat(address)));
  box.appendChild(generateForm(target, "Generate neighborhood",
    [["mode", "slot"], ["slot_radius_pc", NEIGHBORHOOD_RADIUS_PC]].concat(address)));
  box.appendChild(generateForm(target, "Generate column",
    [["mode", "column"], ["column_ring", ring], ["column_slot", slot]]));
  box.appendChild(generateForm(target, "Generate the entire shell (not recommended)",
    [["mode", "shell"], ["shell_ring", ring], ["whole_shell", "1"]],
    "This generates every sector of ring " + ring + " through every layer, often thousands of " +
    "sectors, and can run for hours. Start it anyway?"));
  var hint = document.createElement("p");
  hint.className = "hint";
  hint.textContent = "Starts a background job on the Generate page. The neighborhood reaches about 100 ly.";
  box.appendChild(hint);
  return box;
}
