// static/generatelevers.js
//
// The Generate page's shape dials (web/templates/generate.html, GEN.606): a galaxy lighter than the halo thresholds has
// no dark halo or no stellar halo, so when the stellar mass is typed the Halo options that cannot apply are turned off
// and the reason is shown. The server refuses them as well (galaxy.shape_draw.plan_shape), so without JavaScript a
// plan with such a value stops with the same reason.

function update(mass) {
  const form = mass.form;
  const value = mass.value.trim() === "" ? null : Number(mass.value);
  const noDarkHalo = value !== null && Number.isFinite(value) && value < Number(mass.dataset.noDarkHaloBelow);
  const simple = value !== null && Number.isFinite(value) && value < Number(mass.dataset.simpleBelow);
  const switches = { halo_ratio: noDarkHalo, stellar_halo_share: simple };
  for (const [name, off] of Object.entries(switches)) {
    for (const input of form.querySelectorAll(`[data-lever="${name}"]`)) {
      input.disabled = off;
      if (off) input.value = "";
    }
  }
  for (const note of form.querySelectorAll("[data-halo-reason]")) {
    note.hidden = !(noDarkHalo || simple);
  }
}

for (const mass of document.querySelectorAll('input[data-lever="mass"]')) {
  mass.addEventListener("input", () => update(mass));
  update(mass);
}
