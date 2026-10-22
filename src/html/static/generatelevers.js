// static/generatelevers.js
//
// The Generate page's shape dials (web/templates/generate.html, GEN.606): a galaxy lighter than the halo thresholds has
// no dark halo or no stellar halo, so when the stellar mass is typed the Halo options that cannot apply are turned off
// and the reason is shown. The server refuses them as well (galaxy.shape_draw.plan_shape), so without JavaScript a
// plan with such a value stops with the same reason.
//
// GEN.617: under the dials a line says how dense the galaxy is and how much of it the outline holds, the numbers
// `planetgen plan` prints (the server works them out with the same code, `generation.plan_outline`), with a warning
// when the disc is thin and a refusal when the outline would be empty.

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

function showReport(note, answer) {
  note.replaceChildren();
  const add = (text, className) => {
    const line = document.createElement("span");
    line.textContent = text;
    if (className) line.className = className;
    note.append(line, document.createElement("br"));
  };
  if (answer.error) return add(answer.error, "error");
  add(answer.line + (answer.sample ? " (the dials left on Auto are drawn from a sample seed here)" : ""));
  for (const warning of answer.warnings) add(`Warning: ${warning}`);
  if (answer.refusal) add(`Refused: ${answer.refusal}`, "error");
}

for (const note of document.querySelectorAll("[data-shape-report]")) {
  const form = note.closest("form");
  let timer = null;
  let sent = 0;
  const ask = async () => {
    const query = new URLSearchParams();
    for (const input of form.querySelectorAll('[data-lever], select[name$="_unit"]')) {
      if (!input.disabled && input.name && input.value.trim() !== "") query.set(input.name, input.value.trim());
    }
    const mine = ++sent;
    try {
      const response = await fetch(`${note.dataset.url}?${query}`, { headers: { Accept: "application/json" } });
      const answer = await response.json();
      if (mine === sent) showReport(note, answer);
    } catch (error) {
      if (mine === sent) note.textContent = "";
    }
  };
  form.addEventListener("input", () => { clearTimeout(timer); timer = setTimeout(ask, 400); });
  form.addEventListener("change", () => { clearTimeout(timer); timer = setTimeout(ask, 400); });
  ask();
}
