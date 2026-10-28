// html/static/facilityform.js
//
// The system page's "Place a facility" form (web/templates/system.html,
// web/system_facilities.py). Progressive enhancement over a plain POST
// form, which the server checks in full either way:
//
// - Placement decides the Host list: a host shows only when its
//   `data-placements` holds the chosen placement
//   (planetgen/population/facilities.py `host_placements`).
// - The orbit slider shows only for "in orbit", and reads out the
//   distance, period and speed as it moves. Its steps run on a log scale
//   from the host's `data-min-km` to `data-max-km` (just above the
//   surface to the edge of its sphere of influence), the browser mirror
//   of `facilities.distance_from_step`; the period is Kepler's third law,
//   as `facilities.orbit_for` works it out.
// - The Type list holds only the types of the chosen Kind (`data-kind`), and
//   choosing a type shows that type's fields (`data-facility-type-fields`).
//
// The three lists are Shoelace `sl-select`s (UX.126, templates/partials/ui_fields.html), so the script waits for them
// to be defined, listens to `sl-change`, and reads their `sl-option`s (a hidden option is taken out of the list by
// style.css, a disabled one cannot be picked).

// This module's own `?v=`, so the sibling loads at the same version.
const VERSION_QUERY = new URL(import.meta.url).search;
const { fieldsReady } = await import(`./fieldsready.js${VERSION_QUERY}`);
const { formatDistanceKm, formatDistanceKmSpoken, formatDurationSeconds, formatSpeedKms } = await import(`./unitladder.js${VERSION_QUERY}`);

await fieldsReady();

const G = 6.6743e-11; // m^3 kg^-1 s^-2

function optionsOf(select) {
  return Array.from(select.querySelectorAll("sl-option"));
}

function chosenOption(select) {
  return optionsOf(select).find((option) => option.value === select.value) || null;
}

function offer(option, fits) {
  option.hidden = !fits;
  option.disabled = !fits;
}

function distanceFromStep(step, steps, lowestKm, highestKm) {
  const fraction = Math.min(Math.max(step, 0), steps) / steps;
  return lowestKm * Math.pow(highestKm / lowestKm, fraction);
}

function setUpTypes(form) {
  const kind = form.querySelector("[data-facility-kind]");
  const type = form.querySelector("[data-facility-type]");
  if (!kind || !type) {
    return;
  }
  const sets = form.querySelectorAll("[data-facility-type-fields]");
  const show = () => {
    for (const fields of sets) {
      fields.hidden = fields.dataset.facilityTypeFields !== type.value;
    }
  };
  const filterTypes = () => {
    for (const option of optionsOf(type)) {
      offer(option, !option.dataset.kind || option.dataset.kind === kind.value);
    }
    const chosen = chosenOption(type);
    if (chosen && chosen.disabled) {
      type.value = "";
    }
    show();
  };
  kind.addEventListener("sl-change", filterTypes);
  type.addEventListener("sl-change", show);
  filterTypes();
}

function setUp(form) {
  setUpTypes(form);
  const placement = form.querySelector("[data-facility-placement]");
  const host = form.querySelector("[data-facility-host]");
  const orbit = form.querySelector("[data-facility-orbit]");
  const slider = orbit && orbit.querySelector('input[type="range"]');
  const readout = orbit && orbit.querySelector("[data-facility-orbit-readout]");
  if (!placement || !host || !orbit || !slider || !readout) {
    return;
  }

  const showReadout = () => {
    const option = chosenOption(host);
    if (!option || !option.dataset.minKm) {
      readout.textContent = "";
      return;
    }
    const km = distanceFromStep(Number(slider.value), Number(slider.max),
      Number(option.dataset.minKm), Number(option.dataset.maxKm));
    const meters = km * 1000;
    const periodS = 2 * Math.PI * Math.sqrt(meters ** 3 / (G * Number(option.dataset.massKg)));
    const speedKms = (2 * Math.PI * km) / periodS;
    readout.textContent = `${formatDistanceKm(km)} from its host, one orbit every ${formatDurationSeconds(periodS)}, `
      + `at ${formatSpeedKms(speedKms)}`;
    slider.setAttribute("aria-valuetext", formatDistanceKmSpoken(km));
  };

  const filterHosts = () => {
    const chosen = placement.value;
    let firstShown = null;
    for (const option of optionsOf(host)) {
      const fits = (option.dataset.placements || "").split(" ").includes(chosen);
      offer(option, fits);
      if (fits && firstShown === null) {
        firstShown = option;
      }
    }
    const selected = chosenOption(host);
    if ((!selected || selected.disabled) && firstShown) {
      host.value = firstShown.value;
    }
    orbit.hidden = chosen !== "orbital";
    showReadout();
  };

  placement.addEventListener("sl-change", filterHosts);
  host.addEventListener("sl-change", showReadout);
  slider.addEventListener("input", showReadout);
  filterHosts();
}

for (const form of document.querySelectorAll(".facility-form")) {
  setUp(form);
}
