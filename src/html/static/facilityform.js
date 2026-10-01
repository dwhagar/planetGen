// html/static/facilityform.js
//
// The system page's "Place a facility" form (web/templates/system.html,
// web/system_facilities.py). Progressive enhancement over a plain POST
// form, which the server checks in full either way:
//
// - Placement decides the Host list: a host shows only when its
//   `data-placements` holds the chosen placement
//   (stellarObjects/facilities.py `host_placements`).
// - The orbit slider shows only for "in orbit", and reads out the
//   distance, period and speed as it moves. Its steps run on a log scale
//   from the host's `data-min-km` to `data-max-km` (just above the
//   surface to the edge of its sphere of influence), the browser mirror
//   of `facilities.distance_from_step`; the period is Kepler's third law,
//   as `facilities.orbit_for` works it out.

// This module's own `?v=`, so the sibling loads at the same version.
const VERSION_QUERY = new URL(import.meta.url).search;
const { formatDistanceKm } = await import(`./distance.js${VERSION_QUERY}`);

const G = 6.6743e-11; // m^3 kg^-1 s^-2
const SECONDS_PER_DAY = 86400;
const DAYS_PER_YEAR = 365.25;

function distanceFromStep(step, steps, lowestKm, highestKm) {
  const fraction = Math.min(Math.max(step, 0), steps) / steps;
  return lowestKm * Math.pow(highestKm / lowestKm, fraction);
}

function formatPeriod(seconds) {
  const days = seconds / SECONDS_PER_DAY;
  const fmt = (value) => value.toLocaleString("en-US", { maximumSignificantDigits: 3 });
  if (days < 1) {
    return `${fmt(days * 24)} hours`;
  }
  if (days < DAYS_PER_YEAR) {
    return `${fmt(days)} days`;
  }
  return `${fmt(days / DAYS_PER_YEAR)} years`;
}

function setUp(form) {
  const placement = form.querySelector("[data-facility-placement]");
  const host = form.querySelector("[data-facility-host]");
  const orbit = form.querySelector("[data-facility-orbit]");
  const slider = orbit && orbit.querySelector('input[type="range"]');
  const readout = orbit && orbit.querySelector("[data-facility-orbit-readout]");
  if (!placement || !host || !orbit || !slider || !readout) {
    return;
  }

  const showReadout = () => {
    const option = host.selectedOptions[0];
    if (!option || !option.dataset.minKm) {
      readout.textContent = "";
      return;
    }
    const km = distanceFromStep(Number(slider.value), Number(slider.max),
      Number(option.dataset.minKm), Number(option.dataset.maxKm));
    const meters = km * 1000;
    const periodS = 2 * Math.PI * Math.sqrt(meters ** 3 / (G * Number(option.dataset.massKg)));
    const speedKms = (2 * Math.PI * km) / periodS;
    readout.textContent = `${formatDistanceKm(km)} from its host, one orbit every ${formatPeriod(periodS)}, `
      + `at ${speedKms.toLocaleString("en-US", { minimumFractionDigits: 2, maximumFractionDigits: 2 })} km/s`;
    slider.setAttribute("aria-valuetext", formatDistanceKm(km));
  };

  const filterHosts = () => {
    const chosen = placement.value;
    let firstShown = null;
    for (const option of host.options) {
      const fits = (option.dataset.placements || "").split(" ").includes(chosen);
      option.hidden = !fits;
      option.disabled = !fits;
      if (fits && firstShown === null) {
        firstShown = option;
      }
    }
    const selected = host.selectedOptions[0];
    if ((!selected || selected.disabled) && firstShown) {
      firstShown.selected = true;
    }
    orbit.hidden = chosen !== "orbital";
    showReadout();
  };

  placement.addEventListener("change", filterHosts);
  host.addEventListener("change", showReadout);
  slider.addEventListener("input", showReadout);
  filterHosts();
}

for (const form of document.querySelectorAll(".facility-form")) {
  setUp(form);
}
