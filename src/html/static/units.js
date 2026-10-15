// html/static/units.js
//
// Shows every server-rendered measurement in the reader's unit system
// (UX.78). The server writes each one as <span class="qty" data-q="length"
// data-si="384400000">384,400 km</span> with the Automatic text
// (planetgen/web/lib/fmt.py `quantity_html`); this re-renders it from the SI
// value with the shared unit ladders, and runs the "Units" button in the
// settings menu that cycles Automatic, Metric only and Customary. Without
// script the page reads correctly in Automatic. Panels built by other
// scripts read the saved preference through unitladder.js when they render.

const VERSION_QUERY = new URL(import.meta.url).search;
const { LADDERS, formatQuantity } = await import(`./unitladder.js${VERSION_QUERY}`);
const { LABELS, SYSTEMS, setUnitSystem, unitSystem } = await import(`./unitpreference.js${VERSION_QUERY}`);

function render(root, system) {
  for (const node of (root || document).querySelectorAll(".qty[data-q][data-si]")) {
    const ladder = LADDERS[node.dataset.q];
    const value = Number(node.dataset.si);
    if (ladder && Number.isFinite(value)) {
      node.textContent = formatQuantity(ladder, value, { companions: !node.hasAttribute("data-nc"), system });
    }
  }
}

function label(button, system) {
  button.textContent = `Units: ${LABELS[system]}`;
}

function wire() {
  const system = unitSystem();
  document.documentElement.dataset.units = system;
  if (system !== "auto") {
    render(document, system);
  }
  for (const button of document.querySelectorAll("[data-units-toggle]")) {
    label(button, system);
    button.hidden = false;
    button.addEventListener("click", () => {
      const next = setUnitSystem(SYSTEMS[(SYSTEMS.indexOf(unitSystem()) + 1) % SYSTEMS.length]);
      document.documentElement.dataset.units = next;
      render(document, next);
      for (const other of document.querySelectorAll("[data-units-toggle]")) {
        label(other, next);
      }
      document.dispatchEvent(new CustomEvent("planetgen:units", { detail: { system: next } }));
    });
  }
}

wire();
// Quantities added after load (a panel a script builds) call this.
window.planetgenRenderQuantities = (root) => render(root, unitSystem());
