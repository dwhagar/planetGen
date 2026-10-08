// html/static/navtimes.js
//
// The NAV result page's travel-times table (UX.70): one collapsed section
// with a Warp/Fold switch. The server renders both tables and a hidden
// switch; importing this module shows the switch and one table at a time
// (no inline script, see the Content-Security-Policy in web/__init__.py).
// Without script both tables stay visible, one under the other.

for (const group of document.querySelectorAll("[data-travel-times]")) {
  const buttons = [...group.querySelectorAll("[data-travel-mode]")];
  const tables = [...group.querySelectorAll("[data-travel-table]")];
  const show = (mode) => {
    for (const table of tables) table.hidden = table.dataset.travelTable !== mode;
    for (const button of buttons) button.setAttribute("aria-pressed", String(button.dataset.travelMode === mode));
  };
  for (const button of buttons) button.addEventListener("click", () => show(button.dataset.travelMode));
  const switcher = group.querySelector("[data-travel-switch]");
  if (switcher && buttons.length) {
    switcher.hidden = false;
    show(buttons[0].dataset.travelMode);
  }
}
