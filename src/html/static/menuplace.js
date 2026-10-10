// html/static/menuplace.js
//
// Button menus that open where they can be seen (UX.85). A <details>
// menu whose panel hangs below its button (the Galaxy Map's Menu and Steps,
// the Bookmarks list) used to open out of the bottom of the window, so the
// visitor had to scroll to read it. When one opens, `placeMenu` measures its
// panel and: moves it above the button when there is more room there, caps
// its height to the room that is left and lets it scroll inside, and slides
// it sideways back into the window; a menu that sits in the page's flow is
// scrolled into view. The menus are found by the selectors below, wherever
// they are on the page; every other dropdown is a Shoelace <sl-dropdown>,
// which flips by itself.

const MARGIN = 8;
const MIN_HEIGHT = 120;
const GAP = 6;

// A menu's <details> selector and the selector of its panel inside it.
export const MENUS = [
  ["details.galaxy-menu", ".galaxy-menu-panel"],
  ["details.galaxy-steps", ".galaxy-steps-panel"],
  ["details.bookmarks-menu", ".bookmarks-panel"],
];

function clearPlacement(panel) {
  panel.style.removeProperty("top");
  panel.style.removeProperty("bottom");
  panel.style.removeProperty("max-height");
  panel.style.removeProperty("translate");
  panel.style.removeProperty("margin-top");
  panel.style.removeProperty("margin-bottom");
  panel.style.removeProperty("overflow-y");
}

// Puts `panel` (opened from `menu`) inside the visible window.
export function placeMenu(menu, panel) {
  clearPlacement(panel);
  const view = window.visualViewport;
  const height = view ? view.height : window.innerHeight;
  const width = view ? view.width : window.innerWidth;
  const rect = panel.getBoundingClientRect();
  const button = (menu.querySelector("summary") || menu).getBoundingClientRect();
  const positioned = window.getComputedStyle(panel).position === "absolute";

  if (!positioned) {
    // In the page's flow (the Bookmarks list): scroll it into view.
    if (rect.bottom > height - MARGIN || rect.top < MARGIN) {
      panel.scrollIntoView({ block: "nearest", behavior: "instant" });
    }
    return;
  }

  if (rect.bottom > height - MARGIN) {
    const below = height - MARGIN - button.bottom;
    const above = button.top - MARGIN;
    panel.style.setProperty("overflow-y", "auto");
    if (above > below && above >= Math.min(rect.height, MIN_HEIGHT * 2)) {
      // More room above the button: hang the panel above it. `top` is
      // worked out against the panel's own containing block, which is the
      // controls row on a narrow screen and not the button.
      const room = Math.max(above - GAP, MIN_HEIGHT);
      const used = Math.min(rect.height, room);
      const block = panel.offsetParent ? panel.offsetParent.getBoundingClientRect() : { top: 0 };
      panel.style.setProperty("margin-top", "0");
      panel.style.setProperty("bottom", "auto");
      panel.style.setProperty("max-height", `${room}px`);
      panel.style.setProperty("top", `${button.top - GAP - used - block.top}px`);
    } else {
      panel.style.setProperty("max-height", `${Math.max(below - GAP, MIN_HEIGHT)}px`);
    }
  }

  const placed = panel.getBoundingClientRect();
  let shift = 0;
  if (placed.right > width - MARGIN) shift = width - MARGIN - placed.right;
  if (placed.left + shift < MARGIN) shift = MARGIN - placed.left;
  if (shift) panel.style.setProperty("translate", `${shift}px 0`);
}

function wire(menu, panelSelector) {
  const panel = menu.querySelector(panelSelector);
  if (!panel) return;
  menu.addEventListener("toggle", function () {
    if (menu.open) {
      placeMenu(menu, panel);
      menu.dataset.placed = "1";      // the browser tests wait for this
    } else {
      delete menu.dataset.placed;
    }
  });
}

function wireAll() {
  MENUS.forEach(function ([menuSelector, panelSelector]) {
    document.querySelectorAll(menuSelector).forEach(function (menu) {
      wire(menu, panelSelector);
    });
  });
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", wireAll);
} else {
  wireAll();
}
