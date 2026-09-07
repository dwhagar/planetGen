// html/static/sectormap.js
//
// Click/keyboard handler for the isometric sector-map dots built by
// `lib/starmap.py`: populates the info side panel from the activated
// dot's `data-*` attributes instead of navigating straight to
// `system.py`, so activating a dot shows details first and the panel's
// own link is what navigates away. Built with plain DOM calls (never
// innerHTML/textContent-with-markup) since every data-* value is still
// database content (a system name can contain arbitrary characters via
// `--name`) -- consistent with the rest of `html/`'s dependency-free,
// no-build-step approach.
//
// TODO: gains pointer-drag-to-orbit (update rotateX/rotateY on the CSS
// 3D scene element) and wheel/pinch-to-zoom (a separate scale() on an
// outer wrapper) once `lib/starmap.py` switches to CSS 3D transforms --
// needs a pointerdown->pointerup movement-distance threshold so dragging
// the scene doesn't also fire showSystemInfo() below. See docs/TODO.md,
// "Near-term: interim `../src/html/` browser enhancements".

(function () {
  "use strict";

  function addField(dl, label, value) {
    if (!value) {
      return;
    }
    var dt = document.createElement("dt");
    dt.textContent = label;
    var dd = document.createElement("dd");
    dd.textContent = value;
    dl.appendChild(dt);
    dl.appendChild(dd);
  }

  function showSystemInfo(dot) {
    var panel = document.getElementById("starmap-info");
    if (!panel) {
      return;
    }
    panel.textContent = "";

    var heading = document.createElement("h3");
    heading.textContent = dot.dataset.name || "Unknown system";
    panel.appendChild(heading);

    var dl = document.createElement("dl");
    addField(dl, "Star type", dot.dataset.type);
    addField(dl, "Temperature", dot.dataset.temp);
    addField(dl, "Quadrant", dot.dataset.quadrant);
    addField(dl, "Location", dot.dataset.location);
    panel.appendChild(dl);

    var link = document.createElement("a");
    link.href = dot.dataset.href;
    link.className = "btn";
    link.textContent = "View system →";
    panel.appendChild(link);
  }

  document.querySelectorAll(".star-dot").forEach(function (dot) {
    dot.addEventListener("click", function () {
      showSystemInfo(dot);
    });
    dot.addEventListener("keydown", function (event) {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        showSystemInfo(dot);
      }
    });
  });
})();
