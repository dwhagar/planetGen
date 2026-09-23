// html/static/galaxymap.js
//
// Real interactive zoom/pan for the Galaxy Map (html/lib/galaxymap.py) --
// wires the shared static/mapzoom.js viewBox zoom/pan onto
// #galaxymap-svg, plus the +/-/Reset buttons and a live "~X ly across"
// scale readout. Marker/label clicks (data-nav-target, intercepted by
// static/navform.js) still work unchanged -- mapzoom.js's own drag-vs-
// click gate is what keeps a drag-pan from also firing a stray
// navigation click on whatever dot the pointer happened to release over.
//
// Everything about WHERE dots/labels sit (server-computed px_per_ly, real
// vs. illustrative density, Quadrant/Ring geometry) is unchanged --
// zooming only changes how much of that same fixed drawing is currently
// visible, the same "viewBox is the whole camera" approach
// static/mapzoom.js's own module docstring explains.

(function () {
  "use strict";

  function formatLy(value) {
    if (value >= 100) {
      return Math.round(value).toLocaleString();
    }
    if (value >= 10) {
      return value.toFixed(1);
    }
    return value.toFixed(2);
  }

  function init() {
    var svg = document.getElementById("galaxymap-svg");
    var scaleEl = document.getElementById("galaxymap-scale");
    var controls = document.getElementById("galaxymap-controls");
    if (!svg || typeof window.planetgenInitSvgZoomPan !== "function") {
      return;
    }

    var pxPerLy = parseFloat(svg.dataset.pxPerLy) || 1;

    function updateScale(box) {
      if (!scaleEl) {
        return;
      }
      var spanLy = box.w / pxPerLy;
      scaleEl.textContent = "≈ " + formatLy(spanLy) + " ly across";
    }

    var zoomInBtn = controls ? controls.querySelector('[data-action="zoom-in"]') : null;
    var zoomOutBtn = controls ? controls.querySelector('[data-action="zoom-out"]') : null;
    var resetBtn = controls ? controls.querySelector('[data-action="reset"]') : null;

    var controller = window.planetgenInitSvgZoomPan(svg, {
      zoomInBtn: zoomInBtn,
      zoomOutBtn: zoomOutBtn,
      resetBtn: resetBtn,
      onChange: updateScale,
    });

    updateScale(controller.getBox());
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
