// html/static/phenomenonmap.js
//
// Real interactive zoom/pan for a stellar phenomenon's own AU-scale
// diagram (html/lib/phenomenonmap.py) -- wires the shared
// static/mapzoom.js viewBox zoom/pan onto #phenomenonmap-svg, plus the
// +/-/Reset buttons and a live "~X AU/ly across" scale readout.
//
// Unlike galaxymap.js (which converts a separate px-per-ly ratio),
// render_phenomenon_map_panel draws every shape in real AU coordinates
// directly -- this diagram's own SVG user units ARE astronomical units,
// one-to-one -- so the live viewBox width itself already IS the current
// "how many AU across" figure, no extra scale-factor lookup needed. No
// data-nav-target markers live inside this diagram (it's a passive
// to-scale drawing, not a clickable map), so mapzoom.js's own click-vs-
// drag gate has nothing to protect here -- it's simply inert.

(function () {
  "use strict";

  var AU_PER_LY = 63241.1;
  var LY_DISPLAY_THRESHOLD_AU = 1000; // above this many AU, show ly instead

  function formatSpan(auValue) {
    if (auValue >= LY_DISPLAY_THRESHOLD_AU) {
      var lyValue = auValue / AU_PER_LY;
      return (lyValue >= 10 ? lyValue.toFixed(1) : lyValue.toFixed(2)) + " ly";
    }
    if (auValue >= 10) {
      return Math.round(auValue).toLocaleString() + " AU";
    }
    return auValue.toFixed(2) + " AU";
  }

  function init() {
    var svg = document.getElementById("phenomenonmap-svg");
    var scaleEl = document.getElementById("phenomenonmap-scale");
    var controls = document.getElementById("phenomenonmap-controls");
    if (!svg || typeof window.planetgenInitSvgZoomPan !== "function") {
      return;
    }

    function updateScale(box) {
      if (!scaleEl) {
        return;
      }
      scaleEl.textContent = "≈ " + formatSpan(box.w) + " across";
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
