// html/static/sectormap.js
//
// Drag-to-rotate, scroll/button-to-zoom, and click/keyboard-for-info
// behavior for the 3D sector map built by `lib/starmap.py`. The actual
// rotation math is just two numbers (rotateX/rotateY degrees) fed to
// `#starmap-scene`'s CSS transform -- the browser's own `preserve-3d`
// compositor does the real 3D projection and occlusion, this file only
// tracks drag distance and turns it into degrees. Clicking a dot
// populates the info side panel from its `data-*` attributes instead of
// navigating straight to `system.py`, so a click shows details first and
// the panel's own link is what navigates away. Built with plain DOM
// calls (never innerHTML/textContent-with-markup) since every data-*
// value is still database content (a system name can contain arbitrary
// characters via `--name`) -- consistent with the rest of `html/`'s
// dependency-free, no-build-step approach.
//
// Also keeps the "Galactic Center" compass label facing the camera
// (like a star dot, but selected via the broader `.billboard` class so
// the compass arrow's own straight line -- which must NOT billboard, its
// whole point is showing a real 3D direction -- is left alone), and
// keeps the scale-bar legend (`#starmap-scale-bar`/`#starmap-scale-label`)
// showing the sector's true physical scale at the current zoom level.

(function () {
  "use strict";

  var DEFAULT_ROTATE_X = -18;
  var DEFAULT_ROTATE_Y = -32;
  var ROTATE_SENSITIVITY = 0.4; // degrees per pixel of drag
  var KEY_ROTATE_STEP = 6; // degrees per arrow-key press
  var MIN_ZOOM = 0.5;
  var MAX_ZOOM = 2.5;
  var ZOOM_STEP = 0.15;
  var WHEEL_ZOOM_STEP = 0.08;
  // Below this much total pointer movement, a pointerdown->pointerup on a
  // dot still counts as a click (not a drag) -- without this, the native
  // click event fires for a small in-place jitter too, which is fine, but
  // fires just as readily after the user actually dragged the scene
  // around and happened to release over a dot, which is not fine.
  var DRAG_CLICK_THRESHOLD_PX = 4;

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
    addField(dl, "Octant", dot.dataset.quadrant);
    addField(dl, "Location", dot.dataset.location);
    panel.appendChild(dl);

    var link = document.createElement("a");
    link.href = dot.dataset.href;
    link.className = "btn";
    link.textContent = "View system →";
    panel.appendChild(link);
  }

  function initStarmap(stage) {
    var scene = document.getElementById("starmap-scene");
    var zoomWrap = document.getElementById("starmap-zoom");
    var controls = document.getElementById("starmap-controls");
    if (!scene || !zoomWrap) {
      return;
    }

    var rotateX = DEFAULT_ROTATE_X;
    var rotateY = DEFAULT_ROTATE_Y;
    var zoom = 1;
    var dragging = false;
    var dragDistance = 0;
    var lastClientX = 0;
    var lastClientY = 0;
    // Set true only for the single click event immediately following an
    // over-threshold drag release, then cleared either by that click
    // handler or (if no click follows) a same-tick timeout -- so it can
    // never wedge a *later*, unrelated click closed. Checking
    // `dragDistance` directly in the click handler instead would leave a
    // stale nonzero value in place for any click that arrives without a
    // pointerdown of its own (assistive tech, a programmatically
    // dispatched click), wrongly suppressing it.
    var suppressNextClick = false;
    // Click/keyboard hit-testing only ever targets a star -- kept as its
    // own narrower list (not every `.billboard`) so the compass's
    // "Galactic Center" label near it can't be mistaken for a star dot
    // in `dotAtPoint` below or pick up its own keydown handler.
    var dots = Array.prototype.slice.call(scene.querySelectorAll(".star-dot"));
    // Everything on the map that must keep facing the camera regardless
    // of `.starmap-scene`'s own rotation -- star dots and the compass's
    // text label alike (billboarding a *line*, like the compass arrow
    // itself or a wedge edge, would be wrong: a line's whole visual
    // point is showing its real 3D direction, not facing the viewer).
    var billboards = Array.prototype.slice.call(scene.querySelectorAll(".billboard"));

    var scaleEl = document.getElementById("starmap-scale");
    var scaleBarEl = document.getElementById("starmap-scale-bar");
    var scaleLabelEl = document.getElementById("starmap-scale-label");
    // Exact ly-per-pixel ratio at zoom=1, computed server-side
    // (`lib/starmap.py`'s `_ly_per_px_at_zoom_1`) from the sector's own
    // real `edge_mpc` -- dividing by the live `zoom` factor below keeps
    // the bar's label true to the actual current scale as the user zooms,
    // rather than a value that was only ever right at the default zoom.
    var lyPerPxAtZoom1 = scaleEl ? parseFloat(scaleEl.dataset.lyPerPx) : 0;
    var SCALE_BAR_TARGET_PX = 70;

    // Snaps an arbitrary positive value to the nearest "nice" 1/2/5 * 10^n
    // -- the standard map-scale-bar convention, so the label reads "5 ly"
    // or "20 ly" rather than an ugly "6.283 ly".
    function niceScaleValue(raw) {
      if (!isFinite(raw) || raw <= 0) {
        return 0;
      }
      var magnitude = Math.pow(10, Math.floor(Math.log10(raw)));
      var mantissa = raw / magnitude;
      var niceMantissa;
      if (mantissa < 1.5) niceMantissa = 1;
      else if (mantissa < 3.5) niceMantissa = 2;
      else if (mantissa < 7.5) niceMantissa = 5;
      else niceMantissa = 10;
      return niceMantissa * magnitude;
    }

    function formatLy(value) {
      if (value >= 100) return Math.round(value) + " ly";
      if (value >= 1) return Math.round(value * 10) / 10 + " ly";
      return Math.round(value * 1000) / 1000 + " ly";
    }

    function updateScaleBar() {
      if (!scaleEl || !scaleBarEl || !scaleLabelEl || !lyPerPxAtZoom1) {
        return;
      }
      var lyPerPx = lyPerPxAtZoom1 / zoom;
      var niceLy = niceScaleValue(SCALE_BAR_TARGET_PX * lyPerPx);
      if (!niceLy) {
        return;
      }
      scaleBarEl.style.width = (niceLy / lyPerPx).toFixed(1) + "px";
      scaleLabelEl.textContent = formatLy(niceLy);
    }

    function apply() {
      scene.style.transform = "rotateX(" + rotateX + "deg) rotateY(" + rotateY + "deg)";
      zoomWrap.style.transform = "scale(" + zoom.toFixed(2) + ")";
      // Billboarding: counter-rotate each dot by the algebraic inverse of
      // the scene's own rotation (reverse function order, negated angles)
      // so it keeps facing the camera instead of going edge-on as the
      // scene turns -- see the docstring on lib/starmap.py's `_dot_html`.
      // Only correct because `.starmap-stage` has no `perspective`: with
      // one, this composition stops being pure rotation and a plain
      // inverse no longer cancels it (confirmed the hard way -- it only
      // matched at the one rotation angle it happened to be tested at).
      var counterRotate = "rotateY(" + -rotateY + "deg) rotateX(" + -rotateX + "deg)";
      billboards.forEach(function (billboard) {
        billboard.style.transform = counterRotate;
      });
      updateScaleBar();
    }

    function setZoom(value) {
      zoom = Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, value));
      apply();
    }

    function resetView() {
      rotateX = DEFAULT_ROTATE_X;
      rotateY = DEFAULT_ROTATE_Y;
      zoom = 1;
      apply();
    }

    apply();

    stage.addEventListener("pointerdown", function (event) {
      dragging = true;
      dragDistance = 0;
      lastClientX = event.clientX;
      lastClientY = event.clientY;
      try {
        stage.setPointerCapture(event.pointerId);
      } catch (err) {
        // Pointer capture isn't essential -- dragging still works via
        // ordinary pointermove bubbling if the browser refuses it (e.g.
        // a pointerId that's already gone).
      }
    });

    stage.addEventListener("pointermove", function (event) {
      if (!dragging) {
        return;
      }
      var deltaX = event.clientX - lastClientX;
      var deltaY = event.clientY - lastClientY;
      dragDistance += Math.abs(deltaX) + Math.abs(deltaY);
      lastClientX = event.clientX;
      lastClientY = event.clientY;

      rotateY += deltaX * ROTATE_SENSITIVITY;
      rotateX = Math.max(-89, Math.min(89, rotateX - deltaY * ROTATE_SENSITIVITY));
      apply();
    });

    function endDrag(event) {
      dragging = false;
      if (dragDistance > DRAG_CLICK_THRESHOLD_PX) {
        suppressNextClick = true;
        // Safety net: a pointerup isn't always followed by a click (e.g.
        // pointercancel) -- don't leave this suppressing some unrelated
        // later click if one never arrives to consume and clear it.
        setTimeout(function () {
          suppressNextClick = false;
        }, 0);
      }
      dragDistance = 0;
      try {
        stage.releasePointerCapture(event.pointerId);
      } catch (err) {
        // Already released/invalid -- nothing to clean up.
      }
    }
    stage.addEventListener("pointerup", endDrag);
    stage.addEventListener("pointercancel", endDrag);

    stage.addEventListener(
      "wheel",
      function (event) {
        event.preventDefault();
        setZoom(zoom + (event.deltaY < 0 ? WHEEL_ZOOM_STEP : -WHEEL_ZOOM_STEP));
      },
      { passive: false }
    );

    stage.addEventListener("keydown", function (event) {
      var key = event.key;
      if (key === "ArrowLeft" || key === "ArrowRight" || key === "ArrowUp" || key === "ArrowDown") {
        event.preventDefault();
        if (key === "ArrowLeft") rotateY -= KEY_ROTATE_STEP;
        if (key === "ArrowRight") rotateY += KEY_ROTATE_STEP;
        if (key === "ArrowUp") rotateX = Math.max(-89, rotateX - KEY_ROTATE_STEP);
        if (key === "ArrowDown") rotateX = Math.min(89, rotateX + KEY_ROTATE_STEP);
        apply();
      }
    });

    // Which dot (if any) sits under a given viewport point, found by
    // geometry (`getBoundingClientRect`) rather than the DOM's own hit-
    // testing (`elementFromPoint`/native click dispatch) -- both of those
    // turn out to disagree with each other for an element nested this
    // deep inside a rotated `preserve-3d` hierarchy (confirmed directly:
    // `elementFromPoint` and `elementsFromPoint()[0]` returned different
    // elements for the identical coordinate), so a real click can silently
    // miss a dot that's plainly sitting right there on screen. A dot's
    // `getBoundingClientRect()` remains reliable regardless (it reflects
    // actual rendered position, not hit-test routing), so clicks are
    // resolved against that instead. Ties (overlapping binary dots)
    // go to whichever center is nearest the click.
    function dotAtPoint(clientX, clientY) {
      var best = null;
      var bestDistance = Infinity;
      dots.forEach(function (dot) {
        var rect = dot.getBoundingClientRect();
        if (clientX < rect.left || clientX > rect.right || clientY < rect.top || clientY > rect.bottom) {
          return;
        }
        var centerX = rect.left + rect.width / 2;
        var centerY = rect.top + rect.height / 2;
        var distance = Math.hypot(clientX - centerX, clientY - centerY);
        if (distance < bestDistance) {
          bestDistance = distance;
          best = dot;
        }
      });
      return best;
    }

    stage.addEventListener("click", function (event) {
      if (suppressNextClick) {
        suppressNextClick = false;
        return;
      }
      var dot = dotAtPoint(event.clientX, event.clientY);
      if (dot) {
        showSystemInfo(dot);
      }
    });

    dots.forEach(function (dot) {
      dot.addEventListener("keydown", function (event) {
        if (event.key === "Enter" || event.key === " ") {
          event.preventDefault();
          showSystemInfo(dot);
        }
      });
    });

    if (controls) {
      controls.querySelectorAll("[data-action]").forEach(function (button) {
        button.addEventListener("click", function () {
          var action = button.dataset.action;
          if (action === "zoom-in") setZoom(zoom + ZOOM_STEP);
          else if (action === "zoom-out") setZoom(zoom - ZOOM_STEP);
          else if (action === "reset") resetView();
        });
      });
    }
  }

  var stageEl = document.getElementById("starmap-stage");
  if (stageEl) {
    initStarmap(stageEl);
  }
})();
