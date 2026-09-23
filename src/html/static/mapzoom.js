// html/static/mapzoom.js
//
// Shared viewBox-based zoom/pan for a flat, static, server-rendered
// `<svg>` map -- phenomenonmap.js (a stellar phenomenon's own AU-scale
// diagram) is this module's one caller today (the Galaxy Map used to be
// a second one, static/galaxymap.js, before it became a real 3D scene --
// see lib/galaxymap3d.py); written generically rather than folded into
// that one caller directly, so a future flat SVG map can reuse it the
// same way instead of hand-rolling the same wheel-zoom/drag-pan/
// click-vs-drag logic again.
// Unlike sectormap.js's/galaxymap3d.js's three.js scenes (a real 3D
// camera), there is no camera here at all -- "zooming" is just shrinking/
// growing the SVG's own `viewBox` rect, which the browser already
// re-renders at full vector fidelity (every marker/label was drawn in
// real SVG units server-side, so magnifying the viewBox magnifies them
// too, for free -- no separate level-of-detail logic needed). This is
// also exactly why a flat SVG map's own markers grow steadily *larger*
// relative to the view as you zoom in with no camera to shrink them the
// opposite way -- the effect a real 3D camera (sectormap.js/
// galaxymap3d.js) doesn't have, since a sprite's *world* size stays
// fixed while its on-screen size naturally falls off with distance.
//
// Exponential zoom (not a linear px step) throughout, since a caller like
// phenomenonmap.js's own AU-scale diagram spans a huge dynamic range (a
// full light-year down to 1 AU is a ~63,000x span) where a fixed linear
// step would be unusably coarse at the zoomed-in end or unusably slow at
// the zoomed-out end.

(function () {
  "use strict";

  var DRAG_THRESHOLD_PX = 4;
  var WHEEL_ZOOM_STEP = 1.15;
  var BUTTON_ZOOM_STEP = 1.4;

  function parseViewBox(svgEl) {
    var parts = (svgEl.getAttribute("viewBox") || "0 0 100 100").split(/\s+/).map(Number);
    return { x: parts[0], y: parts[1], w: parts[2], h: parts[3] };
  }

  function setViewBox(svgEl, box) {
    svgEl.setAttribute("viewBox", box.x + " " + box.y + " " + box.w + " " + box.h);
  }

  // Converts a client (mouse/pointer) coordinate to this SVG's own
  // user-space coordinate, via the inverse of its current screen
  // transform -- the standard, exact way to do this regardless of the
  // SVG's own current size/viewBox/CSS scaling (rather than trying to
  // hand-derive the same ratio from getBoundingClientRect()).
  function clientToSvgPoint(svgEl, clientX, clientY) {
    var ctm = svgEl.getScreenCTM();
    if (!ctm) {
      return { x: 0, y: 0 };
    }
    var point = svgEl.createSVGPoint();
    point.x = clientX;
    point.y = clientY;
    var svgPoint = point.matrixTransform(ctm.inverse());
    return { x: svgPoint.x, y: svgPoint.y };
  }

  // Zooms `box` by `factor` (>1 zooms in/shrinks the box, <1 zooms out)
  // around the fixed SVG-space point `anchor` -- the point under the
  // cursor (wheel zoom) or the box's own current center (button zoom)
  // stays visually still while everything else scales around it.
  function zoomedBox(box, factor, anchor, minSize, maxSize) {
    var newW = box.w / factor;
    // Clamp on width alone (since this module's own callers both use a
    // square viewBox) between minSize (the deepest allowed zoom-in) and
    // maxSize (the furthest allowed zoom-out), then rescale height to
    // match, keeping the box's own aspect ratio fixed.
    var clampedW = Math.min(maxSize, Math.max(minSize, newW));
    var scale = clampedW / box.w;
    var clampedH = box.h * scale;

    var ax = (anchor.x - box.x) / box.w;
    var ay = (anchor.y - box.y) / box.h;
    return {
      x: anchor.x - ax * clampedW,
      y: anchor.y - ay * clampedH,
      w: clampedW,
      h: clampedH,
    };
  }

  function boxCenter(box) {
    return { x: box.x + box.w / 2, y: box.y + box.h / 2 };
  }

  /**
   * Wires up wheel-zoom, drag-to-pan, and (if given) zoom-in/zoom-out/
   * reset buttons on `svgEl` -- a flat, static `<svg>` this module treats
   * as a plain 2D camera-less "canvas" (its own `viewBox` IS the camera).
   *
   * @param {SVGSVGElement} svgEl - The map's root `<svg>`. Read once for
   *   its starting `viewBox` (the "reset" target) and its own
   *   `data-min-view-size` (in SVG user units -- the narrowest a
   *   `viewBox` edge may become, i.e. the deepest allowed zoom-in) /
   *   `data-max-view-size` (the widest, i.e. the furthest allowed
   *   zoom-out -- defaults to the starting `viewBox`'s own width when
   *   absent, so a caller that only cares about a zoom-IN limit doesn't
   *   also have to spell out today's already-correct zoom-out bound)
   *   unless overridden by `options`.
   * @param {Object} [options]
   * @param {number} [options.minSize] - Overrides `data-min-view-size`.
   * @param {number} [options.maxSize] - Overrides `data-max-view-size`.
   * @param {HTMLElement} [options.zoomInBtn]
   * @param {HTMLElement} [options.zoomOutBtn]
   * @param {HTMLElement} [options.resetBtn]
   * @param {function(Object):void} [options.onChange] - Called with the
   *   new `{x, y, w, h}` box after every zoom/pan/reset, e.g. to redraw a
   *   scale-bar legend against the live zoom level.
   */
  function initSvgZoomPan(svgEl, options) {
    options = options || {};
    var initialBox = parseViewBox(svgEl);
    var minSize = options.minSize || parseFloat(svgEl.dataset.minViewSize) || initialBox.w;
    var maxSize = options.maxSize || parseFloat(svgEl.dataset.maxViewSize) || initialBox.w;

    function apply(box) {
      setViewBox(svgEl, box);
      if (options.onChange) {
        options.onChange(box);
      }
    }

    function zoomBy(factor, anchor) {
      apply(zoomedBox(parseViewBox(svgEl), factor, anchor, minSize, maxSize));
    }

    svgEl.addEventListener(
      "wheel",
      function (event) {
        event.preventDefault();
        var anchor = clientToSvgPoint(svgEl, event.clientX, event.clientY);
        zoomBy(event.deltaY < 0 ? WHEEL_ZOOM_STEP : 1 / WHEEL_ZOOM_STEP, anchor);
      },
      { passive: false }
    );

    if (options.zoomInBtn) {
      options.zoomInBtn.addEventListener("click", function () {
        zoomBy(BUTTON_ZOOM_STEP, boxCenter(parseViewBox(svgEl)));
      });
    }
    if (options.zoomOutBtn) {
      options.zoomOutBtn.addEventListener("click", function () {
        zoomBy(1 / BUTTON_ZOOM_STEP, boxCenter(parseViewBox(svgEl)));
      });
    }
    if (options.resetBtn) {
      options.resetBtn.addEventListener("click", function () {
        apply(initialBox);
      });
    }

    // Drag-to-pan, plus the click-vs-drag gate every marker's own
    // data-nav-target click (navform.js's document-level delegated
    // listener) needs: a genuine drag must not also fire a navigation
    // click on whatever marker the pointer happened to release over.
    //
    // setPointerCapture is deliberately NOT called on pointerdown: per
    // the Pointer Events spec, once a pointer is captured, its follow-up
    // pointerup *and the click event synthesized from it* both retarget
    // to the capturing element (svgEl) instead of whatever marker was
    // actually under the pointer -- which broke every single marker
    // click (not just drags), since navform.js's delegated listener
    // looks for data-nav-target on the click's own target/ancestors, and
    // svgEl itself never carries it. Capturing only once a real drag is
    // detected (in pointermove, below) keeps a plain click's target
    // exactly as the browser's normal hit-test would have set it.
    var dragState = null;
    svgEl.addEventListener("pointerdown", function (event) {
      if (event.button !== 0) {
        return;
      }
      dragState = {
        pointerId: event.pointerId,
        startClientX: event.clientX,
        startClientY: event.clientY,
        startBox: parseViewBox(svgEl),
        dragged: false,
      };
    });

    svgEl.addEventListener("pointermove", function (event) {
      if (!dragState || event.pointerId !== dragState.pointerId) {
        return;
      }
      var dxClient = event.clientX - dragState.startClientX;
      var dyClient = event.clientY - dragState.startClientY;
      if (!dragState.dragged && Math.hypot(dxClient, dyClient) < DRAG_THRESHOLD_PX) {
        return;
      }
      if (!dragState.dragged) {
        svgEl.setPointerCapture(event.pointerId);
      }
      dragState.dragged = true;

      // Client-pixel delta converted to SVG-user-unit delta via the
      // current screen scale (box width / on-screen width) -- panning
      // needs a *vector* delta, not clientToSvgPoint's absolute point
      // (which would also need to account for the drag itself having
      // already moved the box, double-counting the motion).
      var rect = svgEl.getBoundingClientRect();
      var box = dragState.startBox;
      var svgDx = (dxClient / rect.width) * box.w;
      var svgDy = (dyClient / rect.height) * box.h;
      apply({ x: box.x - svgDx, y: box.y - svgDy, w: box.w, h: box.h });
    });

    function endDrag(event) {
      if (!dragState || event.pointerId !== dragState.pointerId) {
        return;
      }
      if (dragState.dragged) {
        svgEl.addEventListener(
          "click",
          function (clickEvent) {
            clickEvent.stopPropagation();
          },
          { capture: true, once: true }
        );
      }
      dragState = null;
    }
    svgEl.addEventListener("pointerup", endDrag);
    svgEl.addEventListener("pointercancel", endDrag);

    return {
      reset: function () {
        apply(initialBox);
      },
      getBox: function () {
        return parseViewBox(svgEl);
      },
    };
  }

  window.planetgenInitSvgZoomPan = initSvgZoomPan;
})();
