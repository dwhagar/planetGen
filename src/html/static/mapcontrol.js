// html/static/mapcontrol.js
//
// The camera and input controller the Galaxy Map's drill-down
// (galaxystageview.js) and the Sector Map (sectormap.js) share (MAP.64):
// turning an orbit camera, moving its target in the screen's plane,
// zooming within the view's zoom policy, a two-finger pinch, the arrow
// keys, and telling a drag from a click. Each map keeps its own numbers
// (how far a press may wander and still be a click, how fast a drag
// turns) and its own reaction to a click; this module owns the
// mechanics, so both maps behave alike where they meant to and the rules
// live in one place. No visible change from the code it replaces.
//
// A zoom policy is what a view lets the wheel, a pinch and the zoom
// buttons do (MAP.58's rule is the locked one):
//   "free"   -- any distance;
//   "range"  -- a camera distance from `nearest` to `farthest` times the
//               view's own fit (a short range around it);
//   "locked" -- no zoom at all: the view stays at its fit and the wheel
//               is left to scroll the page.

export const ZOOM_FREE = "free";
export const ZOOM_RANGE = "range";
export const ZOOM_LOCKED = "locked";

// A zoom policy of `kind`; `nearest` and `farthest` are only read for a
// range.
export function zoomPolicy(kind, nearest, farthest) {
  return { kind: kind, nearest: nearest, farthest: farthest };
}

// Whether the policy lets the view zoom at all.
export function canZoom(policy) {
  return !!policy && policy.kind !== ZOOM_LOCKED;
}

// A camera distance held to the policy, given the view's fit distance
// (a locked view stays at its fit).
export function clampDistance(policy, distance, fitDistance) {
  if (!policy || policy.kind === ZOOM_FREE) return distance;
  if (policy.kind === ZOOM_LOCKED) return fitDistance;
  return Math.max(fitDistance * policy.nearest, Math.min(fitDistance * policy.farthest, distance));
}

// --- Orbit camera ------------------------------------------------------------------

// Turns an orbit view ({theta, phi}, radians, phi measured from the
// pole) by a drag of dx, dy pixels at `perPx` radians a pixel, the tilt
// held by `clampPhi`. Changes and returns `view`.
export function orbitByDrag(view, dx, dy, perPx, clampPhi) {
  view.theta -= dx * perPx;
  view.phi = clampPhi(view.phi - dy * perPx);
  return view;
}

// The arrow keys' turn of an orbit view: Left and Right swing it around
// by `step`, Up and Down tilt it (held by `clampPhi`). True when the key
// was one of them and the view changed.
export function orbitByKey(view, key, step, clampPhi) {
  if (key === "ArrowLeft") view.theta += step;
  else if (key === "ArrowRight") view.theta -= step;
  else if (key === "ArrowUp") view.phi = clampPhi(view.phi - step);
  else if (key === "ArrowDown") view.phi = clampPhi(view.phi + step);
  else return false;
  return true;
}

// Moves `target` ([x, y, z]) against a drag of dx, dy pixels in the
// screen's plane of `camera` (its matrixWorld up to date), `perPx` world
// units a pixel, so what is under the pointer follows it. Changes and
// returns `target`.
export function panInScreenPlane(THREE, camera, target, dx, dy, perPx) {
  const right = new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld, 0);
  const upward = new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld, 1);
  target[0] -= (right.x * dx - upward.x * dy) * perPx;
  target[1] -= (right.y * dx - upward.y * dy) * perPx;
  target[2] -= (right.z * dx - upward.z * dy) * perPx;
  return target;
}

// --- Pointer, wheel and pinch ------------------------------------------------------

// A wheel event's travel in pixels, lines and pages turned into pixels
// and the whole held to +-`limitPx`.
export function wheelPixels(event, pageHeightPx, limitPx) {
  let deltaPx = event.deltaY;
  if (event.deltaMode === 1) deltaPx *= 33;
  else if (event.deltaMode === 2) deltaPx *= pageHeightPx || 400;
  return Math.max(-limitPx, Math.min(limitPx, deltaPx));
}

// The pointer half of the controller, as handlers for the canvas's
// pointerdown, pointermove, pointerup (and pointercancel), pointerleave
// and click events; with `options.attach` it adds them to `canvasEl`
// itself. `options`:
//   dragClickPx     how far a press may travel and still be a click;
//   measure         "offset" (straight from where it went down, the
//                   default) or "path" (every move added up);
//   buttons         the mouse buttons that count (default any);
//   isPan(event)    whether a press moves the view rather than turning it
//                   (a pan press never clicks);
//   turnAtOnce      the view follows the pointer from its first move,
//                   before it has travelled dragClickPx;
//   canDrag()       whether a press past dragClickPx may start a drag
//                   (default yes);
//   onDragStart()   a drag has begun;
//   onDrag(dx, dy, pan)   the pointer moved while dragging;
//   pinch           {canStart(), start(), move(ratio)} for two fingers on
//                   a touch screen: ratio is the fingers' first spread
//                   over their spread now (more than 1 when they close);
//   clickOn         "pointerup" (the default: a press released close to
//                   where it went down, not with the right button) or
//                   "click" (the browser's click, skipped after a drag);
//   onClick(event, pointerType)   a click;
//   onHover(event)  the pointer moved with no press down;
//   onLeave(event)  it left the canvas.
export function createPointerControl(canvasEl, options) {
  const o = options || {};
  const clickPx = o.dragClickPx || 0;
  const byPath = o.measure === "path";
  let pointer = null;
  // Touch points down on the canvas, for a pinch: id -> {x, y}.
  const touches = new Map();
  let pinch = null;
  // After a drag, the browser's click that follows is not a click.
  let suppressNextClick = false;

  function spread() {
    const p = Array.from(touches.values());
    return Math.hypot(p[0].x - p[1].x, p[0].y - p[1].y) || 1;
  }

  function travelled(event) {
    return byPath ? pointer.path : Math.abs(event.clientX - pointer.x0) + Math.abs(event.clientY - pointer.y0);
  }

  function down(event) {
    if (o.buttons && o.buttons.indexOf(event.button) < 0) return;
    if (event.pointerType === "touch" && o.pinch) {
      touches.set(event.pointerId, { x: event.clientX, y: event.clientY });
      if (touches.size === 2 && o.pinch.canStart()) {
        pinch = { d: spread() };
        o.pinch.start();
        pointer = null;
        return;
      }
    }
    pointer = {
      x0: event.clientX, y0: event.clientY, x: event.clientX, y: event.clientY, path: 0,
      id: event.pointerId, type: event.pointerType,
      pan: !!(o.isPan && o.isPan(event)), dragging: false,
    };
    try {
      canvasEl.setPointerCapture(event.pointerId);
    } catch (err) {
      // Pointer capture isn't essential: dragging still works through
      // ordinary pointermove bubbling if the browser refuses it.
    }
  }

  function move(event) {
    if (event.pointerType === "touch" && touches.has(event.pointerId)) {
      touches.set(event.pointerId, { x: event.clientX, y: event.clientY });
      if (pinch && touches.size === 2) {
        o.pinch.move(pinch.d / spread());
        return;
      }
    }
    if (pointer && event.pointerId === pointer.id) {
      const dx = event.clientX - pointer.x;
      const dy = event.clientY - pointer.y;
      pointer.x = event.clientX;
      pointer.y = event.clientY;
      pointer.path += Math.abs(dx) + Math.abs(dy);
      if (!pointer.dragging && (o.turnAtOnce || (travelled(event) > clickPx && (!o.canDrag || o.canDrag())))) {
        pointer.dragging = true;
        if (o.onDragStart) o.onDragStart();
      }
      if (pointer.dragging && o.onDrag) {
        o.onDrag(dx, dy, pointer.pan || !!(o.isPan && event.shiftKey));
      }
      return;
    }
    if (o.onHover) o.onHover(event);
  }

  function up(event) {
    if (event.pointerType === "touch" && o.pinch) {
      touches.delete(event.pointerId);
      if (pinch) {
        if (touches.size === 0) pinch = null;
        pointer = null;
        return;
      }
    }
    if (!pointer || event.pointerId !== pointer.id) return;
    const moved = travelled(event);
    const type = pointer.type;
    // A turn-at-once press is "dragging" from its first move: whether it
    // was a drag is the distance it travelled.
    const wasDrag = (pointer.dragging && !o.turnAtOnce) || pointer.pan || moved > clickPx;
    pointer = null;
    try {
      canvasEl.releasePointerCapture(event.pointerId);
    } catch (err) {
      // Already released: nothing to clean up.
    }
    if (o.clickOn === "click") {
      if (wasDrag) {
        suppressNextClick = true;
        // A pointerup isn't always followed by a click (a pointercancel):
        // don't leave this eating some later, unrelated click.
        setTimeout(function () {
          suppressNextClick = false;
        }, 0);
      }
      return;
    }
    if (wasDrag || event.button === 2) return;
    if (o.onClick) o.onClick(event, type);
  }

  function click(event) {
    if (o.clickOn !== "click") return;
    if (suppressNextClick) {
      suppressNextClick = false;
      return;
    }
    if (o.onClick) o.onClick(event, "mouse");
  }

  function leave(event) {
    if (o.onLeave) o.onLeave(event);
  }

  const handlers = { down: down, move: move, up: up, leave: leave, click: click };
  if (o.attach) {
    canvasEl.addEventListener("pointerdown", down);
    canvasEl.addEventListener("pointermove", move);
    canvasEl.addEventListener("pointerup", up);
    canvasEl.addEventListener("pointercancel", up);
    canvasEl.addEventListener("pointerleave", leave);
    canvasEl.addEventListener("click", click);
  }
  return handlers;
}
