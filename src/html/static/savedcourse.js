// html/static/savedcourse.js
//
// A saved course (NAV.17, docs/TODO.md NAV.4): the record kept for one course, with both of its forms, and what
// opening it again says has changed. The record is an entry of the `course` kind of the saved-list store
// (savedlist.js); the Courses list and the Save button that fill and show it are NAV.18 and NAV.4.
//
// A record is made from what `GET /api/nav?from=...&to=...` answers:
// - from, to: the two end references (`system:12`, `nebula:3`; the printed forms the API answers with), and
//   fromName, toName, as they were named then;
// - name: what the user called it; saved: when (an ISO time);
// - direct: the straight line, {distance_ly, bearing_deg, mark_deg} (NAV.1);
// - route: the system-to-system route, {path, distance_ly, hops: [{distance_ly, bearing_deg, mark_deg}]}, or null
//   where there is none (both ends in one system, or no route through known space);
// - adjusted: the path bent clear of the gravity wells it passes (NAV.6), as the API answers it, or null until a
//   course has one.
// A course is known by its two ends, in order.
//
// `changesSince(record, result)` plans the course again from the current galaxy's answer and lists what differs, in
// words: the correlative update moves systems along their galactic orbits and a regenerated system can vanish.
// An empty list means the record still holds.

const VERSION_QUERY = new URL(import.meta.url).search;
const { createStore, localBackend } = await import(`./savedlist.js${VERSION_QUERY}`);
const { formatDistanceLy } = await import(`./unitladder.js${VERSION_QUERY}`);

export const COURSE_KIND = "course";
export const MAX_COURSES = 100;
const KEY_PREFIX = "planetgen.courses.";

// How far a number may move before it counts as changed: a thousandth of a light-year, a thousandth of a degree.
const DISTANCE_EPSILON_LY = 1.0e-3;
const ANGLE_EPSILON_DEG = 1.0e-3;

export function isCourse(entry) {
  return !!entry && typeof entry === "object" && typeof entry.from === "string" && typeof entry.to === "string"
    && typeof entry.name === "string" && entry.from !== "" && entry.to !== "";
}

export function courseKey(entry) {
  return isCourse(entry) ? entry.from + ">" + entry.to : null;
}

function round(value) {
  return typeof value === "number" && Number.isFinite(value) ? value : null;
}

function leg(course) {
  if (!course) return null;
  return { distance_ly: round(course.distance_ly), bearing_deg: round(course.bearing_deg), mark_deg: round(course.mark_deg) };
}

function routeOf(route) {
  if (!route) return null;
  return {
    path: (route.path || []).map(String),
    distance_ly: round(route.distance_ly),
    hops: (route.hops || []).map(leg),
  };
}

// The record for the course `result` (the answer of `/api/nav`) plans, named `name` (default: "<from> to <to>"),
// saved at `now` (a Date; default: the time now).
export function recordOf(result, name, now) {
  const from = result.origin.ref;
  const to = result.destination.ref;
  const fromName = result.origin.name || from;
  const toName = result.destination.name || to;
  return {
    kind: COURSE_KIND,
    name: String(name || "").trim() || fromName + " to " + toName,
    from: from,
    to: to,
    fromName: fromName,
    toName: toName,
    saved: (now || new Date()).toISOString(),
    direct: leg(result.direct),
    route: routeOf(result.route),
    adjusted: result.adjusted || null,
  };
}

function formatCourse(course) {
  const pad = function (degrees) { return String(Math.round(((degrees % 360) + 360) % 360) % 360).padStart(3, "0"); };
  return pad(course.bearing_deg) + " mark " + pad(course.mark_deg);
}

function differs(a, b, epsilon) {
  if (a === null || b === null) return a !== b;
  return Math.abs(a - b) > epsilon;
}

function legChanges(label, was, now) {
  if (!was || !now) return [];
  const changes = [];
  if (differs(was.distance_ly, now.distance_ly, DISTANCE_EPSILON_LY)) {
    changes.push(label + " is now " + formatDistanceLy(now.distance_ly) + ", it was " + formatDistanceLy(was.distance_ly) + ".");
  }
  if (differs(was.bearing_deg, now.bearing_deg, ANGLE_EPSILON_DEG) || differs(was.mark_deg, now.mark_deg, ANGLE_EPSILON_DEG)) {
    changes.push(label + " heads " + formatCourse(now) + ", it was " + formatCourse(was) + ".");
  }
  return changes;
}

// What differs between `record` and `result`, the answer of `/api/nav` for the same two ends now: a list of
// sentences, empty when nothing has changed.
export function changesSince(record, result) {
  const changes = [];
  const now = recordOf(result, record.name, new Date(0));
  changes.push(...legChanges("The direct course", record.direct, now.direct));
  if (record.route && !now.route) {
    changes.push("The route through known space is gone: no route is found now.");
  } else if (!record.route && now.route) {
    changes.push("A route through known space exists now: " + now.route.path.length + " stops, " + formatDistanceLy(now.route.distance_ly) + ".");
  } else if (record.route && now.route) {
    const was = record.route.path.join(" ");
    if (was !== now.route.path.join(" ")) {
      const stops = now.route.path.length === record.route.path.length
        ? "still " + now.route.path.length + " stops"
        : now.route.path.length + " stops; it was " + record.route.path.length;
      changes.push("The route now passes through different systems (" + stops + ").");
    }
    if (differs(record.route.distance_ly, now.route.distance_ly, DISTANCE_EPSILON_LY)) {
      changes.push("The route is now " + formatDistanceLy(now.route.distance_ly) + ", it was " + formatDistanceLy(record.route.distance_ly) + ".");
    }
  }
  if (record.adjusted !== null && now.adjusted !== null
      && JSON.stringify(record.adjusted) !== JSON.stringify(now.adjusted)) {
    changes.push("The adjusted path around gravity wells has changed.");
  }
  return changes;
}

// What opening a course says when the course can no longer be planned (the API's error for the same two ends).
export function unavailableChange(message) {
  return "The course can no longer be planned: " + (message || "an end is gone") + ".";
}

// The saved courses of the database named by `suffix()` (a function giving its name), in the browser.
export function courseStore(suffix) {
  return createStore({
    backend: localBackend({ prefix: KEY_PREFIX, suffix: suffix }),
    keyOf: courseKey,
    max: MAX_COURSES,
    words: {
      invalid: "That course can't be saved.",
      full: function (most) { return "You have " + most + " saved courses, the most there can be; delete one first."; },
      blocked: "This browser isn't letting the site save courses.",
    },
  });
}
