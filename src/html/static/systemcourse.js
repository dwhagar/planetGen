// html/static/systemcourse.js
//
// NAV.22: a NAV course inside one system, as the 3D system view draws it. The server (`planetgen/web/nav_page.py` `system_course`) sends
// each leg of the course that lies in the system as a polyline in kilometers in the system's own frame: from the body to the heliopause
// (an "out" leg), from the heliopause to the body (an "into" leg), or between two bodies ("within"), bent clear of what it would pass
// too close to (NAV.27). This module turns those kilometers into the units the view draws in (systemscale.js), so the line meets the
// heliopause sphere and the bodies wherever they are drawn. Plain data in, plain data out: no three.js here.

// In the compressed scale the distance from the middle is on a log scale, so a straight line in kilometers is a curve in the view; it is
// cut into this many pieces for each of its segments so the curve is smooth. At true scale a straight line stays straight.
export const PIECES_PER_SEGMENT = 24;

// The points of `points` (kilometers) as the layout draws them: [[x, y, z], ...] in world units.
export function drawnLine(points, layout) {
  const pieces = layout.mode === "true" ? 1 : PIECES_PER_SEGMENT;
  const line = [];
  for (let n = 0; n + 1 < points.length; n += 1) {
    const a = points[n];
    const b = points[n + 1];
    for (let k = 0; k < pieces; k += 1) {
      const t = k / pieces;
      line.push(layout.place("barycenter", [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]));
    }
  }
  if (points.length) line.push(layout.place("barycenter", points[points.length - 1]));
  return line;
}

// What to draw for `course` (the server's `system_course`) under `layout`: one entry for each leg, `{kind, label, from, to, line, start,
// end, aimed}`; `start` and `end` are the line's two ends and `aimed` is where a moving destination is aimed at (or null).
export function drawnCourse(course, layout) {
  if (!course || !Array.isArray(course.legs)) return [];
  return course.legs.filter((leg) => Array.isArray(leg.points) && leg.points.length >= 2).map((leg) => {
    const line = drawnLine(leg.points, layout);
    return {
      kind: leg.kind, label: leg.label, from: leg.from, to: leg.to, line: line, start: line[0], end: line[line.length - 1],
      aimed: leg.aimed ? layout.place("barycenter", leg.aimed) : null,
    };
  });
}
