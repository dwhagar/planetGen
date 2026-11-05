**Headline:** The course planner can now bend a straight course around the objects it would pass too close to.

### Added
- `galaxy/avoid.py` takes a start, an end and the keep-out spheres between them and returns a short detour: a straight line when nothing is in the way, otherwise tangent lines and arcs around the spheres it crosses. A sphere that holds the start or the end is ignored, overlapping spheres are treated as one larger sphere, and a course it cannot bend safely is given as the straight line with a warning. Nothing shows it yet; the route and the NAV page call it next (NAV.6).

generation-output: unchanged
