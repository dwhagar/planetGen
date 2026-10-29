**Headline:** Scrolling the wheel toward a dense crowd of stars now brings the camera into the crowd and stops there instead of flying through it to a star behind.

### Fixed

- **Wheel zoom heads for the most populated layer in front of the zoom direction (MAP.301).** The wheel zooms toward the thing under the pointer, and when that was a lone star behind a crowd (a cluster, a dense region) the camera passed through the crowd without stopping. The pick now also looks along the ray under the pointer for the densest layer of drawn stars in front of that thing (a bounded search over log-spaced depth slices, a narrow cone around the ray, at least twelve stars) and, where it finds one, the approach slows and ends in that layer; the next scroll carries on past it. A pointer over empty space, a lone star, a body, or a crowd smaller than twelve stars is unaffected. No stored value or generated output changes.
