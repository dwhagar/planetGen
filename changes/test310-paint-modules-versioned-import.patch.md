**Headline:** The solid-world and star painters load their shared helper with the release stamp, so a new release cannot serve a stale copy.

### Changed

- `solidpaint.js` and `starpaint.js` import `bodytexture.js` with the page's `?v=` stamp, as every other script does (TEST.310; the page-shell test caught it).
