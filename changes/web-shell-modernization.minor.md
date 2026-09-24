### Added
- **A light/dark/system theme switch** at the bottom of the side nav. The choice is
  remembered in the browser (`static/theme.js`) and applied before the page first
  paints, so there is no flash of the wrong theme. `style.css` gains the matching
  `:root[data-theme="dark"]` token block.
- **A favicon** (`static/favicon.svg`), a `<meta name="description">`, and
  `<meta name="theme-color">` for light and dark.

### Fixed
- **Phones got the desktop layout zoomed out.** The page shell now has
  `<meta name="viewport" content="width=device-width, initial-scale=1">`, so the
  existing narrow-screen rules apply. Those rules are also fixed: the side nav
  becomes a wrapping row of links instead of full-width buttons, the page gets an
  even gutter, and wide tables (including ones in wiki notes) scroll inside their
  own box instead of widening the page.

### Changed
- **Static files are versioned and cacheable.** Every CSS/JS/icon URL carries
  `?v=<release version>` (`fmt.static_url`), and the map scripts pass the same
  version on to three.js and `bodyRendering.js`. The Apache example now tells
  browsers to keep versioned `static/` files for a year, revalidates unversioned
  ones, and compresses HTML, CSS, JS, JSON and SVG with `mod_deflate` (three.js
  goes from about 740 KB to about 190 KB). Needs `a2enmod headers deflate`.

### Security
- **One source for the pages' security headers, and a stronger CSP.**
  `lib/page.py`'s `SECURITY_HEADERS` is the only place the HTML pages get them;
  the Apache example no longer repeats them vhost-wide (remove those three
  `Header always set` lines from an existing site config). The CSP adds
  `base-uri 'self'; form-action 'self'; frame-ancestors 'none'; object-src 'none'`
  to `default-src 'self'`.
