### Added

- **Shoelace components for buttons and menus (UX.40, first step).** The site now ships a subset of [Shoelace](https://shoelace.style) 2.20.1 as vendored ES modules under `static/vendor/shoelace/` (no CDN, no bundler, same-origin files under the existing Content-Security-Policy), registered on every page by `static/components.js` and coloured from the site's own light and dark tokens (`static/shoelace-theme.css`). `scripts/vendor_shoelace.py` rebuilds the vendored folder from the npm package.

### Changed

- **The header's Menu and settings gear are Shoelace dropdowns (UX.2).** They close on an outside click or Escape and give focus back to their button, replacing the `<details>` script in `theme.js`; each panel is as wide as its longest entry (at most 22 rem, never past the screen) instead of a fixed 14 rem minimum.
