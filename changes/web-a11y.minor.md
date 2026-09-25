### Added
- **Browser checks for every Flask page.** A new test
  (`src/tests/test_web_a11y.py`, and its own `browser-a11y` CI job) loads
  each page in headless Chromium at phone and desktop widths, in light and
  dark, and fails on serious or critical axe-core (WCAG 2.1 AA)
  violations, horizontal page scroll, console or CSP errors, or a missing
  skip link or `aria-current` marker. The page list comes from the app's
  routes, so pages moved off CGI later are checked automatically.
  axe-core 4.13.0 is vendored under `src/tests/vendor/axe-core/`; the new
  `browser` extra installs Playwright.

### Fixed
- **The current section in the header was below WCAG AA contrast** in the
  light theme (4.45:1); it now uses the link colour.
