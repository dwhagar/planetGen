# Map interface and front-end libraries

How the browser side of planetGen should look and behave next, and what the
vendored front-end libraries cost to keep: the health of Shoelace, three.js,
Xterm.js and TanStack, the form-field work, the pilot-console look and its
measured colour tokens, soft nebula edges, Select mode and picking, the mini
map, bookmarks and charted regions, job warnings and logs, and the admin
density control. It adds measured numbers to the plan in
[library-migration.md](library-migration.md) and to the map plans in
[galaxy-drilldown-navigation.md](galaxy-drilldown-navigation.md).

Informs: UX.3, UX.40, UX.41, UX.42, UX.43, UX.45, UX.46, UX.47, UX.48, UX.49,
UX.76, MAP.75, MAP.102, MAP.119, MAP.121, MAP.122, MAP.132, MAP.139, MAP.140,
MAP.141, MAP.142, MAP.143

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result or a raw file (URL in Sources), [C]
computed or measured with a script or a registry read on 2026-10-09, [R]
recalled and unconfirmed. The research environment could read search-result
text and raw files but not papers or the full vendor sites, so every [R] item
is listed under Evidence notes to check when paper access is allowed.

## Decisions already taken

- Boss (2026-10-03 05:38Z): move to third-party libraries; Shoelace, TanStack,
  Xterm.js and three.js are vendored as ES module builds served by Flask, with no
  bundler and no CDN at runtime ([library-migration.md](library-migration.md),
  section 5).
- Boss (2026-10-07 11:47Z), UX.43: "All design elements should be built from the
  ground up to look like, as much as possible, a starmap and navigational system
  that a space pilot might use."
- Boss (2026-10-07 11:47Z), MAP.122: a "Select Mode" at all galactic views,
  defaulting to Galaxy (wedges, blocks, slabs), with a Star mode where stars are
  selectable and can be bookmarked or made a waypoint.
- Boss (2026-10-07 11:47Z), UX.46 and UX.48: wiping a galaxy wipes all bookmarks;
  each cluster of contiguous generated content gets a bookmark that zooms to show
  all its sectors and highlights the region's shape.
- Boss (2026-10-03 05:38Z), MAP.121: sectors between the camera and the sector in
  view fade, "the more they obstruct the view the more invisible they become";
  neighbours are dimmed without comets or rogue planets. Boss (2026-10-07 16:27Z):
  one map engine for sectors, galaxy and systems.
- Boss (2026-10-03 05:38Z), MAP.119: expected star density is editable by admins,
  static for users and guests, and an edit applies to the current view (a slab or
  block). Boss (2026-10-08 01:59Z), MAP.132: blocks are never fully opaque.
- Boss (2026-10-09 00:25Z to 00:44Z), MAP.140 to MAP.143: double-click goes to the
  selected object; faint context around the selection; fuzzy, gaseous nebula
  edges; colouring by habitable locations.
- Boss (2026-10-01), UX.3: a banner to every visitor while a background job
  changes the galaxy, with an ETA rounded up to the next hour.

## 1. Front-end library health

### 1.1 What is vendored

| Library | Vendored | Latest on npm (published) | Status | Note |
|---|---|---|---|---|
| @shoelace-style/shoelace | 2.20.1 | 2.20.1 (2025-03-11) | README: "Shoelace is sunset. There is no active development" [S]; no registry deprecation [C] | Final release; MIT; `npm audit` on the pinned set (Shoelace, three, xterm, TanStack) found 0 vulnerabilities [C] |
| three | 0.186.0 (r186) | 0.186.1 (2026-09-24) | none | About one minor every 2 months (0.183 Feb 2026, 0.184 Apr, 0.185 Jun, 0.186 Sep) [C] |
| @xterm/xterm, @xterm/addon-fit | 6.0.0, 0.11.0 | 6.0.0 (2025-12-22), 0.11.0; beta 6.1.0-beta.304 (2026-08-30) | old `xterm` and `xterm-addon-fit` packages are deprecated [C] | Vendored names are already the current ones |
| @tanstack/table-core | 8.21.3 | 9.2.8 (2026-10-08) | v8 is no longer `latest` [C] | Stay on 8 until a feature needs 9 |
| @tanstack/virtual-core | 3.13.12 | 3.17.11 (2026-09-14) | none | Four minors behind; low risk |
| three-mesh-bvh | not vendored, not imported anywhere in `src/html/static` | 0.9.16 (2026-10-08) | none | Dropped for picking, section 5.2 |

Source for every row: `registry.npmjs.org/<package>` JSON [C]. The Shoelace
README also says issues belong on Web Awesome and that the package stays
available under MIT; the Font Awesome blog posts agree that 2.20.1 is the last
release [S]. The Shoelace 2.20.1 files bundle Lit 3 (lit-html 3.2.1, lit-element
4.1.1, reactive-element 2.0.4 in the vendored chunks) plus floating-ui and
tinycolor [C]. Their licences are [R] (MIT or BSD-3); `vendor/shoelace/LICENSE.md`
covers Shoelace only, so the notices file should list the bundled ones.

### 1.2 Shoelace is finished but nothing breaks

Do not migrate now. The successor, Web Awesome (`@awesome.me/webawesome`
3.14.0, 2026-09-24), is active and MIT in its free tier (copyright Fonticons;
every version since 3.0.0-beta.1 is MIT in the registry [S][C]; 3.0.0 final
2025-10-28, then 3.1 to 3.14 every 2 to 4 weeks). It has 73 free components
including input, select, checkbox, switch, textarea, dialog, drawer, dropdown,
tab-group, tree, toast and progress-bar, a `shoelace` theme file to ease
migration, and a `native.css` for plain elements [C]. The paid Pro tier adds
components such as a date picker and data grid; the site needs none, and
current Pro prices were not found. Costs found by reading the 3.14.0 tarball:

- Size: `dist` is 13 MB unpacked against 8 MB for Shoelace [C]; only imported
  components ship, as now.
- Icons break the CSP. The default icon library fetches from
  `https://ka-f.fontawesome.com/releases/v7.3.0/svgs/...` [C], and the "system"
  icons are `data:` URIs. The site sends `default-src 'self'`, so both must be
  saved as local files, which is what `shoelaceicons.js` and
  `scripts/vendor_shoelace.py` already do for Shoelace.
- Themes use CSS relative colour syntax (`oklch(from currentColor ...)`) [C],
  which needs a recent browser (about Chrome 119, Safari 16.4, Firefox 128 [R]).
  The site's browser baseline is not written down.
- Tags become `wa-*`, "primary" becomes "brand", and the token set is new
  (`--wa-*`) [R for exact renames; the migration guide could not be fetched].
  `shoelace-theme.css` derives every `--sl-color-*` from the site tokens and
  would be rewritten, not just renamed.

Size of a migration today: 39 `<sl-*>` tags in 9 files under `src/planetgen`
(14 sl-button, 7 sl-menu-item, 7 sl-dialog, 5 sl-dropdown, 4 sl-menu, 1
sl-icon-button, 1 sl-divider) plus 2 in JavaScript [C]. UX.49 as written adds
about 14 selects, 5 textareas and 20-odd visible inputs. The migration is cheap
now and gets dearer with each added form field.

Alternatives looked at [C unless noted]: Material Web 2.5.0 (README: "maintenance
mode pending new maintainers" [S]), Lion `@lion/ui` 0.21.1 (unstyled, small
community [R]), Spectrum Web Components 1.12.4 and Carbon web components 2.65.0
(heavy, strong house looks), Pico CSS 2.1.1 (native elements only, no menus or
dialogs), Open Props 1.7.23 (variables only, useful as a token source). None
beats staying on Shoelace or moving to Web Awesome; native elements styled by
the site tokens remain the best choice for accessibility and no-JS behaviour.

Recommendation: keep Shoelace for the components in use; re-evaluate Web
Awesome when Shoelace develops a real problem (browser break, advisory) or when
the UX.43 restyle makes Shoelace's token model the obstacle. If the move is
wanted at all, do it inside UX.43: one restyle, one rename.

### 1.3 three.js

- **Imports.** `three.module.js` imports `three.core.js` since r171 [S]. The
  repo bundles both into one file with esbuild and loads it by relative path
  with `?v=`, with no import map. That is right here: an import map is inline
  JSON, which `default-src 'self'` blocks (the CSP consequence is inferred; the
  inline-JSON rule is [S]). Keep relative imports.
- **Reproducible [C].** `npm pack three@0.186.0`, then
  `esbuild --bundle --minify --format=esm package/build/three.module.js` with
  esbuild 0.28.2 gives SHA-256
  `fba62dfd6f0cdf58c9647d9cfe6e83cec5cf98eb24a7677ad116f4c0fca33902`, identical
  to `static/vendor/three.module.min.js` (742,188 bytes). Stay on r186; 0.186.1
  is a patch.
- **Tree shaking is not worth a build step [C].** The code names 41 `THREE.*`
  symbols. A bundle of only those is 562,098 bytes (140 KB gzip) against 742,188
  (189 KB gzip), because `WebGLRenderer` dominates. Gzip or brotli on the front
  end matters more.
- **Deprecations in r186 [C].** `Clock` (since r183; "This module has been
  deprecated. Please use THREE.Timer instead." in the vendored file),
  `Source` renamed `TextureSource` (r186), `Matrix3.scale/rotate/translate`
  (r185), `ColorManagement.fromWorkingColorSpace/toWorkingColorSpace` (r177),
  `PCFSoftShadowMap` (unused). Only `Clock` is used, at
  `phenomenonrender.js:528`. Replace it with `THREE.Timer` before the next bump.
- **WebGPU.** The manual describes `WebGPURenderer` as the modern renderer with
  automatic WebGL 2 fallback, and says larger features are no longer added to
  `WebGLRenderer` [S]. The engine is raw GLSL `ShaderMaterial` code (stars,
  cloud volumes, glow) that would need a rewrite to TSL. No benefit justifies
  that for 70,000 points and a few volumes. Revisit only if a target browser
  drops WebGL.
- **Notes for new shaders.** The Galaxy Map renderer uses
  `logarithmicDepthBuffer: true` (`galaxymap3d.js:592`) and caps pixel ratio at
  2 (line 593). Every custom shader needs the `logdepthbuf_*` chunks as the star
  shader has, and log depth disables early-z, which matters for overdraw.
  `Data3DTexture` and `readRenderTargetPixelsAsync` exist in r186 [C].

### 1.4 Vendored-version lock file (recommended)

Facts [C]: `three` has registry signatures but no provenance attestation (as do
Shoelace, xterm and addon-fit); the two TanStack packages carry SLSA
provenance. `vendor/THIRD_PARTY_NOTICES.txt` records versions and rebuild
commands but no hash, no esbuild version and no check that the files on disk
match.

Add `static/vendor/VENDORED.json`, written by `scripts/vendor_shoelace.py` and a
new small `scripts/vendor_other.py` (three, xterm, tanstack): per package, the
version, npm `dist.integrity`, `built_with` (for example `esbuild 0.28.2`),
licence and the SHA-256 of each shipped file. A pytest hashes every file under
`vendor/` and fails on a mismatch or an unlisted file; a vendor bump regenerates
the lock; the notices file names the esbuild version. Subresource Integrity adds
nothing for same-origin files; the lock catches accidental edits and shows a
reviewer exactly what changed. Optional CI step: `npm pack <pkg>@<version>`,
compare `dist.integrity`, rebuild, compare hashes (done by hand for three), and
run `npm audit` on a throwaway lock built from `VENDORED.json` (0 findings
today).

## 2. Form fields (UX.49)

Size is not the issue. Added to `COMPONENTS` [C, same import-closure rule as
`scripts/vendor_shoelace.py`]:

| Added | Extra files | Extra raw | Extra gzip |
|---|---|---|---|
| input | 6 | 28.5 KB | 6.9 KB |
| select and option | 12 | 52.0 KB | 13.0 KB |
| checkbox | 6 | 16.3 KB | 4.6 KB |
| textarea | 6 | 20.8 KB | 5.5 KB |
| radio-group and radio | 11 | 25.2 KB | 7.4 KB |
| input, select, option, checkbox, textarea together | 25 | 108.0 KB | 26.9 KB |

The current set is 63 files, 244 KB raw, 69 KB gzip. The risks are behavioural:

1. **No-JS and first paint.** An `sl-input` that is undefined, or with scripts
   blocked, does not exist, and `shoelace-theme.css` already hides undefined
   components. The site is server-rendered and forms post normally; a vanished
   field is worse than a plain input.
2. **Autofill and password managers** work less reliably inside a shadow root
   [R]. Sign-in, account, 2FA (`autocomplete="one-time-code"`) and new-password
   fields should stay native.
3. **Server-side validation errors** (ADM.21's Pydantic errors list every field)
   must show against the right field. `sl-input` has `setCustomValidity` and
   `help-text`, but the Jinja `field()` macros (`generate.html`,
   `generate_system.html`) need a second code path.
4. **Border contrast.** Shoelace's input border is `--sl-color-neutral-300`;
   with the site's token mapping that is 1.78:1 against the page in the light
   theme and 2.70:1 in dark [C]. WCAG 1.4.11 needs 3:1. Set
   `--sl-input-border-color` to the new `--line` token (section 3.3), 3.4:1 or
   better.
5. **Other Shoelace states against the site tokens [C].** Default button hover
   5.5 (light) and 5.9 (dark); primary button 5.08 and 5.84, its hover 4.23 and
   5.02; placeholder 4.41 and 5.76; menu-item focus 5.08 and 5.84; icon button
   5.71 and 8.04. Two light-theme values fail 4.5: primary hover and the
   placeholder. Fix by lowering the `primary-500` mix or darkening the
   placeholder token. These were measured before UX.76 landed (section 3.5) and
   should be rechecked.
6. **Tests.** The browser test UX.49 asks for fills each kind and submits;
   Playwright pierces open shadow roots [R].

Recommended re-scope (for Boss to confirm, since his wording is "every text
box"): one token-styled look for all fields; `sl-select` and
`sl-checkbox`/`sl-switch` where their features (clearable, listbox) earn it;
hidden, submit, password and one-time-code fields stay native. If UX.49 is done
on Shoelace, accept a later `sl-` to `wa-` rename; if as native controls, the
question is deferred.

## 3. The pilot's starmap console (UX.43, UX.42, UX.76)

### 3.1 References and what to take

- **Elite Dangerous** [S]: one dominant HUD hue (orange) with blue highlights, a
  colour matrix to recolour; avoid hue-only status (players with colour
  deficiency report hostile, neutral and friendly look alike).
- **Glass-cockpit alerting** [S, FAA AC 25-11, EASA]: red is warning (act now),
  amber is caution, advisory is any colour except red or green; white for general
  data. Magenta, cyan and green symbology meanings are [R].
- **NASA Open MCT** [S, raw SCSS, Apache-2.0]: dark neutral chrome, one key colour
  `#0099cc`, status colours kept for status (info `#60ba7b`, alert `#ffb66c`,
  error `#da0004`), 12 px type, 2 to 4 px radii. Avoid its low-contrast greys
  (`#aaa` on `#393939` is 4.4:1 [C]).
- **Game and desktop planetarium maps** (Star Citizen, Starfield, Stellaris,
  Celestia, SpaceEngine, Gaia Sky, Stellarium Web) [R, no sources retrieved]:
  dark field, thin line art, glowing selected markers, haloed labels, scale
  readouts; avoid ornament that costs legibility and endless animation.
- **LCARS**: nothing literal; the name, colour set, elbow bars and font are
  protected trade dress [R, not legal advice].
- **B612 Mono**: designed for cockpit displays, SIL OFL 1.1, `@fontsource/b612-mono`
  5.3.0, Latin woff2 19 KB per weight [C]; the Airbus and ENAC origin is [R].

Adopt the aviation meanings (red critical, amber caution or selected, green
nominal, cyan interactive, magenta planned route, white general data) and always
pair a colour with a symbol or word.

### 3.2 Design rules

- **Layout.** A chart-table instrument panel: a dark "space" field (the map canvas,
  already `#05070c` or `#05070d` in both themes) inside chrome of 1 px hairline
  borders, a corner tick or chamfer (not rounded 10 px cards) and a small-caps
  header strip. Keep the functional layout; add a few shared components (panel
  header, readout row, status chip, tick scale).
- **Type.** `system-ui` for sentences; readouts, coordinates, designations and
  table numerals in a monospace with `font-variant-numeric: tabular-nums
  slashed-zero` (candidate B612 Mono self-hosted, `ui-monospace` fallback). Labels
  11 to 12 px uppercase, 0.06 to 0.1em tracking, under about 24 characters; form
  inputs 16 px so iOS does not zoom.
- **Lines and glow.** 1 px hairlines for panels and grids; 1.5 to 2 px for selection
  and route; a 2 px focus ring with 2 px offset. Decorative grids may be low
  contrast (WCAG 1.4.11 exempts decoration); any line that conveys state must reach
  3:1. Glow (`drop-shadow(0 0 6px color-mix(in srgb, var(--accent) 55%,
  transparent))`) only on interactive state, never on text, and contrast is judged
  without it.
- **Motion.** Transitions 120 to 200 ms; pulses only under
  `prefers-reduced-motion: no-preference`, period at least 1.5 s, stopping after a
  few cycles (WCAG 2.3.1). `nebulaview.js` already gates its auto-turn; make that
  the shared pattern.
- **Icons and density.** Keep the 24x24 stroke-2 `currentColor` sprite (`icons.svg`,
  UX.28), `aria-hidden` with a text label on the button. Label:value pairs, numbers
  right-aligned, units muted at 85% size, at most 7 +/- 2 readouts per panel before
  a Details fold (MAP.139), one primary action per panel.
- **Wording (UX.42).** The same panel header carries the in-universe names
  ("Chart", "Charted only"); schedule UX.42 with UX.43 so each page is touched once.
  Admin pages may stay plain.

### 3.3 Token table (measured)

The site's current `style.css` tokens are `--bg`, `--bg-elevated`, `--bg-subtle`
(dark: `#23253470`, with alpha), `--text`, `--text-muted`, `--border`,
`--accent`, `--accent-contrast`, `--link`, `--link-visited` and `--error`
(current dark values: `--bg #14151e`, `--accent #8891f5`). The table keeps those
names so `shoelace-theme.css` follows automatically, and adds the console
tokens. `--line` is new and stricter than `--border`: `--border` stays for
decorative dividers, `--line` is for control borders and state lines. Contrast
is the WCAG 2.x relative-luminance ratio [C]; "checks" give text on bg, panel,
raised.

| Token | Dark | Light ("daylight chart") | Role | Dark checks | Light checks |
|---|---|---|---|---|---|
| `--bg` | `#05080d` | `#e9eef3` | page, map space | n/a | n/a |
| `--bg-elevated` | `#0b121b` | `#f6f8fa` | panels | n/a | n/a |
| `--bg-raised` | `#121c28` | `#ffffff` | menus, dialogs, inputs | n/a | n/a |
| `--text` | `#d9e6f2` | `#12202c` | body | 15.8 / 14.8 / 13.5 | 14.2 / 15.6 / 16.6 |
| `--text-muted` | `#93a9bd` | `#44596b` | secondary | 8.3 / 7.8 / 7.1 | 6.2 / 6.8 / 7.3 |
| `--text-faint` | `#7a93a9` | `#4e6377` | units, captions | 6.3 / 5.9 / 5.4 | 5.3 / 5.8 / 6.2 |
| `--line` | `#4b6b88` | `#6f8599` | control borders, state lines | 3.37 on panel | 3.59 on panel |
| `--grid` | `#16222f` | `#d3dde6` | decorative grid | 1.25 (decor) | 1.18 (decor) |
| `--accent` / `--link` | `#4dd0e6` | `#00667d` | interactive, hover ring | 11.0 / 10.3 / 9.4 | 5.6 / 6.2 / 6.6 |
| `--accent-contrast` | `#04121a` | `#ffffff` | text on accent fill | 10.4 on accent | 6.6 on accent |
| `--link-visited` | `#bfa8ff` | `#6a3fa5` | visited links | 9.8 / 9.2 | 6.3 / 6.9 |
| `--sel` | `#ffb347` | `#8a4b00` | selection, caution | 11.3 / 10.6 / 9.7 | 5.8 / 6.4 / 6.8 |
| `--route` | `#ff7fd0` | `#9d2b84` | planned course | 8.8 / 8.2 / 7.5 | 5.8 / 6.3 / 6.8 |
| `--ok` | `#5fd68f` | `#16703f` | nominal, charted | 11.0 / 10.3 / 9.4 | 5.3 / 5.8 / 6.1 |
| `--error` | `#ff6b6b` | `#b42318` | critical | 7.2 / 6.8 / 6.2 | 5.6 / 6.2 / 6.6 |
| `--focus` | `#ffd166` | `#7a3e00` | focus ring (2 px, 2 px offset) | 13.0 on panel | 7.8 on panel |

Every text pair is at least 4.5:1 and every control border at least 3:1 in both
themes [C]. Further measurements [C]:

- Hover fill `color-mix(accent 14%, panel)` keeps body text at 11.3 (dark) and
  12.7 (light). A pressed fill at 24% keeps text at 8.8 and 10.8, but
  accent-coloured text on it is only 4.27 in light, so pressed states use
  `--text`, not `--accent`.
- On the always-dark map background `#05080d`, in both themes: selection 11.3,
  route 8.8, hover cyan 11.0, ok 11.0, critical 7.2, label text 15.8, muted
  label 8.3.
- A focus ring over a cyan-filled button is 1.27, which is why the 2 px offset
  gap is required.
- Colour-blind check of the semantic set (CAM02-UCS delta E; below about 10 is
  hard to tell apart [R]): minimum pair differences 21 normal, 11.5 protan (cyan
  against text), 6.4 deutan (cyan against magenta), 4.7 tritan (cyan against
  green); red against green under deutan is 9.3. So the route is also dashed or
  thicker; critical is red plus an icon and a word; selection is a ring; the
  legend lists names, not only swatches. Amber is shared by selection and
  caution; they never appear in the same place (a ring on the map, a banner with
  icon and text elsewhere), as in aviation's single amber.

### 3.4 WCAG 2.2 AA checklist applied

[S] for the criterion list (Deque, BOIA: 2.4.11, 2.5.7, 2.5.8 and 3.3.8 are AA,
4.1.1 was removed); exact wording [R].

- 1.4.3 and 1.4.11: every text pair above 4.5, control borders (`--line`) at 3.4 or
  better, icons in `currentColor` of text tokens.
- 1.4.1 not colour alone: dashed route, icon plus word for critical, selection
  ring, named legend.
- 1.4.12 text spacing: no fixed-height label boxes. 2.3.1: pulses at 1.5 s or
  slower; reduced motion gates all decorative motion.
- 2.4.11 focus not obscured: sticky map controls must not cover focused fields.
- 2.5.7 dragging: move up and down buttons for bookmark order; a number field
  beside any slider. 2.5.8 target size: the 44 px coarse-pointer rule stays.
- 3.3.8 accessible authentication: unchanged.

### 3.5 UX.76 status

UX.76 (icons and highlights follow the theme) is done in PR #800. The cause was
a CSS cascade fault, not an icon-fill or token fault: `a.btn:visited` is as
specific as the secondary-button rules, so a visited secondary link kept the
primary button's dark text on a clear background. The fix pins
`a.btn.btn-secondary` (and its `:visited`) to `var(--accent)`, the active
variant to `var(--text)`, and makes the secondary hover an inset accent ring.
`src/tests/test_web_browser_colors.py` renders one of each button look under all
four combinations of OS scheme and `data-theme`, forces hover and focus through
the DevTools protocol, and checks computed contrast (text 4.5:1, icons 3:1).
The research suspects about the alpha-bearing dark `--bg-subtle`, runtime
elements with no background and a missing `currentColor` were not needed. Three
items remain worth carrying into UX.43: the vendored-SVG check (every file under
`vendor/shoelace/assets` and `icons.svg` contains `currentColor`) is not in that
test; the two light-theme Shoelace failures in section 2 are not exercised by
it; and `--bg-subtle` still carries alpha, so a background that depends on what
lies behind it can recur.

## 4. Nebulae with soft boundaries (MAP.142)

### 4.1 How nebulae are drawn today

- **Galaxy Map (`galaxymap3d.js`).** Each nebula or remnant is a `THREE.Sprite` with
  a 128 px radial-gradient texture (`cloudTexture`; diffuse alpha 0.47 at the centre,
  0.14 at 0.75 R, 0 at R). At 20 px radius or more (`NEBULA_MESH_MIN_PX`) it fetches
  `GET /galaxy/nebula/<id>/shape?lod=low` and draws the marching-cubes mesh instead
  (`buildNebulaMesh`: `MeshBasicMaterial`, `DoubleSide`, `depthWrite: false`,
  `depthTest: false`, opacity = edge opacity x 0.55). That flat shell is the hard edge
  Boss reported.
- **Sector Map (`sectorscene.js`).** A `BackSide` sphere with `alpha = (1 -
  exp(-density * path / (2R))) * maxAlpha`, `path` being the chord through a uniform
  sphere: soft toward the middle, a sharp limb, blind to the irregular shape.
- **Nebula page (`nebulaview.js`).** The translucent shape mesh with dimmed stars.
- **Colours.** `NEBULA_LOOKS` (core/edge opacity): diffuse `#e3a6c8` 0.47/0.14, emission
  `#ff6f91` 0.69/0.25, reflection `#6fa8ff` 0.63/0.22, planetary `#5be8c9` 0.66/0.24,
  dark `#1c1c24` 0.91/0.56, shared with `starmap.py`.
- **Shape.** `planetgen/galaxy/nebula_shape.py`: metaballs with kernel `(1 - d^2/r^2)^2`,
  domain-warped, marching cubes (`MESH_LODS` low 10, full 28 cells per axis). The scalar
  field `F` exists server side but is not served.

### 4.2 Techniques compared

- Instanced soft-noise sprites or a point cloud with radial falloff: good edges, fill-rate
  bound, need a scatter inside the shape. Celestia's nebulae work this way; early add-ons
  suffered huge overdraw, and additive blending cannot make dark dust [S].
- Fullscreen 3D-noise ray-march: excellent, but 128 steps x 7 evaluations per pixel is the
  literature's quality target [S]; Gaia Sky does it on desktop [S]. Not for this site.
- Ray-march a bounded proxy (front faces of the bounding mesh) through a small 3D texture
  of the real field: excellent, covered pixels only, 16 to 24 taps. Recommended for the
  focused nebula.
- Nested translucent iso-shells (iso x 1.0, 0.7, 0.5) with an edge-fade shader: good,
  about 6 fragments per pixel worst case, reuses marching cubes and the mesh loader.
- Today's gradient sprite: perfectly soft, wrong shape; keep it below about 20 px.

Cost [C]: pixel ratio is capped at 2, so the worst case is 8.29 Mpx (4K at ratio 1, or
1080p at ratio 2) and 1.32 Mpx for a 390x844 phone. A 16-step ray-march over every pixel
is 133 M samples per frame on 4K (8.0 G/s at 60 fps) and 21 M on the phone (1.3 G/s). No
GPU was available, so "heavy on integrated and phone GPUs" is [R]. Ray-march only the
covered pixels of the proxy, cap at 24 steps, stop at alpha 0.98, and render at half
resolution when one nebula covers over half the canvas (quarter resolution shows [S]).

### 4.3 Keeping the extent readable

Projected opacity for a radial profile is `1 - exp(-k N(b))`, `N(b)` the column density at
impact parameter `b`. Central alpha 0.45, `b` in units of the bounding radius R [C]:

| Profile | b at 50% of peak | b at 10% | Visible edge (alpha 5%) |
|---|---|---|---|
| Uniform sphere (Sector Map today) | 0.905 | 0.998 | 0.997 (sharp limb) |
| smoothstep 0.6R to R | 0.706 | 0.875 | 0.869 |
| smoothstep 0.3R to R | 0.588 | 0.821 | 0.812 |
| Metaball kernel (1 - r^2)^2 | 0.538 | 0.801 | 0.791 |

- A 0.6R to R falloff shrinks the visible edge to 0.87R, so a cloud whose sphere is its
  stated `radius_ly` looks about 13% smaller. Define the extent as the 50% contour
  instead, with a soft band symmetric about it (50% at R, 10% at about 1.1 R, 90% at
  about 0.85 R), so today's hard outline is the 50% contour and soft and hard views agree.
- The metaball kernel has a closed-form line integral, `N(b) = (16/15)(1 - b^2)^2.5`,
  checked against quadrature [C]: one `pow` in GLSL, no loop, no lookup. For an irregular
  shape use the real field.
- Serve a 24^3 or 32^3 grid of `F` as 8-bit (13.8 KB or 32.8 KB) from a new endpoint beside
  `/shape` (`NebulaShape.field` already evaluates it; derived data, so reproducibility is
  unaffected), as a `THREE.Data3DTexture`. Density: `d = smoothstep(iso*0.55, iso*1.45, F)`.
- Fade shells by view, alpha x `smoothstep(0.0, 0.35, abs(dot(n, v)))`, so the grazing
  silhouette dies away; the mesh has no normals yet (`computeVertexNormals` on load, or
  `dFdx`/`dFdy`).
- Keep the near fade (`CLOUD_NEAR_FADE [1.2, 2.0]` R) and the Sector scene's `BackSide`
  handling of the inside view.

### 4.4 Colours per class, colour-blind safety, dark nebulae

Composited core colour against the map background `#05070c` [C]: diffuse 2.92 (under
3:1), emission 4.02, reflection 3.82, planetary 6.03, remnant 3.40, dark 1.17
(effectively invisible on an empty map). Pairwise colour difference (CAM02-UCS, 100%
severity): lowest 16.1 in normal vision, 6.8 under protanomaly (diffuse and emission), 8.7
and 10.8 under deuteranomaly, 8.3 under tritanomaly [C]. An Okabe-Ito-based set (diffuse
`#CC79A7`, emission `#E8602C`, reflection `#56B4E9`, planetary `#2FD9A8`, remnant `#F0E442`)
lifts the normal minimum to 24.1 and protanomaly to 17.1 but leaves deuteranomaly at 8.8 and
tritanomaly at 7.5 [C]. Hue alone cannot carry the class; add a cue that survives greyscale:

- Emission: warm pink-red with faint filament texture. Reflection: blue, smoother, lower
  contrast. Planetary: a bright thin shell, dim interior, central point. Remnant: the thin
  shell with an orange or cyan rim. Diffuse: pale mauve, the faintest, with a 3:1 floor
  where selectable.
- Dark (absorbing): not a near-black translucent shell. Draw an occluder in the density pass
  (multiplicative blending `dst * (1 - a)`, darkening background nebulae and stars behind)
  plus a thin dashed `--text-muted` outline so it can be found. It is the only class with a
  different blend mode.
- Emission is additive in galaxy view (no sorting); other classes use normal alpha;
  `depthWrite: false` throughout. Selection and hover use a ring, not a colour change.

### 4.5 Recommendation

1. Galaxy Map: keep the sprite until 20 px; above that draw three nested iso-shells with the
   edge-fade shader (the server serves `lod=low` at iso 1.0, 0.7 and 0.5, about three times
   the current mesh bytes), fading in over 20 to 40 px so nothing pops.
2. Sector Map and nebula page: for the focused nebula, a bounded ray-march of the
   `Data3DTexture`, 16 steps, with the soft band; neighbours stay sprites or shells. Keep the
   analytic sphere shader as the first-frame and no-3D-texture fallback.
3. MAP.132's markers cover the far zoom; sprites under 2 px stay hidden (`CLOUD_MIN_PX`).
4. Tests: a node test of the closed-form profile (the 50% radius is exactly R for a symmetric
   band) and a browser test that screenshot alpha at 0.85 R, R and 1.2 R is above, near and
   below 50% of peak.

## 5. Interaction

### 5.1 Select mode (MAP.122)

- **Control.** A two-option segmented control, "Galaxy | Stars", as
  `role="radiogroup"` (Shoelace `sl-radio-group` with `sl-radio-button`, or two
  native radios styled as a segment). A toggle button hides which mode is on; a
  radiogroup shows both. Default Galaxy (Boss).
- **State.** In the URL (`select=stars`), so a view is shareable and back and
  forward restore it; "every galaxy view" then reads one parameter. A
  single-key shortcut, if added, fires only while the map has focus, as keys 1
  to 9 do now (WCAG 2.1.4).
- **Feedback.** A cursor change, a pill on the map ("Selecting stars") and the
  hint line. On touch, a tap selects and shows the info panel.
- **Reach.** 10 px for a mouse, 22 px radius for touch (a 44 px target; WCAG 2.5.8
  needs 24 px). When two or more candidates fall inside the reach, show a short
  disambiguation list (the stage view already handles "choices") instead of
  guessing the nearest.
- **Boundary (MAP.141).** Drop candidates outside the galaxy outline in the same
  filter that applies `enabled()` to picker layers (`createPicker` in
  `mappick.js`).
- **Layers.** Star mode adds a `points` layer at higher priority; Galaxy mode
  keeps today's order, in which stars deliberately cannot take the click
  (MAP.101). Star picks offer Bookmark and Waypoint (Boss); waypoints are the
  ones [course-routing.md](course-routing.md) section 4a describes.

### 5.2 Picking 70,000 to 1,000,000 points

Node 22, one thread, Xeon 2.8 GHz; a browser is of the same order. The
benchmark varies up to about 30% between runs [C, `pickbench.mjs`].

| Method | 70,000 pts | 100,000 | 1,000,000 | Notes |
|---|---|---|---|---|
| Current style: `Vector3.project(camera)` per entry plus `Math.hypot` (`mapcore.nearestOnScreen`) | 18 ms | 25-33 ms | 262-318 ms | Allocation-free but object-heavy |
| Typed-array loop, pre-multiplied matrix, no objects | 5.4 ms | 5.8-8.3 ms | 107 ms | Easy win |
| `Raycaster` on `Points` | 5.1 ms | 7-13 ms | 56-91 ms | Threshold in world units, not pixels |
| `PointsBVH` (three-mesh-bvh 0.9.16) | 2.0 ms | 1.3-1.6 ms | 3.6-3.9 ms | Build 0.6 s at 70k, 0.5-0.7 s at 100k, 3.5-4.4 s at 1M, on the main thread |

- With the 70,000-star cap per view (MAP.109), 5 ms on a click is fine, and the
  pick runs on a press, not every move. For hover, throttle to one query per
  frame while the pointer is still, or show the tooltip for the nearest star
  within reach after the pointer rests for 80 ms.
- If hover must be smooth: when the camera changes, project all stars once into
  a `Float32Array` and bin into 32 px cells; a query then touches about 9 cells.
  The 5 ms projection is paid per camera change, not per query.
- Raycaster and BVH use a world-space threshold; stars here are drawn a fixed
  number of pixels wide (`sizeAttenuation` off), so the equivalent threshold is
  `reachPx * distance * 2 tan(fov/2) / canvasHeight` and changes with distance.
  A screen-space test is simpler and exact.
- **Drop `three-mesh-bvh`.** The build costs more than the pick saves, PointsBVH
  supports neither worker build nor serialisation [S, its README], tiles stream
  in and out (`setStars` rebuilds), and the stars are not pickable today. MAP.102
  (done, PR #574) shipped tiles and camera-relative rendering with no BVH import
  anywhere in `src/html/static`. Use a plain `Raycaster` for the few-thousand-
  triangle nebula meshes.
- **GPU colour-id picking** (id pass to a 1x1 target through
  `camera.setViewOffset`, read back with `readRenderTargetPixelsAsync`, as in
  the three.js `webgl_interactive_cubes_gpu.html` example [S]) costs one extra
  draw of the points and an asynchronous readback, with no CPU loop. It wins
  above about 300,000 points or when picking must respect occlusion or exact
  sprite shape. Not timed (no GPU); the id shader needs the same `logdepthbuf`
  chunks as the star shader. Not needed below the cap.

### 5.3 Double-click (MAP.140)

- Do not depend on `dblclick`; touch behaviour varies (iOS Safari double-tap
  zoom, `touch-action`) [R]. The canvas already uses `touch-action: none` and
  its own pointer events (`createPointerControl` in `mapcontrol.js`,
  `clickOn: "pointerup"`). Add a recogniser there: a second `pointerup` within
  350 ms and 24 px of the first, on the same object, fires
  `onDoubleClick(event, pointerType)`. 24 px is the WCAG minimum target size, so
  a double-tap on a small target tolerates the same slop.
- Conflict with MAP.113 (PR #584: Escape or a click on the selected object
  clears the selection). In a double-click the first click clears and the second
  re-selects, so "double-click on a selected object" never sees a selected
  object. Fix: when a click hits the selected object, wait 350 ms before
  clearing; a second click in the window cancels the clear and goes. Selecting
  an unselected object stays instant. Cost: 350 ms delay when deselecting by
  click. The no-delay alternative (clear at once, remember the object for 350
  ms, reselect on the second click) flickers the ring and panel.
- Offer the same action without a gesture: a "Go to" button in the info panel
  (MAP.139's Details link) and Enter on the focused selection. A double-tap is
  a single-pointer gesture and allowed by WCAG 2.5.1, but nothing should be
  reachable only by it. The single click keeps selecting; the double-click goes
  to the next stage down or the object's page.

### 5.4 Context rendering (MAP.121, MAP.141)

- **Faint neighbours.** Alpha by distance from the focus, `a = a0 * exp(-d /
  lambda)` with `lambda` about one region width, floor 0.08 so the ring stays
  visible, ceiling 0.3; additive for stars, normal alpha for fills. A sector's
  26 neighbours (3x3x3 minus itself) are 26 small draws. Comets, rogue planets
  and asteroid fields are left out (Boss).
- **Above and below.** The same rule on the z axis, in the same style as the
  sides, pickable with a cheaper hit test (a box test per neighbour region, not
  the star list).
- **Blockers** (Boss: "the more they obstruct the view the more invisible they
  become"). With `s` the segment from the camera to the focus centre, `p` a
  neighbour's centre, `t` its projection parameter on `s` and `dist` the distance
  from `p` to `s`: if `0 < t < 1` and `dist < r_p + r_focus`, `fade =
  smoothstep(r_p, r_p + r_focus, dist)`, and opacity multiplies by `fade`
  (`fade = 0` means the neighbour lies exactly on the line of sight). Pure CPU,
  once per frame, tens of regions. Use `depthWrite: false` and sort by distance
  from the camera; with so few translucent layers the ordering error is small.
- One shared function, `contextOpacity(region, focus, camera)`, in the common
  engine piece MAP.121 calls for, not three copies.

### 5.5 Mini map as a second camera (MAP.75)

- One renderer, two passes: `renderer.autoClear = false`; main pass over the full
  canvas; then `setScissorTest(true)`, `setViewport(x, y, w, h)`,
  `setScissor(x, y, w, h)` (CSS pixels, origin bottom left), `clearDepth()`, render
  with the mini camera, `setScissorTest(false)`. The three.js `webgl_camera.html`
  example does this [S], and `systemmap.js` already scissors sub-viewports.
- Limit the mini camera with `camera.layers` (galaxy outline, a coarse star layer,
  current block and slab outlines) and skip the pass when neither camera nor
  selection changed (dirty flag).
- Hit testing: a transparent `div` with `role="group"` over the region runs its own
  pick against slab geometry; do not route clicks through the main camera. The slab
  buttons stay as the keyboard path (MAP.54).
- A 2D SVG or canvas mini map is simpler to make accessible and costs no GL pass; it
  is the fallback if the scissor path proves fiddly (log-depth buffer, DPR,
  resize). On phones the mini map sits below the map, as a second small canvas.

### 5.6 Star point rendering

The pipeline is sound (one `THREE.Points`, custom `ShaderMaterial`, pixel
`gl_PointSize`, camera-relative positions). Per-star attributes are 13 floats =
52 bytes [C], 3.6 MB at 70,000 stars; quantising colour, scalars and flags gives
about 20 to 24 bytes (1.5 MB), an optional saving for phones and tile swaps. Above
about 1e6 points use level of detail by tile (MAP.102, MAP.116), not a faster
picker. Overdraw is the real cost: cap the glow quad for faint stars.

Filed as MAP.147 (Boss, 2026-10-09 22:41Z): it starts with a measurement of the current payload and a comparison of options by Research Lane 3.

### 5.7 Habitable-worlds colouring (MAP.143)

Use a single-hue sequential ramp that survives greyscale (cividis or viridis style
[R]) in 4 to 5 bins with a legend; a sector with zero habitable worlds has no fill,
not the darkest bin. Blocks are never fully opaque (Boss, MAP.132).

## 6. Bookmarks and charted regions (UX.45 to UX.48)

Behaviour that stays valid is in [galaxy-drilldown-navigation.md](galaxy-drilldown-navigation.md)
section 8.2 (MAP.23, NAV.40, MAP.81): the breadcrumb star and page buttons, the Bookmarks
`<details>` menu, in-place rename (Enter saves, Escape cancels), delete, keys 1 to 9 while the map
has focus (not from a text box or select), the NAV page's Bookmarks select, and bookmarks that
carry a course pick (`data-pick`, `data-keep-name`, `data-keep-value`, `data-nav-url`).
`static/bookmarks.js` keeps per-browser `localStorage` under `planetgen.bookmarks.<db name>`, up to
100 entries `{name, kind, value, created, url?, sectorId?}`, identified by kind and value, with
every storage access wrapped, and wires whatever a page has without inline script (CSP). The new
items build on it.

- **UX.46 stamp.** The TODO says "stamped with the galaxy seed". The seed alone is not enough: an
  admin may regenerate the same seed, and `kind: "system"` bookmarks store database ids
  (`system:12`) that restart after a wipe, so old ids would point at different objects. Stamp with
  the seed plus a per-wipe instance id (for example `seed_hex + ":" + created_at_epoch`, or a random
  64-bit id written at wipe or creation). `galaxy_shape` has `galaxy_seed BINARY(16)` but no
  `created_at` or instance column (checked in `models.py`), so one small migration is needed. A page
  load that finds a mismatch drops the list (Boss's rule); offering a "recovered" export first is
  optional.
- **UX.47 manager page.** Sort by name, kind, date added and distance from the current position;
  inline rename; groups; delete with a 5-second "Deleted. Undo" toast rather than a confirm; open
  on map; export and import.
  - Order: a plain list with move up and move down buttons meets WCAG 2.5.7 (TanStack Table is
    overkill for 100 rows); drag handles are an optional enhancement (SortableJS 1.15.7, MIT, last
    published 2026-02-11 [C], or Pointer Events). Announce moves in an `aria-live="polite"` region.
  - Groups: one level of named `<details>` sections, a bookmark in at most one, a "Move to group"
    menu per row.
  - Export and import: JSON `{"format": "planetgen-bookmarks", "version": 1, "galaxy": "<stamp>",
    "groups": [...], "bookmarks": [...]}`, a Blob download and an `<input type="file">`. Validate on
    import, because `urlOf` returns `entry.url` straight into an `href`: cap size (1 MB), count and
    string lengths, take `kind` from a whitelist, check `value` against `ENDPOINT_RE` or a
    stage-URL pattern, accept `url` only as a same-origin relative path (`/galaxy`, `/system`...).
    CSP blocks a `javascript:` URL, but an off-site link is still an open-redirect nuisance. A stamp
    mismatch on import asks the user.
  - Show the count against the 100 cap (localStorage allows about 5 MB; 1,000 entries are about
    150 KB).
- **UX.48 charted regions.** Derived from the database, so not in localStorage: serve them from
  `GET /api/galaxy/regions` (cached, invalidated when a generation run ends) and list them in a
  read-only "Charted regions" section, with a local rename alias only.
  - Contiguity: face adjacency (6-connected) on the cylindrical sector lattice (ring, layer, slot;
    slot wraps; neighbouring rings can have different slot counts), so adjacency comes from
    `planetgen.galaxy` helpers, not index arithmetic. Edge-only or corner-only contact is not one
    region (26-connectivity is the alternative).
  - Maintenance: incremental union-find as sectors are charted; regions only merge, keeping the
    older id and any rename. A wipe drops everything with the galaxy.
  - Outline: collect faces between a charted sector and a non-charted neighbour, emit each face's
    four edges, drop edges whose two boundary faces are coplanar (interior seams); the rest is the
    silhouette. Draw as `LineSegments` in `--sel`, 2 px, dashed on the far side.
  - Fit: bounding sphere `R`, camera distance `R / sin(fov/2) * 1.15` along the current view
    direction, clamped to zoom limits; reuse MAP.78's per-level fit. With many regions draw only the
    focused one and the one under the pointer.

## 7. Job logs, progress and warnings (ADM.22, ADM.40, UX.3)

Built and still correct: a running job's output streams over Server-Sent Events
into an Xterm.js terminal (`static/joblog.js`; scrollback 20,000 lines, `<pre>`
fallback), with progress in a native `<progress>` (`partials/job_status.html`).
The stream route (`generate_page.generate_job_stream`) sends `log` events whose id
is the byte offset reached, `state` events when the job view changes, and `done`.
It closes after `STREAM_SECONDS = 40` because Apache's `WSGIDaemonProcess
request-timeout` (60 s in `examples/apache/planetgen.conf.example`) restarts the
daemon for a longer request; `EventSource` reconnects with `Last-Event-ID`, so
nothing is lost. A first connection replays at most the last 256 KB; keep-alives go
out every 15 s. Since the browser reports each deliberate close as an error,
`connectionnote.js` (`quietReconnect`, ADM.40) shows "Connection lost;
reconnecting..." only after a grace period. "Copy log" fetches the full output from
the download route (ADM.25).

UX.3 (the banner for every visitor, ETA rounded up to the hour) is not built. Its
2026-10-07 plan reads the RQ job's published progress (ADM.22), not
`progress.json`. Carry from the job pages: single-line messages, a failed job's log
stays open until the user continues, tracebacks in full with a Copy button, and
ETA wording that follows ADM.39's time-left estimate so the banner and the queue
page agree. The TODO's open points (every page or only map and list pages;
command-line runs; wording with no ETA yet) are unchanged by this research.

## 8. Admin-editable expected density (MAP.119)

`galaxy_shape.expected_system_count_at_density_1` (E, 6.31 systems for a 4 pc
sector at relative density 1) calibrates the whole galaxy: a sector qualifies when
`E * relative_density(center) >= 1` ([galaxy-disk-density.md](galaxy-disk-density.md)),
and a qualifying sector's own `relative_density` is its `--density` multiplier at
fill. [reproducible-galaxies.md](reproducible-galaxies.md) requires that seed plus
version key plus recorded admin changes rebuild the same galaxy, with changes stored
as a net difference by stable address, not database id. So:

1. Never edit E or the skeleton. The admin edits a per-region multiplier `m`
   (default 1.0) applied to the `--density` used when sectors there are filled.
2. Overrides are append-only data: `density_overrides(id, scope_kind, scope_address,
   multiplier, created_at, admin_id, reason, supersedes_id)`. A change inserts a
   row, undo inserts the reverse, nothing is updated in place. The scope is an
   address (block or slab by ring or layer range), never a database id.
3. Every filled sector records `density_multiplier_used` and `override_id`;
   reproduction replays the override rows in order, then regenerates.
4. Only unfilled sectors are affected ("later fills use" it); the UI says how many
   sectors in view are filled and unaffected.
5. The settings file carries the rows (reproducibility doc, section 7) as part of
   the net difference: a new admin change kind for GEN.59 and GEN.61.
6. Range 0.1 to 10 on a log scale, beyond that a typed override with a warning;
   store a decimal string or rational.
7. Write an entry to the admin activity log (`planetgen.admin.activity_log`).

Presentation: guests and users see a static readout ("Expected density: 6.3 systems
per 4 pc sector (x1.0)"). Admins also get a number field, a native log-scale
`<input type="range">` (`sl-range` is not vendored; the number field is the WCAG
2.5.7 non-drag alternative), "Reset to 1.0" and quick steps (x0.5, x2). The map
re-colours with the draft multiplier at once, since the client computes density
analytically for the block colours (`galaxymap3d.js`), with a preview line ("Expected
systems in view: 1.9 million to 3.8 million (+100%); 0 filled sectors affected, 412
unfilled"). Applying is a deliberate `sl-dialog` confirm listing scope, old and new
value and counts; undo by a session toast plus a history list with a Revert button.
An existing override shows as a chip ("x2.0 by Boss, 2026-10-09"); overlapping scopes
resolve narrowest-wins and the dialog says so.

## 9. Boss's uploaded UX documents: what still holds

"Web UX Development Notes.md" (architecture proposal) and "Web UX and Job Management
Guide.md" (migration guide) are generic. Neither mentions Select mode, bookmarks or
job warnings, so the procedures below are those that bear on libraries, logs,
progress and the 3D view. The source documents are not edited.

Still valid and in use:

- Xterm.js fed by Server-Sent Events instead of polling, with a buffer limit and a
  plain fallback (Notes section 3, Guide Phase 2; section 7 above).
- A native `<progress>` with ARIA instead of hand-built DOM updates (Notes section 4).
- Jinja renders a structural shell with `data-*` attributes and static scripts hydrate
  it (Guide "Integration"): `data-bookmark-*`, `data-job-*`, `data-bookmarks-menu`,
  the data tables. Keep the pattern for new components.
- TanStack Table with Virtual (UX.41, PR #606) and Shoelace (UX.40) as proposed.
- Camera-relative rendering against float32 jitter (Notes section 2, Guide "Coordinate
  frame shift"): MAP.102, `starFrame`.
- Standalone ES module builds rather than a runtime bundler (Guide "Asset bundling"):
  esbuild is used only offline to bundle three and TanStack into single files.

Superseded or done differently:

- Guide Phase 1 and its SQLite section (`ProcessPoolExecutor`, WAL, a status queue):
  Boss chose RQ on Redis ([library-migration.md](library-migration.md) section 3), and
  the database is MySQL or MariaDB, so SQLite locking and "SQLite query overhead" do
  not apply.
- Guide Phase 2 step 1 (a `queue.Queue` per job id in `log.py`): an in-process queue
  does not cross the processes Apache and the RQ workers run in. The log is a file the
  stream route reads by byte offset.
- Guide "WSGI starvation" (eventlet or gevent workers): deployment is Apache with
  threaded mod_wsgi daemons plus the 40 s stream cycle.
- OGC 3D Tiles, `3d-tiles-renderer`, `three-mesh-bvh` and a BVH export from
  `galaxyGeometry.py` (Notes section 2, Guide Phase 4): tiles with level of detail
  and camera-relative rendering shipped as MAP.102 without either library; a BVH was
  measured and rejected for picking (section 5.2).
- Notes claim that Xterm.js gives "WebGL/Canvas virtual scrolling ... steady 60 FPS
  regardless of throughput": only `xterm.mjs` and the fit addon are vendored, with no
  WebGL renderer addon, and the 20,000-line scrollback is the real limit. Unverified.
- File names (`sectormap.js`, `workQueue.py`, `galaxyGeometry.py`) predate the
  reorganisation; use the `planetgen.*` paths in
  [library-migration.md](library-migration.md) section 6.

## Evidence notes

[C] results rest on scripts and registry reads dated 2026-10-09. The Node
benchmark is single-thread V8, not a browser, and varies by up to about 30%
between runs. Contrast numbers use the WCAG relative-luminance formula and sRGB
`color-mix` emulation for Shoelace's derived tokens; they are not screenshots.
The environment could read search-result text and raw files only, not papers or
the Web Awesome migration guide.

[R] items to verify when paper access is allowed:

- iOS Safari and Android Chrome `dblclick` behaviour on touch with
  `touch-action: none` (the design avoids relying on it).
- Aviation meanings of magenta, cyan, green and white (AC 25-11 / 25.1322 text
  not retrieved; red, amber and advisory rules were seen).
- B612's Airbus and ENAC origin (licence and file sizes are [C]).
- Design traits of Star Citizen, Starfield, Mass Effect, Stellaris, Celestia,
  SpaceEngine, Gaia Sky, Stellarium Web and Universe Sandbox.
- GPU timings: ray-march and id-pick cost, mobile texture rates, the 64 px
  point-size limit.
- Browser baselines for CSS relative colour syntax (Web Awesome themes) and
  import-map integrity.
- Web Awesome Pro pricing and the exact `sl-` to `wa-` renames.
- Licences of the packages inside the Shoelace bundle (Lit, floating-ui,
  tinycolor).
- Okabe and Ito (2008) as the palette origin, and "below delta E 10 is hard to
  tell apart".
- Autofill behaviour inside shadow roots; Playwright piercing of open shadow
  roots.
- The cividis or viridis ramp as colour-vision safe, for MAP.143.

## Sources

- npm registry JSON read 2026-10-09 (`https://registry.npmjs.org/<name>`) for
  Shoelace, Web Awesome, three, @xterm/xterm, @xterm/addon-fit, the old xterm
  packages, @tanstack/table-core and virtual-core, three-mesh-bvh, @material/web,
  @lion/ui, Spectrum, Carbon, Pico, Open Props, sortablejs and
  @fontsource/b612-mono; tarballs of three 0.186.0 and 0.186.1, Web Awesome
  3.14.0, Shoelace 2.20.1 and b612-mono 5.3.0 unpacked and read.
- Raw GitHub: `shoelace-style/shoelace` README (sunset notice),
  `shoelace-style/webawesome` README and LICENSE.md, `material-components/material-web`
  README, `gkjohnson/three-mesh-bvh` README, `mrdoob/three.js` examples
  `webgl_interactive_cubes_gpu.html` and `webgl_camera.html`, `nasa/openmct`
  `_constants-espresso.scss` and `_constants-snow.scss`.
- Search results: https://blog.fontawesome.com/introducing-web-awesome/ ;
  https://blog.fontawesome.com/how-does-web-awesome-stack-up-against-shoelace/ ;
  https://threejs.org/manual/en/webgpurenderer ;
  https://discourse.threejs.org/t/three-js-now-imports-three-core-js/75187 ;
  https://sbcode.net/threejs/importmap/ ;
  https://celestiaproject.space/forum/viewtopic.php?p=86447 ;
  https://tonisagrista.com/blog/2024/rendering-aurorae-nebulae/ ;
  https://arxiv.org/pdf/1609.05344 ; Frontier forum HUD colour threads ;
  https://www.faa.gov/sites/faa.gov/files/AC25-11.pdf and EASA alerting documents ;
  https://dequeuniversity.com/resources/wcag-2.2 ;
  https://www.boia.org/blog/whats-new-in-wcag-2.2-new-success-criteria-for-digital-accessibility ;
  https://ntrs.nasa.gov/api/citations/20160006385/downloads/20160006385.pdf .
- Repo files read: the map, bookmark, job-log and theme files under
  `src/html/static/`, `vendor/THIRD_PARTY_NOTICES.txt`, `scripts/vendor_shoelace.py`,
  `planetgen/galaxy/nebula_shape.py`, `planetgen/web/__init__.py` (CSP),
  `planetgen/web/generate_page.py`, `planetgen/db/models.py`,
  `src/tests/test_web_browser_colors.py`. The research scripts (contrast, colour
  blindness, nebula profile, Shoelace size, picking benchmark, tree-shake test)
  lived in the research session's scratchpad, not in the repo.
