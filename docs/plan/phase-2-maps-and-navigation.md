# Phase 2: Maps and navigation

This is one phase of the plan drawn up on 2026-10-01 from Boss's list
of that evening and his research notes. [docs/TODO.md](../TODO.md) is
the master file: it holds each item's full text (what's wrong, where to
look, what "done" means), and its "Plan: phases" section indexes every
phase. This file holds what the phase needs beyond that: its goal, the
order and dependencies, how to split it into build threads, and the
research notes that apply. Where the two disagree, TODO.md wins; when an
item ships, it leaves TODO.md and its row here is deleted in the same PR.

## Goal

Rebuild the Galaxy Map's selection around the arc pick (a 3D galaxy with
no sector lines, an arc, then a slab, then segments down to a sector),
color sectors and blocks by what is in them, join the Galaxy, Sector and
System maps on one engine, and build the shared picker and courses on
top. Starts once phase 1's map groundwork (MAP.63, MAP.64, TEST.70) and
object references (NAV.7) are in.

## Items

### Galaxy Map selection: the arc pick

| ID | Item | Parent |
|---|---|---|
| MAP.60 | Galaxy Map scale readout: one scale line |  |
| MAP.55 | Galaxy Map buttons: a menu, with only back, forward, up, reset and bookmark showing |  |
| MAP.85 (new) | The galaxy pick is an arc, on a 3D galaxy with no sector lines |  |
| MAP.52 | Galaxy Map highlights the wrong area; pick a 40-degree wedge around the cursor (bug) |  |
| MAP.56 | Drop the 3x3 block pick: select a slab, zoom in, select a segment (bug) |  |
| MAP.53 | Rotate a zoomed-in wedge, and zoom it to fit the window (bug) |  |
| MAP.58 | Galaxy Map zoom limits: a short manual range on the galaxy wedge, locked below it |  |
| MAP.78 | Zooming into a wedge must show the whole wedge at every drill-down level (bug) |  |
| MAP.54 | Slab leader lines instead of the slab slider (bug) |  |
| MAP.76 | Leader-line layout | MAP.54 |
| MAP.59 | Make it plain that a zoomed-in slab is a slab, not a wedge |  |
| MAP.75 | The mini map as a second engine view | MAP.59 |
| MAP.77 | Galaxy Map draws block divisions inside a picked slab before zooming to it (bug) |  |
| MAP.80 | Sector-level zoom on the Galaxy Map should show almost every star in the sector (bug) |  |

One build thread, in this order, because every item rewrites
`galaxystageview.js` and `galaxystages.js`: MAP.60 and MAP.55 first
(small, they free space and remove the Wedges button and lines); then
MAP.85 with MAP.52 (the arc pick replaces the wedge pick; MAP.52's
bearing width and snapping carry over); MAP.56 (arc, slab, segment
ladder); MAP.53 and MAP.58 with MAP.78 (rotate, fit, zoom limits, all on
MAP.64's controller); MAP.54 with MAP.76 (slab buttons and leader
lines); MAP.59 with MAP.75 (ghost, mini map, header); MAP.77 (which
lines show at each level); MAP.80 (sector-level detail). Each of these
items carries an "Arc pick (MAP.85)" note in TODO.md saying how the arc
changes it.

### Sector and block colors

| ID | Item | Parent |
|---|---|---|
| MAP.86 (new) | Sector and block colors from what is in them: filled sectors translucent (bug) |  |

After MAP.85 removes the lines, color carries the structure. Needs a
per-sector color, saturation and lightness computed when a sector is
saved and served with the tiles (an API field, and maybe a column), then
the averaging in `galaxyblocks.js`.

### One map engine

| ID | Item | Parent |
|---|---|---|
| MAP.61 | One map engine and control set for the Galaxy Map and the Sector Map |  |
| MAP.65 | One picking, hover and info-panel layer | MAP.61 |
| MAP.66 | The sector as the drill-down's last stage, on the same page | MAP.61 |
| MAP.67 | One URL and history scheme for every level | MAP.61 |
| MAP.68 | Remove the old Sector Map code | MAP.61 |

After the selection work: MAP.65 (picking, hover, info panel), MAP.66
(the sector as the drill-down's last stage, needs TEST.70), MAP.67 (one
URL scheme), then MAP.68 removes the old Sector Map code.

### Sector Map objects

| ID | Item | Parent |
|---|---|---|
| MAP.79 | Rogue planets clog the Sector Map: dim them, and a show/hide button per kind of object (bug) |  |
| MAP.82 (new) | Unmarked rogue planets barely visible (bug) | MAP.79 |
| MAP.83 (new) | The "Mark rogue planets" button shows when it is on (bug) | MAP.79 |
| MAP.84 (new) | Marked rogue planets grow and become clickable; unmarked ones stay small (bug) | MAP.79 |
| MAP.81 | Ctrl+1 to Ctrl+9 bookmark keys clash with the browser's tab switching (bug) |  |
| MAP.87 (new) | Stars on the Sector Map need to be brighter, most of all the dim ones (bug) |  |

MAP.79's per-kind toggles and dimming, with its three rogue planet
subitems (dim when unmarked, the Mark button's highlight, size and
clickability when marked). They touch `sectormap.js` and
`lib/starmap.py`; if MAP.61 is under way, build them on the shared
engine. MAP.81 (bookmark keys) needs Boss's pick of keys. MAP.87 (brighter dim
stars) changes only `_star_light` in `lib/starmap.py` and can go any
time.

### 3D system view

| ID | Item | Parent |
|---|---|---|
| MAP.62 | A full 3D star system view with a free camera |  |
| MAP.69 | A system scene endpoint with 3D orbits | MAP.62 |
| MAP.70 | Positions at any time | MAP.62 |
| MAP.71 | Scale modes that keep everything visible | MAP.62 |
| MAP.72 | Rendering at system scale | MAP.62 |
| MAP.73 | Free camera on the shared engine | MAP.62 |
| MAP.74 | The 3D view on the system page, the flat diagram kept | MAP.62 |

Its own thread once MAP.61's engine exists: MAP.69 (scene endpoint) and
MAP.70 (positions at any time) first; NAV.27 needs MAP.70.

### Picker and courses

| ID | Item | Parent |
|---|---|---|
| NAV.3 | One shared picker for the Galaxy, Sector and System displays |  |
| NAV.13 | A picker module: select, step out, step in, step sideways | NAV.3 |
| NAV.14 | One breadcrumb for every level | NAV.3 |
| NAV.15 | Pick mode everywhere | NAV.3 |
| NAV.16 | NAV endpoints can be any object | NAV.3 |
| NAV.29 | Replace "Nav from here" and "Nav to here" with "Start Here" and "End Here" while picking (bug) |  |
| NAV.30 | Hide "View phenomenon" and "View system" links while picking a course (bug) |  |
| NAV.31 | Galaxy wedges don't highlight on the navigation screens (bug) |  |
| NAV.32 | Every Galaxy and Sector Map control works on the navigation screens (bug) |  |
| NAV.33 | After picking one end of a course, stay at that zoom level (bug) |  |
| NAV.5 | Show a course on the Galaxy Map |  |
| NAV.20 | Draw the direct line and the route apart | NAV.5 |
| NAV.21 | Fit the view to the whole course | NAV.5 |
| NAV.22 | Courses inside a sector and a system | NAV.5 |
| NAV.23 | Open a saved course on the map | NAV.5 |
| NAV.4 | Save a course |  |
| NAV.17 | A saved course record with both forms | NAV.4 |
| NAV.18 | Save, list, open, rename and delete, per browser | NAV.4 |
| NAV.10 | Routing that scales past a few thousand systems |  |
| NAV.11 | Travel times for the system-to-system route too | NAV.10 |
| NAV.12 | A maximum hop length (open question) | NAV.10 |
| NAV.6 | Courses that steer clear of gravity wells |  |
| NAV.24 | A keep-out radius for every kind of object | NAV.6 |
| NAV.25 | Find the obstacles along a path | NAV.6 |
| NAV.26 | Bend the path around keep-out spheres | NAV.6 |
| NAV.27 | Moving bodies inside a system | NAV.6 |
| NAV.28 | Show and save the adjusted course | NAV.6 |

NAV.3 (the shared picker, `picker.js`) with NAV.13 to NAV.16 first, on
NAV.7's references; the pick-mode bugs NAV.29 to NAV.33 come with
NAV.15. Then drawing courses (NAV.5, NAV.20 to NAV.23), saving them per
browser (NAV.4, NAV.17, NAV.18; NAV.19 waits for accounts in phase 4),
routing at scale (NAV.10 to NAV.12, NAV.12 needs Boss's answer), and
courses that bend around gravity wells (NAV.6, NAV.24 to NAV.28).

### Interface cleanup

| ID | Item | Parent |
|---|---|---|
| UX.21 | Clean up the web interface: overlapping buttons and dead controls (bug) |  |

After MAP.55 and MAP.60, which already remove some dead controls.

## Research notes

Boss's research notes (kept in the project's shared files under `todo-tasks/research/`) proposed fixes and numbers. They were checked against the code on 2026-10-01; where they were wrong about the code, the correction is given. Their numbers are starting points to tune, not requirements.

- **The arc pick (MAP.85).** The notes lay it out in four stages: (1)
  the whole galaxy as a particle field of its stars and spiral arms,
  with no grid, block or sector lines; (2) hovering projects an arc
  through the full height of the disk, highlighted, with faint
  boundaries of the neighboring arcs; (3) clicking flies the camera to
  the arc, where the user picks a height band (a slab); (4) the slab's
  sectors appear for the next pick. Note the word: "arc" was MAP.19's
  name for one cell of the old 3x3 region pick
  ([galaxy-drilldown-navigation.md](../design/galaxy-drilldown-navigation.md));
  it now means this first pick, and the design doc needs the same
  change when MAP.85 is built.
- **Sector colors (MAP.86).** The notes' starting values:

  | State | Opacity | Saturation | Hue from |
  |---|---|---|---|
  | Unfilled | 0.03 | 0 | neutral dark grey |
  | Filled, empty | 0.15 | 0.10 | the site's accent |
  | Sparse (1 to 10 stars) | 0.20 | 0.35 | mostly M dwarfs, red-amber |
  | Dense (50+ stars) | 0.45 | 0.85 | G and A stars, yellow-white |
  | Densest | 0.70 | 1.00 | O and B stars, blue |

  Saturation scales with star density, hue is the luminosity-weighted
  average of the stars' colors by temperature, lightness scales with the
  log of total luminosity, and a block's color, saturation and opacity
  are the plain average of its sectors'. Boss asked for filled sectors
  to stay translucent, "just a hair more solid" than unfilled, so the
  top of that opacity range is a ceiling to try, not a target.
- **Rogue planets on the Sector Map (MAP.82 to MAP.84).** The notes'
  values:

  | Property | Unmarked (default) | Marked |
  |---|---|---|
  | Radius | 1.5 px | 5 px |
  | Opacity | 0.20 | 1.0 |
  | Color | muted grey (#4A5568) | bright with a glow ring (#00E5FF) |
  | Hit radius | 3 px | 12 px |
  | Button | plain | highlighted, `aria-pressed="true"` |

  The button exists today (`lib/starmap.py`, `sectormap.js`) but starts
  pressed and has no pressed style.

## Build threads

Each thread is briefed with its exact item IDs and takes no others.

1. Galaxy Map selection (the one ordered thread above), then MAP.86.
2. Sector Map objects: MAP.79, MAP.82 to MAP.84, MAP.87.
3. One map engine: MAP.61, MAP.65 to MAP.68 (after thread 1).
4. 3D system view: MAP.62, MAP.69 to MAP.74 (after thread 3 starts).
5. Picker and courses: NAV items in the order above.
6. UX.21 at the end.

## Open questions for Boss

- MAP.85: is an arc a third of the radius by about 40 degrees (the
  default), the whole center-to-edge wedge, or smaller?
- MAP.52, MAP.53, MAP.54, MAP.55, MAP.56, MAP.58: the open questions in
  their TODO.md entries.
- MAP.81: which bookmark keys?
- NAV.12: the maximum hop length.
