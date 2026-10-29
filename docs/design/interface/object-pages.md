# Object pages (UX.168 to UX.178)

Every star, planet, moon, asteroid belt and comet has a page of its own at `/object/<kind>:<Object ID>` (`web/item_pages.py`).
The system page and the sector page list their contents as summary trees whose rows link to those pages. This note says how
the parts fit and the rules a change has to keep.

## The page

In order: the Location block (galaxy, sector, system, the star the body orbits and, for a moon, its planet), links to the
previous and next object at the same level (the planets that share a star in their order out, a planet's moons, the system's
stars, the belts or comets of a star), the rendered example (star, planet and moon), the property sheet, an "About" panel
(the description's paragraph, a life world's dominant species, the facilities the body hosts) and "Belonging to it", the
summary tree of the body's children. A planet or moon also shows its PHI-6 symbols and habitability chips under the title,
and an admin has that body's Admin menu in the action bar (Regenerate, Change class, Delete; the page takes the POST).

## The region (not built, UX.231)

The at-a-glance block will say "in Vang Zarh Region" for a sector or object whose position lies in a charted region (a 9-block with a `region_stats` row, `../galaxy/named-regions.md`), and nothing when it lies in none. Until UX.231 the pages say nothing of regions.

## The property sheet

`describe/sheet.py` is the one registry of what a page shows about a body: per kind (`star`, `planet`, `moon`, `belt`,
`comet` and the eight standalone phenomenon kinds, UX.175) an ordered list of `PropertySpec(key, label, group, quantity, column, solar)`. The page (`item.html`) and
`GET /api/objects/<ref>/sheet` both read it, so the website, the API and the tests agree. A value the row does not hold
is returned as "not stored", never left out, so an empty line means no data and not a bug.

**A new column needs a sheet row.** `tests/test_property_sheet.py` walks every column of `stars`, `planets`, `moons`,
`asteroid_belts`, `comets` and the eight phenomenon tables (from `db.models.metadata`, with the unit from `db/column_units.py`) and fails for a column
that is neither a `PropertySpec` of its kind nor on that kind's `NOT_SHOWN` list with a reason (the ids, links and
bookkeeping a reader never needs). A lane that adds a column to one of those tables adds its sheet row, or its `NOT_SHOWN`
reason, in the same pull request.

## Phenomena

A nebula, asteroid field, black hole, neutron star, supernova remnant, rogue planet, interstellar comet or quasar keeps its
own page, `/phenomenon/<type>/<id>` (`web/system_pages.py`, `phenomenon.html`): its view, a curated data table with the
derived values (Einstein radius, galactic speed) and class links, the Location block, what a nebula or remnant holds, and
since UX.175 an "All stored properties" panel that is the same property sheet as a body's, served by the same API. Not yet
on the phenomenon pages: the paragraph text (the phenomenon description is built from a freshly generated object and
cannot yet be rebuilt from a stored row), neighbour links and the tree. Facilities, star clusters and regions are not
object kinds with an Object ID yet, so they have no page of their own.

## The summary tree

`web/lib/objecttree.py` draws the rows and `static/objecttree.js` makes them live. A row is a toggle button, the kind's
mark, the name (a link to the object's page), a summary line (star: spectral type, role, mass, luminosity; planet and
moon: class, distance, period and the PHI-6 symbols; belt: density and range; comet: kind, activity and period; system:
single or binary and the star types) and, for a row with children, the counts below it. Rules:

- No link sits inside a button or a `<summary>` (assistive technology reads that as a nested control). The toggle is its
  own button, hidden until the script runs; without script every level is simply listed and a row whose children are not
  on the page holds a link to its own page.
- A level loads when its row opens, from `GET /object/<ref>/tree` (the children API, `db/children.py`, keyset paging, 50
  rows at a time) and a longer level ends in a "Show more" button. A system's stars are drawn open onto their own bodies.
- The keys are the WAI-ARIA tree pattern: Up and Down, Right (open, then first child), Left (close, then the parent),
  Home and End; one name is in the Tab order at a time.
- The sector page's Contents table puts the same tree under a system's name (`static/datatable.js` calls
  `objecttree.js`'s `embedTree`); an opened row is kept across the table's redraws.

## The rendered example

`web/lib/bodyview.py` builds the look record of the body (`web/maps/bodylook.py`, the record the System Map paints from) and
`static/bodyview.js` paints it with `solidpaint.js` or `starpaint.js` and shades it as a lit sphere on a 2D canvas, with
the atmosphere's limb glow or a star's corona. A giant, which the painters do not cover yet, is its tinted disc; a belt
and a comet have no example yet. The caption is `look_note`, a stand-in until the colour account of UX.158 replaces it. See
[body-rendering.md](body-rendering.md).

## Links

Every link to a body is `/object/<ref>`: the Document references (`item_pages._ref_link`), the search results
(`system_pages.object_url`), the map popups' View button (`systemnav.js`), and the tree rows. The old `/system/<id>#planet-12`
anchors, `body_anchor`, `static/bodyanchor.js` and the row-highlight styles are gone; `tests/test_item_pages.py` fails for
any source that links a body to a system anchor.

## Tests

`test_property_sheet.py` (the sheet), `test_object_children.py` (the children API), `test_object_tree.py` and
`test_item_pages.py` (markup and pages), `test_object_page_walk.py` (every kind on a generated database),
`test_web_browser_object_tree.py` (the tree in Chromium, with axe-core), `test_web_browser_bodyview.py` (the picture) and
`test_web_a11y.py` (every page at phone and desktop width).
