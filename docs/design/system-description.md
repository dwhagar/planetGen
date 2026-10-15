# System description: one model, many renderers

How planetGen turns a star system, a planet, a moon, a belt or a comet into the text and tables people read, so that
the web page, the API, the command line and the wiki exports all say the same thing from one place.

Status: design, 2026-10-15. Decisions marked "Boss" are his; everything else is a recommendation. Written by the
Foundations Wording lane. Informs: UX.30, the work list in `/mnt/project-files/wording/text-output-audit-2026-10-15.md`
(W2 to W12).

Boss (2026-10-15 04:21Z): "sweep all planet and system text output ... Make sure wording matches the habitability index
and vice versa ... we're stripping out the markdown from the web site display in favor of more advanced tables like we're
using, and generating markdown on the fly as needed; we can output markdown / wikitext for the user on request but it is
no longer the default." Earlier (UX.30, 2026-10-01): "Rework planet information displays to pull away from the markdown
render and instead fits with the modern web UI." The system overview stays prose (Boss, 2026-10-10 19:12Z).

## 1. What is wrong today

- Every model class writes its own prose and tables: `Star`, `BinaryStarProxy`, `WideBinary`, `Planet` (planets and
  moons), `AsteroidBelt`, `Comet`, `StarSystem.summary_paragraph`, and six phenomenon classes, each in
  `to_paragraph_list()`. The output format is a mutable flag on the shared config (`SystemConfig.MARKDOWN`), flipped
  around each render by `db/render.py`, and saved in the run-config row (`markdown` column).
- The web's expanded body rows are that Markdown, split per body by `/api/systems/<id>/sections`, then converted to HTML
  by `web/lib/mdrender.py` (the Python `markdown` package cut down to headers, pipe tables and `<sup>`).
- Units, header levels, wikitext template names and sentences live in the model classes. Quantities in sentences are
  formatted by four different helpers and bypass the unit ladders (`util/ladders.py`).
- Life timelines (`*_evolutionary_paragraphs`), `flavor_text` and `hab_note` are stored as finished prose, and the code
  reads the stage back out of the text.

The audit file lists every surface and every finding with `file:line`.

## 2. The shape

```
model object graph / stored rows
        |  describe()          (pure, no format, no HTML, no Markdown)
        v
   Document  (typed, JSON-serialisable)
        |
        +--> HTML renderer      (web: tables, chips, .qty spans)
        +--> Markdown renderer  (on request)
        +--> wikitext renderer  (on request, wiki upload)
        +--> text renderer      (command line)
        +--> JSON               (API, `GET /api/systems/<id>/description`)
```

`describe()` lives in a new package `planetgen/describe/` (generation side: it may import the model and `physics`). It is
a pure function of the system. The web app does not import it (the web app never imports `numba`, and the generation
package pulls it in); the web gets the Document as JSON from the API and renders HTML from that, the way it already
gets rows. The one-off system page, which today runs the CLI in a subprocess, asks the CLI for the JSON form and
renders it the same way. The Markdown, wikitext and text renderers are plain Python in `planetgen/describe/` because the
API (`/text`), the CLI and the wiki upload all run on the generation side.

The old per-class `to_paragraph_list()`, `__str__` renderers, `properties_to_string`, `to_scientific_notation`'s
format switch, `db/render.py` and `SystemConfig.MARKDOWN` are deleted when the new path lands (no shims; Boss, 2026-10-07).

## 3. The Document

A Document is a tree of sections. Every node is a plain dict or list in JSON; in Python they are small frozen
dataclasses in `describe/model.py`.

```
Document   {"kind": "system", "ref": "system:12", "title": "40 Eridani", "sections": [Section, ...]}
Section    {"kind": "system|star|pair|planet|moon|belt|comet",
            "ref": "planet:31",            # the object reference the site already uses (anchors, links)
            "title": "40 Eridani II",
            "blocks": [Block, ...],
            "children": [Section, ...]}    # a planet's moons; a wide pair's stars hold their own bodies
Block      one of
  Facts     {"type": "facts", "caption": "Orbit", "rows": [{"key": "distance", "label": "Distance", "value": Inline}]}
  Ratings   {"type": "ratings", "overall": {...}, "rows": [{"domain": "pressure", "tier": "Red", "score": 0.0,
                                                            "text": [Inline, ...]}], "gear": [Inline, ...]}
  Composition {"type": "composition", "label": "Atmosphere", "parts": [{"name": "nitrogen", "share": 0.78}, ...]}
  Prose     {"type": "prose", "role": "overview|summary|life|timeline|flavour", "text": [Inline, ...]}
Inline     a string, or one of
  {"qty": {"kind": "pressure", "si": 0.0}}          a quantity: kind (the ladder family), SI value
  {"ref": "planet:31", "name": "40 Eridani II"}     a link to another object
  {"strong": "Pressure: Red"}                       emphasis, rendered per format
```

Rules:

1. **No format in the model.** No `#`, `**`, `<sup>`, `{{Exp|..}}`, `'''` and no Unicode superscripts. Emphasis and
   links are Inline nodes the renderers turn into whatever the format has.
2. **Quantities are values, never strings.** A quantity names its family (`distance`, `speed`, `period`, `mass`,
   `temperature`, `pressure`, `gravity`, `dose`, `luminosity`, `density`...) and its SI value. Each renderer calls the
   one ladder for that family (`util/ladders.py` for the text forms; `unitladder.js` for the browser, which is how the
   unit preference UX.78 can later re-render every `.qty` span from `data-si`). Nothing in `describe/` formats a
   number itself, and a test fails on a raw `:.2f` or `f"{x} km"` in that package.
3. **A fact row has a stable `key`.** The key is the old wikitext parameter name where there was one (`type`, `radius`,
   `mass`, `temp`, `lum`, `hab`, `orbit`, `loc`, `wobble`, `class`, `distance`, `period`, `speed`, `gravity`,
   `moon_wobble`, `separation`, `mutual_orbit`, `eccentricity`, `primary_limit`, `secondary_limit`), so the wikitext
   renderer writes the same templates and a wiki that already holds those pages keeps working. The label is the display
   name (the Markdown key map, `Property | Value`).
4. **Text is short and complete.** A Prose or Ratings sentence is a finished sentence in the site's voice (section 6).
   Wording is chosen in `describe/`, not in the model classes and not in the renderers.
5. **Everything is derivable.** A Document is a function of the stored rows (or a fresh system). Nothing is cached in
   the database, so rewording changes every page with no regeneration.

### 3.1 What goes in each section

| Section | Facts | Other blocks |
|---|---|---|
| system | the overview prose only (stays prose: counts, habitable worlds, heliosphere, edge of the system; the cloud sentence) | `Prose(role=overview)`, `Prose(role=flavour)` when present |
| star | type, radius, mass, temperature, luminosity, habitable zone, galactic orbit, planetary wobble, location | `Prose(role=summary)`: age, lifespan, the one sentence about what that means for life |
| pair (close binary) | type, mass, luminosity, habitable zone, separation, mutual orbit, barycenter offset, galactic orbit, location | `Prose`: the pair's age |
| pair (wide binary) | separation, eccentricity, mutual orbit, barycenter offset, each star's planetary limit | none |
| planet, moon | class, distance, period, speed, radius, gravity, moon wobble; surface temperature, pressure (so the 0 kPa of an airless world is a fact, not a gap) | `Composition` (atmosphere, surface), `Ratings` (rocky bodies), `Prose(life)`, `Prose(timeline)`, `Prose(flavour)` |
| gas giant | the same facts; internal pressure and temperature | `Composition`, `Prose` for the depth sentences |
| belt | inner and outer edge, density | `Composition` |
| comet | type, perihelion, eccentricity, period, current distance and speed, nucleus size, activity | `Composition` |

The table rows and prose are the ones the old text carried; W2 adds a coverage test that every old fact key and every
old sentence topic has a home, so nothing is lost silently. Phenomena (nebula, quasar, rogue planet, asteroid field,
compact remnant, supernova remnant) are described by the same model from W10; the web's phenomenon page already shows
their facts as tables (`FIELD_SPECS`) and only the CLI shows their prose.

### 3.2 The rating table

The six PHI-6 domains are a table, not a paragraph:

| Domain | Rating | Why |
|---|---|---|
| Pressure | Red (0.00) | The air pressure is 0 kPa, so a pressure suit is needed. |
| Temperature | Red (0.00) | The air temperature is -227 C, so an environment suit is needed. |
| Air | Red (0.0) | There is no atmosphere (0.0 kPa). |
| Water | Red (0.00) | There is no liquid water to support life. |
| Radiation | Yellow (0.45) | The radiation dose is a mean of 742 mSv a year, inside the Yellow band (up to 1 Sv a year). |
| Suitability | Red (0.20) | There is no atmosphere (0.0 kPa), no food or water. Shelter is potentially available. |

with the overall line above it ("PHI-6 habitability is 0.00, and a visitor needs a sealed suit.") and the equipment label.
These rows are Boss's text (GEN.216, 2026-10-15 03:46Z). The colour is a chip, using the symbols of UX.101; the colour is
never the only signal (the word is in the cell). Markdown and wikitext render the same rows as lines of bold text, as
GEN.216's examples do, because a table is not needed there. One module builds the rows from `physics/habitability` and
`physics/protection_limits` (the single limits table), so the colour, the score, the band the sentence names and the gear
come from one calculation. Until GEN.215 and GEN.216 merge, W2 calls the existing functions and the rows carry the
current wording; W6 replaces the wording.

## 4. Renderers

| Renderer | Where | Output |
|---|---|---|
| HTML | `web/lib/` (from JSON) | A body's detail is: a facts table (`<table>` in the site's `table-scroll` wrapper), the ratings table with chips and symbols, the composition bars, prose paragraphs, then its moons. Quantities are `<span class="qty" data-q="..." data-si="...">`. Same tables, chips and both themes as the rest of the site; no header levels or Markdown, so no `mdrender`. |
| Markdown | `describe/markdown.py` | Today's layout: `# system`, `## planet`, `### moon`, `Property | Value` pipe tables, bold rating lines. Quantities through the text ladder; superscripts as `<sup>` as now. |
| wikitext | `describe/wikitext.py` | Today's layout: `= system =`, `{{Planet Data|...}}`, `{{Class Data|...}}` for moons, `{{Exp|a|b}}`, `[[Category:Star Systems]]`. |
| text | `describe/text.py` | Plain tables and sentences for the terminal (Unicode superscripts, no markup). The command line default. |
| JSON | the API | The Document as above. |

The web page keeps its row list (name, chips, a `<details>` per body); only what is inside a row changes, from converted
Markdown to the Document's HTML. The page-level Wikitext/Markdown buttons keep working: they ask `/text?format=` and show
the result in the code box (and the wiki upload posts the same text). They are on request only.

## 5. Inputs

`describe(system)` takes a `StarSystem` object graph. That is what `db/store.load_star_system` already returns for a stored
system and what generation builds for a new one, so the stored page, a fresh `planetgen system --output -` and the
one-off page all go through the same code with no second row shape. Each body keeps its `db_id` (set by the loader) so a
Section's `ref` matches the row the page links to. The habitability rows read the same column-named values the stored
scores come from (`habitability_world.values_from_body`).

Cost: the page already loads the graph for `/sections`; the Document replaces that call, so there is no extra load. The
Document is not cached in the database; the page cache (`web/lib/pagecache.py`) already covers repeat views.

## 6. Wording rules (the style the sentences are written to)

- Short, plain, in-universe: the voice of a chart or a survey, not of the generator. Public pages use the UX.42 words
  (chart, charted, uncharted); admin pages keep plain wording. Boss's lines are the model: "The air pressure is 0 kPa, so
  a pressure suit is needed."
- One fact per sentence; a sentence names its value (as a quantity) and, where it helps, the band it falls in.
- No filler joins ("of a mix of", "which is expected to", "approximately" before a computed value). Lists use "A, B and C".
- Ages and durations use the ladder (Gyr, Myr after UX.81); "Billion Years" goes.
- Moons and planets share templates with a noun parameter; no regex rewriting of planet wording into moon wording.
- "Habitable" in a sentence means the PHI-6 rating (overall Green or Blue is "potentially habitable", default, Boss may move
  the line). The class flag is called "habitable class". The zone keeps its name ("habitable zone").
- The index is PHI-6 everywhere a person reads it. Stored column names (`phi4`, `phi4_<domain>`) are internal and stay
  (default; renaming is an optional migration).

W10 rewrites the fixed text (class descriptions, star notes, flavour) to these rules; this note does not.

## 7. Stored prose becomes structure (needs a migration)

| Stored now | Becomes | Why |
|---|---|---|
| `planet_evolutionary_paragraphs`, `moon_evolutionary_paragraphs` (finished paragraphs with ages already formatted) | `life_pace`, `life_stage` (a `MILESTONE_KEYS` key) and the milestone ages, or just the pace and stage, with the ages recomputed from the star | The ages follow the unit preference; the stage no longer has to be parsed back out of text (`evolution.life_stage_from_paragraphs`); wording can change |
| `flavor_text` | a flavour key into `tuning.*_FLAVOR` | Rewording needs no regeneration; "Sensors show ..." stays in `describe/` |
| `hab_note` (VARCHAR 255) | removed; the chip tooltip is built from the stored scores | It is derivable |
| run-config `markdown` flag | removed | A presentation choice made per request |

Generation draws are unchanged, so galaxies stay reproducible; only what is stored for the text changes. Only Foundations
Database adds migrations, so this is a request to that lane (work-list item W12); the describer reads the structured
fields when they exist. Until W12, W2 reads the stored paragraphs and `flavor_text` as they are and carries them as
`Prose` blocks. Because there is no backward compatibility to keep, W12 changes the schema in one step and deletes the
old columns, tables and the parsing helpers with it.

## 8. API and command line

- `GET /api/systems/<id>/description` returns the Document (JSON). `GET /api/systems/<id>/text?format=markdown|wikitext|text`
  renders it. `GET /api/systems/<id>/sections` and `render_system_sections` are removed with `mdrender.py` and the
  `markdown` dependency of the web app. `docs/api.md` documents both.
- CLI: `--markdown` on `system`, `sector`, `galaxy` and the phenomenon commands becomes `--format text|markdown|wikitext`
  (default `text`; Decided default). `--output FILE` writes whichever format was asked. `-` is stdout.
- The wiki upload (`wiki/`) renders wikitext or Markdown from the Document when an admin uploads, as now.

## 9. Tests

- Coverage: every fact key and sentence topic of the old render has a home in the Document (a fixture over the example
  systems in `examples/systems/`).
- Renderer parity: Markdown and wikitext for the same Document carry the same facts and the same sentences; quantities
  agree with the ladder text (`util/ladders.py`) and the browser's `unitladder.js` (existing Python/JS parity test).
- No-format rule: nothing in `describe/model` or the builders contains Markdown, wikitext or HTML (a scan test), and no raw
  numeric formatting (`:.Nf`, `toFixed`) in `describe/`.
- Web: the system page contains no output of `markdown` and no `<sup>`; every quantity is a `.qty` span with `data-si`.
- Index and wording agree (W9): colour, score, band named, gear and equipment label from the same world, for random worlds
  and the example systems; an airless body is never better than Red in Air; text that claims life agrees with the life
  guard.
- The Markdown tests (`test_mdrender.py`, `test_fuzz_text_rendering.py`, `test_system_render.py`,
  `test_bughunt_render_edges.py`) are rewritten against the Document and the renderers, not patched.

## 10. Order of work

1. This note (W1).
2. `describe/` builds the Document for every section kind from a `StarSystem`, with the coverage test; the old
   renderers still exist (W2). Quantities not yet on a ladder (mass, temperature, pressure: UX.22 and UX.23, Foundations UX)
   are carried as a plain number with a unit symbol until their ladders land, then switched by changing one table.
3. The three consumers move over and the old code is deleted in the same PR as each move: the web page and one-off page (W3,
   which is UX.30), the API and exporters with the removal of `SystemConfig.MARKDOWN` (W4), the CLI (W5).
4. Wording (W6, W8, W10) and units at the remaining call sites (W11), after GEN.215 and GEN.216 for the habitability text.
5. The stored-prose migration (W12), when Foundations Database has room.

## Open questions (defaults apply unless Boss says otherwise)

1. The command line prints plain text by default; Markdown and wikitext by `--format`. Default: yes.
2. "Potentially habitable" in the overview means overall Green or Blue. Default: yes.
3. The one-off system page shows the HTML page with a download button per format, not a Markdown code box. Default: yes.
4. Keep `phi4*` column names. Default: yes.
