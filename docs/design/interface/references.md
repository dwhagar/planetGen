# Help > References: the out-of-universe section (DOC.87 to DOC.94)

Status: R1 (DOC.87) and R2 (DOC.88) built; R3 to R8 planned. Boss asked for this on 2026-10-31 (2:49 PM EDT) and approved the design.

## The idea

The whole site speaks in universe: a working navigation system that simply knows things. One Help entry, **References**, is the only
place that speaks out of universe: where each number came from, what we have research for and what we do not, every stand-in (the GEN.725
registry) and how it was decided, what was fitted to a result, and how scientific each part is. Every other page keeps only a quiet
"Where this comes from" link to the exact References entry.

## What is built

* **Addresses and links (DOC.87).** `web/lib/reference_urls.py` is the one place that knows the layout; a page calls
  `reference_link(kind, id)` (`web/helpers.py`, also a template global) and never writes a `/references` address. See "Links into
  References" in [text-style.md](text-style.md).
* **The word rule (DOC.87).** `plainwords.OUT_OF_UNIVERSE` lists the words that tell a reader the site is a model; they appear only on the
  References pages and the admin pages (`exempt(path)`). `check_sentence(sentence, references=True)` reads a sentence of a References page,
  where those words are allowed and every other rule still holds. `tests/test_references_crawl.py` crawls the pages and fails on a word
  outside References that its `KNOWN` list does not name under the item that removes it.
* **The shell (DOC.88).** Help > References is the last entry of the Help menu, after a divider. `/references` is a static page (no
  database, no script) that opens with the Out of universe banner (`id="out-of-universe-banner"`, the contract every References page
  keeps) while the site header takes another tint (`body.references-mode`, `--references-bg` and `--references-border`). Its numbers are
  read from `docs/references/overview.json`, which `scripts/build_reference_pages.py` builds from the reference list
  (`docs/reference/references-data.json`) and the stand-in registry (`gap-filling-plan.md`, read by `build_stand_in_values.py`'s parser).
  `--check` fails when the committed file differs, and `tests/test_references_pages.py` holds that; `update.sh` runs the build in step 3 and
  a failure only warns. Later items add their files to `BUILDERS` in the script and their pages to `web/references_pages.py`.

## The page tree

1. **Overview** `/references` (built): what this section is, the headline numbers computed from the registry, the tier and tweak legends.
2. **Papers and books** `/references/papers` (DOC.89): by letter, area and field, with search; each entry links to its DOI page, else arXiv, else ADS, else a search; back-links to where it is used.
3. **Datasets, catalogs and archives** `/references/data` (DOC.89).
4. **Stand-in registry** `/references/stand-ins` and one page for each stand-in (DOC.90), and **Fitted to a result** `/references/fitted`.
5. **Method by area** `/references/method/<area>`, **what we have and lack** `/references/gaps` and the computed **scorecard** `/references/scorecard` (DOC.91); `/help/properties` moves here with no redirect.
6. **Speculative chemistry** `/references/speculative` (DOC.92), with one section for each life class and one for a world with several living classes (GEN.1209, GEN.1213).
7. **Credits and version** `/references/about`.

## Rules

* The in-universe pages never say speculative, stand-in, tweak, fitted, tier, model, seed, generated, version and the like (rule 16 of the
  [plain-words guide](text-style.md)). A survey-style word ("estimated distance", "unsurveyed") stays.
* Every number on a References page is built from a registry, never typed.
* The API field `speculative`, the developer docs, the Debug view and the changelog keep their words.
* DOC.93 strips the leftovers (color note, sector text, life class "reaches", sheet source lines to links); DOC.94 locks the crawl once
  DOC.93 lands.

## Citations (DOC.97 to DOC.100)

Boss (2026-10-31, 7:22 PM EDT): all citations and references live in the central registry (`docs/reference/references-data.json`) and `docs/REFERENCES.md`, nowhere else; and a reader of Help > References can click a citation to reach the reference list and the document information, in that section only.

* **DOC.100 (built)** gives every record a stable citation ID and an anchor on /references/papers (`#paper-<id>`). The ID is `R` and at least four digits (`R0123`). `scripts/build_references.py` (`assign_ids`) gives each paper, book, dataset, archive and setup dataset its ID once, in the order the list prints on the first build, and writes it to the registry (`ids`, `id_next`, `id_aliases` in `docs/reference/references-data.json`). An ID is never reused, even when its source is removed. A source goes by every name it has had (its DOI, its arXiv id, its key), so a preprint that gets a DOI keeps its ID, and two IDs found to be one paper leave the older as the ID and the newer as an alias that still resolves. `docs/REFERENCES.md` and `.html` print the ID in bold before each entry, with the anchor `paper-<id>` (`R0123` is `#paper-r0123`). `web/lib/reference_urls.py` (`is_paper_id`, `paper_anchor`, `reference_path("papers", "R0123")`) is the one place the site learns the shape.

### The in-line citation style (Boss, 7:27 PM EDT 2026-10-31: "use whatever in line citation style is standard for NASA documentation")

NASA technical documents (NASA Technical Memoranda and standards, and the NASA Systems Engineering Handbook) cite by number in square brackets, `[12]`, and list the numbered references at the end. This project does the same with one change: the number is the source's own fixed ID, not its order in the document, so a citation survives moving the sentence and one search finds every place a source is used.

* **One source:** `[R0123]`, after the words it supports and before the closing mark: "the slope is 1.3 [R0123]."
* **Several sources in one place:** one bracket, IDs in ascending order, a comma and a space between: `[R0123, R0456]`. No ranges, no "and".
* **A place inside the source** (a table, a page, a section) follows its ID after a colon: `[R0123: Table 2]`. The page numbers the notes carried stay with it where they were printed; they are not required.
* **Where:** the same in a design document, a report, the TODO text, a code comment, a docstring and a data-file header.
* **No author, no year, no journal, no DOI, no arXiv number** appears outside the registry, `docs/REFERENCES.md` and the References pages. Whatever the old text said about the source ("read", "recalled", "Table 2 of", the value taken) goes to the registry entry of the stand-in it supports, or stays as plain words about the value with no name in them.
* A page a visitor reads has no citation at all. It keeps only the quiet "Where this comes from" link.

* **DOC.97** makes each citation written on a References page a link to that anchor, and fails the build on an unknown or ambiguous one.
* **DOC.98 (built for the design documents, scripts, source code, data files and tests; the TODO text, plan and generated pages follow through the TODO thread)** moves citations out of everything else. `scripts/rewrite_citations.py` turns each author and year, DOI and arXiv number into the source's ID, using the ties the settle step recorded (DOC.95): `report` counts, `apply` rewrites, `check` fails while one is left. In Python it rewrites comments and docstrings (the syntax tree minus docstrings must be unchanged), and string literals and data files with `--strings`; the files built from them (the data files of `scripts/build_*.py`, `stand-in-values.md`, the dials block of `docs/reference/cli.md`) are rebuilt after. Where two papers share a first author and a year and nothing beside the mention chooses (volume and page, the letter `2023a`, an identifier, title words in the sentence), `--all-candidates` cites both: those places are the ones for Sources 1 to narrow. The journal, volume and page that followed a citation go with it; names written without a year stay (they cite nothing).
  * **To cite a new source:** `python scripts/build_references.py add <doi or arXiv number>` looks it up (Crossref, DataCite, arXiv), lists it and prints its ID; then `python scripts/build_references.py render`. A source those services do not know is refused: it is not listed, so it is not cited. Do not run `harvest` or `match` any more: they read citations out of the text, which the text no longer holds, and would drop the record of how each old citation was tied.
  * **Rebasing a branch across the sweep:** merge main, take the sweep's side for any conflict that is only a citation, then run `python scripts/rewrite_citations.py apply --all-candidates --strings --scope <your changed files>` and rebuild the generated files you touch.
* **DOC.99** locks it: a citation written outside the registry fails the gate for changed files.

