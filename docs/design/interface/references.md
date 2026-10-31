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
