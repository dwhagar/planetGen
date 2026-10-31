**Headline:** The pages that will explain where the numbers come from have one address scheme, one link helper, and a word check that keeps the doubt words out of the rest of the site.

### Added
- **`web/lib/reference_urls.py`: the addresses of the References pages (DOC.87, References R1).** `reference_path(kind, id)` gives the address of an entry: `/references`, `/references/papers#paper-<doi>`, `/references/data`, `/references/method/<area>#<value>`, `/references/stand-ins/<id>`, `/references/fitted`, `/references/speculative#class-<id>`, `/references/gaps`, `/references/scorecard` and `/references/about`. An id becomes an anchor by `anchor(text)` (lower case, one hyphen for each run of other characters). A kind or an id that does not exist raises `ValueError`, so a typo fails the page instead of leaving a dead link. `AREAS` is the one list of the eleven method areas.
- **`reference_link(kind, id, label)` in `web/helpers.py`, and a template global of the same name.** It returns the quiet "Where this comes from" link (`class="reference-link"`) and is the only way a page points at References; the layout lives in one module.
- **`plainwords.OUT_OF_UNIVERSE`, `out_of_universe(text)` and `exempt(path)` (rule 16 of the plain-words guide).** The words that tell a reader the site is a model (speculative, stand-in, tweak, fitted, tuned, educated guess, tier, generator, generated, seed, simulated, model, version, and their forms) belong only on the References pages and the admin pages. "Level" is not on the list; "sea level" is in universe.
- **`tests/test_references_crawl.py`: the page crawl.** It renders the Help pages and the class pages, reads the Document of ten seeded systems, every sheet label and every sheet source line, and fails on an out-of-universe word that its `KNOWN` list does not name. The words that exist today are listed under the item that removes them: the properties page (DOC.91), the life help pages (DOC.92) and the sheet source lines (DOC.93); an entry that no longer happens fails, so its owner deletes it. It also fails when a template or module writes a `/references` address by hand, when a link into `/references` on a crawled page does not resolve, and (from the first References page, DOC.88) when a References page lacks the Out of universe banner.

### Changed
- **`docs/design/interface/text-style.md`:** rule 16 now names the list, the exemption, the banner and the link helper, and a new section "Links into References" gives the address layout. `architecture.md` lists the new module.

No page text changes and no generator output changes: the check only reads text, and nothing calls `reference_link` yet (DOC.88 builds the pages, DOC.92 and DOC.93 and the rewrite items add the links).

generation-output: unchanged

TODO items resolved: none (DOC.87 closes when this merges; the TODO thread owns the file).
