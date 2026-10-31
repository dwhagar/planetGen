**Headline:** Help now ends with a References page, the one place that says where the numbers come from and how sure each part is, built from the registries and marked as out of the story.

### Added
- **Help > References, `/references` (DOC.88, References R2).** The last entry of the Help menu, after a divider. The start page opens with the Out of universe banner (`id="out-of-universe-banner"`), and the site header takes another tint (`body.references-mode`, new `--references-bg` and `--references-border` tokens in light and dark), so a reader always knows they have left the story. It shows what the section is, the headline shares computed from the stand-in registry (read from a source, unverified, educated guess from physics, chosen for the output, fitted to a result), how many papers, books and data sets are on the reference list and how many have a DOI page or an arXiv copy, and the nine source tiers and five tweak levels with their counts. It is a static page: no database and no script. It is in the sitemap lists (`seo.TOP`).
- **`scripts/build_reference_pages.py`.** Writes `docs/references/overview.json` from `docs/reference/references-data.json` and the stand-in registry of the gap-filling plan (read with `build_stand_in_values.py`'s parser). `--check` fails when the committed file differs; later pages add their files to `BUILDERS`. The same bytes every build.
- **`update.sh` step 3 builds the pages** (`build_reference_pages` in `deploy-common.sh`): on a current checkout it changes nothing; a failure only warns, now and in the closing warnings, and `/references` then says its pages are not built.
- **`web/references_pages.py`, `web/lib/references_site.py`, `templates/references.html`, `docs/design/interface/references.md`** (the design, the page tree, the rules; indexed in `docs/design/README.md`), rows in `architecture.md`.

### Changed
- **`plainwords.check_sentence(sentence, references=True)`** reads a sentence of a References page: the process words (stand-in, tweak, epoch, placeholder) are allowed there and every other rule holds.
- **`HELP_MENU` and `HELP_SECTIONS`** gain References; `base.html` gains a `body_class` block.
- **`tests/test_references_crawl.py`:** the banner test no longer skips (a page exists), and the crawl reads `/references` too.

`tests/test_references_pages.py` (18 tests) holds the build (same bytes, `--check`, counts equal the registry counted a second way), the page (banner, tint, menu order, guide), and the setup step. No generator output changes.

generation-output: unchanged

TODO items resolved: none (DOC.88 closes when this merges; the TODO thread owns the file).
