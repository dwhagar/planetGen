**Headline:** The nine reference pages (papers, stand-ins, method and the rest) are now classed in the indexing matrix, so the drift test passes again.

### Fixed
- `web/seo.py` lists `web.references_data_page`, `references_fitted`, `references_gaps`, `references_method`, `references_papers_page`, `references_scorecard`, `references_speculative_page`, `references_stand_in` and `references_stand_ins` as top pages beside `web.references_page`, so `test_every_web_page_is_classified_on_purpose` (TEST.586) is green. No generated output changes.
