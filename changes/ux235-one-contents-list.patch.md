**Headline:** A page that already has its own contents list no longer shows a second "On this page" bar.

### Changed
- **One contents list per page (UX.235).** The automatic "On this page" bar now skips a page that already has a contents list of its own (a navigation block with three or more links to sections on the same page, or any element flagged `data-own-toc`). Of the site's pages, only the admin Settings page showed two (its group list and the bar); it keeps its group list. Pages with none keep the bar.
