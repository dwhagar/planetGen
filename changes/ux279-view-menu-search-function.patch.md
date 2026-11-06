**Headline:** The View menu's search is written and tested, ready for the search box.

### Added
- **Search over the View menu's entries (UX.279).** `static/viewmenusearch.js` takes the layers' entries and the text typed and says which match
  (every typed word must start a word of the name, the group's name or a keyword, in any case and order), in the menu's order, with the
  headings and groups to open and a screen-reader count ("5 matches"). Nothing shows it yet: the search box is UX.281. No look change.
