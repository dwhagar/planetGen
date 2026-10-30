**Headline:** The search can now find systems by how many stars they hold: a Stars menu on the search page, `stars:3` in the search box, and `stars=` in the API.

### Added
- **The stars filter (UX.239, first slice).** On the search page a Stars menu (Single, Binary, Trinary, Quaternary, Quintenary and beyond) keeps the systems with that many stars; the Systems panel lists them with their star count, whether or not a name is typed. In the search box, `stars:3` finds triples, `stars:3+` three or more, and `stars:binary`, `trinary`, `quaternary` or `quintenary` work too (the rest of the box stays a name search). `GET /api/search?stars=` takes the same values (a count, a count and a plus, or a name) and answers `400` for anything else. The count is read from the stored stars of each system. Nearby, the command line and the sort key are the later slices.

No generator output changes.

generation-output: unchanged
