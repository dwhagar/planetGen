### Added

- The Species list, the Polities list and a polity's Systems are data tables
  like the others (UX.41): sort from the headers, filter species by
  spacefaring and era and polities by government and era, and scroll through
  every row. The old "All / Spacefaring / Not spacefaring" links are now the
  Spacefaring menu (`?spacefaring=yes|no`).
- `GET /api/species`, `GET /api/polities` and `GET /api/polities/<id>` take
  `sort`, `order`, filters and `facets=1`.
