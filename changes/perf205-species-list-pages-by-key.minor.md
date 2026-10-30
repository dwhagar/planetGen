**Headline:** The Species list opens any page as fast as the first one, whatever the number of species, when sorted by name, era, tech level or spacefaring.

### Changed
- **The Species table pages by key.** The page links carry the species the page before ended on, so a page deep in the list costs what the first does; before, each page read and threw away every species before it. At a million species a deep page by name went from 8 seconds to 4 milliseconds, and the first page by era or tech level from 5 to 13 seconds to 4 milliseconds. `GET /api/species` takes `after=<species id>` for the same.
- **Sorted by name, era, tech level or spacefaring, the list reads in index order** (new indexes, built online by the update). Sorted by homeworld or system it still sorts every species, since those names live in other tables (the polity sort is quick too, as few species have a polity).
- **Within equal values, a descending sort now lists the names in reverse too** (it kept the names in order before).
