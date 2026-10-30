**Headline:** A search with a stray life class value no longer fails, and the search tests match the new filters.

### Fixed
- **Stray life class in a search.** `/api/search?life=:` (a value that is no life class) returned a server error; it now matches nothing, as a stray value of the other filters does. The search tests also know the new "supports a life class" filters and the reference every planet row carries.
