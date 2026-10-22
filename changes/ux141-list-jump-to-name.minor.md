**Headline:** The list pages for stars, planets, moons, belts, comets and phenomena jump to a name.

### Changed
- **Start a list at a name (UX.141, third part).** Under the name sort, the "Starting at" box on each `/lists/<kind>` page begins the list at the first name from the typed text on, like the Systems list, and the page after it still starts at the last object of the one before. `GET /api/object-lists/<kind>` takes `from=<name>` for the same. A start past the last name shows the last page with a note.

TODO items resolved / Added: none / Added: none
