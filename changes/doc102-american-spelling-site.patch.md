**Headline:** The rest of the site and the command-line help now use American spelling too ("color", "center", "neighbors", "canceled", "disk").

### Changed
- **American spelling everywhere a visitor reads.** Command-line help, progress labels, error and status messages, the admin settings, the Generate page, the sector page, the galaxy settings help, and the text of the star, system, sector, nebula, and phenomenon pages no longer use British spellings. About 150 sentences changed.
- **A check keeps it that way.** A new test reads every sentence the program shows and every template and page for British spelling, and the checker's word list is longer ("canceled", "catalog", "disk", "toward", "modeled", "labeling", and more).
- Stored values and names that carry a British spelling ("cancelled" as a job state, "vapour" as a hydrosphere value, "centre" as a key in the API) are not changed here, because renaming them changes stored data or the API.
