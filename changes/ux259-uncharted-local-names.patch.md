**Headline:** An uncharted star or scattered object is now shown by its local name and birth sector, with its real Object ID, instead of "Uncharted G2V star 123".

### Changed
- **Uncharted objects carry their local name and Object ID on every page and list (UX.259).** The uncharted star page, the uncharted phenomenon page, the Uncharted stars table, the sector page rows, the Nearby list and the scattered-phenomena lists now name a scatter row by its local name (with its sector where the name stands alone) and show its Object ID, read from the row's own `uid` (GEN.1254, GEN.1255). A row made before the Object IDs existed, or a galaxy with no naming key, keeps the older label. Display only: nothing stored changes. Map labels and the map infobox are Lane 12's side and follow once the API fields land.

generation-output: unchanged
