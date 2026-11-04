**Headline:** A generated system now keeps the plan its planets were built from, so the bodies it threw out can be read from it: the first step of making rogue planets come from real ejections.

### Changed
- **The architecture of each star is kept (GEN.1111).** `StarSystem` holds, for each star, the planet architecture its kept planets were built from (the last attempt's), and `ejected_losses()` gives the bodies that architecture threw out (the `thrown` losses of GEN.1119 that end ejected, with the scatterer, the distance, the formation distance and the Hill radius). It is not stored: a system rebuilt from the database has none, because a donor's rogues are made and saved when it is generated. Nothing reads it yet; the ejection events and the rogue rows of GEN.1111 do.

generation-output: unchanged

TODO items resolved: none (GEN.1111 stays open).
