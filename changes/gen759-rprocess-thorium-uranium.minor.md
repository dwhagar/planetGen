**Headline:** The thorium to uranium ratio of a rocky body now follows the r-process history of its host's galaxy instead of a fixed scatter, so a younger star's planets are thorium-richer.

### Changed
- **Uranium from the r-process (GEN.759).** `radiogenic` now reads `rprocess_enrichment` (GEN.663): a body's uranium is its thorium over the ratio the r-process leaves, one production ratio drawn from 1.667 to 2.0 (Frebel 2007) times the enrichment of its host's birth gas in the disc of its galaxy (the host's age and the disc's star formation history). The typical star gives a Th/U of 3.2 to 3.9 in rock of the Solar System's age, around the bulk silicate Earth's 3.83. The old uranium scatter of 0.05 dex (a placeholder) and its constant are deleted. Potassium, a multiple of the uranium, and the heat and ground dose that read them move with it.
- Generator epoch bump (numbered at merge): every rocky body's uranium, potassium, heat and ground dose can change a little.

generation-output: changed
TODO items resolved: GEN.759. Added: none.
