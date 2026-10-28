**Headline:** A rocky body's stored heat flow and Urey ratio are now those at its star's final age, not the provisional younger age they were first made at.

### Changed
- **The heat is run again once the star's age is settled (GEN.763).** A star's age only rises after its planets and moons are made, but each body stored its heat flow and Urey ratio at the provisional age. Many bodies therefore showed a Urey ratio of exactly 1.0 (the mantle held at its ceiling) and a small moon showed 4.8. The life pass now runs `interior_heat.update_heat` on every planet and moon at the final age (no draws, so nothing else moves). Over 24 sampled systems the medians run from 0.78 (under 1e-4 Earth masses) to 0.49 (2 to 5 Earth masses), and one body ends above 1.0, a small moon whose mantle is still warming.
- **The design note covers small and large bodies (GEN.763).** Section 3.2 of `docs/design/composition/interior-heat.md` has the sweep from 0.0004 to 10 Earth masses, why a ratio above 1 appears only while a mantle warms, and the bounds a test now holds.
- Generator epoch bump (numbered at merge): the stored `urey_ratio` and `internal_heat_flux_w_m2` of most rocky bodies change.

generation-output: changed
TODO items resolved: GEN.763. Added: none.
