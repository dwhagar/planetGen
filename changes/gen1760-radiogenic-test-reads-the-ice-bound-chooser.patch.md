**Headline:** A test now picks a body's tide the way the generator does after GEN.1760, so it no longer expects a warm water body to have none.

### Changed
- GEN.1760 follow-up: `test_radiogenic` rebuilt the heat of a stored body with no tide for any icy body (`tidal_heating.is_icy`); the generator now leaves the tide out only for an ice-bound one (`interior_heat.is_ice_bound`), so a warm water body carries the rock mantle's tide. The test reads the same chooser. No generated value changes.
