**Headline:** The test that checks a young red dwarf's escape strips surface water from its planets now looks at the surface water itself, so it passes again on worlds that hold many oceans.

### Changed
- **`test_runaway_history` (TEST.504).** The M3V case read the ocean cover, which is 1.0 on the water-rich planets the delivery law now draws, so a loss of a few thousandths of the planet's mass could not show. It now compares `hydrosphere.draw_surface_water` with and without the loss, under the same bound draw. No generator code changed.

generation-output: unchanged
