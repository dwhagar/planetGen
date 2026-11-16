**Headline:** The air-history test of the life gases no longer depends on how big the first generated world is, so it stays green when a generation change moves that world.

### Fixed

- `test_the_history_of_two_classes_steps_at_each_start_and_at_the_rise` was red on Train 40 after GEN.1709 (the inner-system period spread) changed the first world the test's fixture finds (0.78 Earth masses with an air of 3.1e20 mol before, 3.5 Earth masses and 5.2e21 mol after). Its fixed methane fluxes (1e12 and 5e12 mol/yr) left the first class's methane under `biogenic_oxygen.MINIMUM_PA` on the larger air, so the step at 0.5 Gy read 0. The fluxes now go in proportion to the air the world holds; the assertions are unchanged.

generation-output: unchanged
generation-affects: none
