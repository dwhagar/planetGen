**Headline:** `generate-sectors` declares its write to the rogue planets' interior layers, so `test_command_touches` is green again.

### Fixed

- A rogue planet is a planet and stores its interior layers in `rogue_planet_interior_layers` (migration 0207; the writer came in with the radiation-multiplier fix), but `queue/touches.SECTOR_WRITES` did not list the table, so `test_command_touches` failed with "generate-sectors writes undeclared rogue_planet_interior_layers". The write is declared beside the other rogue-planet tables; the check is not loosened (TEST.812).

generation-output: unchanged
