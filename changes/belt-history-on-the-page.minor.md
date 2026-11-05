**Headline:** a belt now shows what its star's evolution did to it, in plain words on its page.

### Added
- Schema v200 (migration 0200, GEN.1232): `asteroid_belts.history_text`, one paragraph of the sentences `physics/belt_history.py` wrote for the events a belt went through (the part the giant branch engulfed, the share the mass loss threw inward). A belt nothing happened to stores nothing, and so does a row made before v200, which gets it when the galaxy is regenerated.
- The belt's page shows that paragraph after its composition, and its property sheet has a "What the star did to it" row.
- `belt_history.history_text` is now called by the store, so its wiring register line is gone. The giant event and companion limit (`giant_event`, GEN.1232b) are still not called by generation; that waits for the belt birth.
