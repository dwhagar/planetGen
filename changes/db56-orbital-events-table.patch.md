**Headline:** The database now has a place to record every collision, merger, split, ejection, capture, unbind and warning of the orbital update.

### Added
- Migration 0118 (DB.56) adds `orbital_events`, one row per orbital event with its kind, outcome code, regime, impact and detection time, the seed string that replays it, the numbers it was decided from, the conservation check, the products and the flags of docs/design/collisions-and-mergers.md section 7 (including `numeric_warning`), and `orbital_event_bodies`, the participants and products of an event, one row each. Nothing writes the tables yet (GEN.197's `record_event` does) and no view reads them (GEN.143).
