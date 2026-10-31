**Headline:** The orbital event writer exists: a collision found between two bodies can now be turned into a stored event, with the object-ID rules.

### Changed
- **The orbital event writer (GEN.197).** `planetgen/db/orbital_events.py` asks the pair closest-approach library whether two bodies' stored paths come within the sum of their radii in the update's window (`find_contact`), resolves the hit with the collision outcome table (`resolve_contact`) and writes it as one `orbital_events` row with its `orbital_event_bodies` rows (`record_event`): the event kind, the numbers it was decided from, a mass and momentum check worked out from the outcome, the participants and the products, and one line in the activity log. The ID rules are in it: a merger keeps the heavier body's ID and retires the other, an ejected planet keeps its ID, the field a collision makes takes a new run-time ID, and a destroyed body's ID is never handed out again. No migration; nothing calls it yet, because the orbital update that finds the pairs is a later item, so generation output does not change.

generation-output: unchanged
