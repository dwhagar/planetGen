**Headline:** A new table, the effect ledger, records for every body of a system and every effect that acts on it (light, UV, tide, heating, Kozai and secular forcing, precession, eccentricity perturbation) how large each source's contribution is compared with the largest source. It changes no result.

### Added

- `effect_ledger` table (migration 0236): a row for each source of an effect on a body, its share of the largest source (1.0 for the largest). Every source is computed; a light or UV source under a share of 1e-3, 0.1 percent (`tuning.EFFECT_LEDGER_FLOOR`) is added to a summary row of the body and effect, with how many there were and the largest source's size, so the tail is never lost and a floor under 1e-37 stores everything. The five gravity effects keep only the largest source and the summary row and are recomputed on demand.
- Generation fills the ledger for light and far ultraviolet: after the bodies are final, each planet and moon notes the bolometric flux and the far-ultraviolet ratio each star of its system puts on it (`radiation_field.field_at`, GEN.1441), by source. Nothing is drawn and no value changes.
- `system.effect_ledger` takes `add(effect, receiver, source, value)` from the code that computes an effect; `store.insert_star_system` writes it with the system, and deleting the system deletes its rows.

generation-output: unchanged
generation-affects: none
