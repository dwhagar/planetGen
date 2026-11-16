**Headline:** The wiring register no longer lists five built-and-unread items as waiting; each now says why nothing in the app reads it yet.

### Changed

- GEN.1559: five rows of `docs/design/research/wiring-register.md` change from an open item to `by design`, each with the reason: the Harris cluster catalogue (`harris_clusters`; Boss, 2026-10-16: a guide to the laws, never loaded into a galaxy), the `dynamical_state` column of planets and moons (every generated body is bound; the unbind and capture test of GEN.291, release 8.3, is the only writer of `unbound`), `asteroid_fields.formed_at` and `origin_event_id` (empty for every generated field; written by the orbital events of GEN.346, release 8.3), `object_uid.is_on_demand` (the rule the belt member page already applies by looking the rank up; the field member ledger of GEN.1105 reads it) and `region_run.next_units` (a region run has no stage 3 yet; the command reads it when GEN.1120 turns it on).

generation-output: unchanged
generation-affects: none
