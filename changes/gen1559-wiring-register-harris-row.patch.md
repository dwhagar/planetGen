**Headline:** The wiring register marks the Harris cluster catalogue as a reference that is never loaded into a galaxy.

### Changed

- GEN.1559: the `harris_clusters` row of `docs/design/research/wiring-register.md` changes from an open item to `by design`, with the reason (Boss, 2026-10-16: its distributions guide the laws and no real cluster is placed in any galaxy; the cluster generator's laws were calibrated against it and are typed in). The other rows that wait on a feature (`dynamical_state`, `formed_at`, `origin_event_id`, `object_uid.is_on_demand`, `region_run.next_units`) keep the open item each waits on, so they come back when it is built.

generation-output: unchanged
generation-affects: none
