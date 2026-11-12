**Headline:** An edited world inside a cloud keeps the cosmic-ray boost of the squeezed heliosphere.

### Changed
- **GEN.1607 confirmed not a bug; a test now holds it.** `radiation.generate_dose` does set `dose_helio_mult` back to 1.0 each call, but the only paths that re-apply life or air to a stored body are the admin edits, and `db.edits.save_system_edits` ends in `store.refresh_system_dose`, which works the multiplier and the surface dose out again from the container the system is stored in. `test_an_edited_world_in_a_cloud_keeps_the_cosmic_ray_multiplier` places a system in a dense cloud (multiplier 2.5), re-applies life and the dose to each grounded planet (multiplier 1.0), saves the edit and checks the stored multiplier is 2.5 again.

generation-affects: none
