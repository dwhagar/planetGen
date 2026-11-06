**Headline:** A rocky world with a hydrogen-rich air stores the helium and heavy-element shares of that air as its envelope's details, as a giant does, so the sheet and the columns no longer show an envelope with nothing in it.

### Changed

- **The envelope of a hydrogen-air dwarf has its details (GEN.1303).** `envelope_share.fill` (GEN.1226) stores a rocky planet's or moon's envelope share from its hydrogen-rich air, but left `envelope_helium_fraction` and `envelope_z_fraction` empty, so five bodies of the audit sweep (envelope shares of 2e-8 to 2e-5) failed `detail-ranges` of `test_sweep_bulk` ("envelope_helium is None with an envelope share of ..."). `envelope_share.air_details` now works both out from the stored air: the helium as a share of the air's mass and the elements heavier than helium (counted by element, so the oxygen of water is heavy and its hydrogen is not) as the other share; they are `None` again when the air stops being hydrogen-rich. A body whose air holds no helium stores a helium share of 0, which is what its air holds. Only bodies with a hydrogen-rich air change.

generation-output: changed
