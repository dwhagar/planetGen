**Headline:** The density identity test is recorded again for the new star formation history, and passes.

### Changed

- `tests/fixtures/density_identity.json` is recorded again (TEST.631). `test_the_spiral_written_as_components_gives_the_density_it_gave_before` failed in all six cases on main since train 9: GEN.1230 and GEN.1242 (epoch 131) made the thin disc form from an inside-out infall, which moved the ages of each population and so the mean system mass behind the model's units. Checked point by point before the fixture was rewritten: every total density and the halo are the old ones times one number per shape (1.0010 to 1.0014, and 1.0444 for the large four-armed disc), the total mass is unchanged, and only the per-population units and the mix of a place between age groups move. The fixture's `nuclear_disc` is now `nuclear_disk`, the spelling DOC.103 gave the stored name.

generation-output: unchanged
