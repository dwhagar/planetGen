**Headline:** A sector can now hold an asteroid field left by two colliding rogue planets, at the physical rate: about one field for every 1e20 cubic parsecs.

### Changed
- **The density of asteroid fields leaves 0 (GEN.1103).** The sector roll's density of the kind was 0, so no field was ever rolled. It is now the density of the fields that two colliding rogue planets leave (`generation/collision_field.py`): the rate of the collisions that fall in the rubble window of the pair's speeds (`field_rates`), times how long the pieces take to pass a light year (`field_dispersal`), about 7e-21 standing per cubic parsec. A field a sector rolls is one of those: a pair of terrestrial rogues, a speed in the window, an impact that is not grazing, the layers and the ledger of what flies free (`field_composition`, `field_mass`), with the makeup of the layers, the speed of the pieces and a history sentence of the impact (its peak shock and the second remnant). It replaces the hand-made draw for the fields the roll makes. Fields of supernova remnants are hosted by their remnants (GEN.1104) and not rolled, and the Rogues top-up adds its own. Physical sources only, no game rate (Boss, 4:54 PM EDT 2026-10-30), so the roll all but never fires; the builder is tested directly.

generation-output: unchanged

TODO items resolved: GEN.1103.
