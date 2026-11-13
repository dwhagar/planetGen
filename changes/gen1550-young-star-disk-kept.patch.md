**Headline:** A young star that keeps a belt now keeps its own protoplanetary disk instead of losing it.

### Changed

- The cold debris-disk step (`StarSystem._set_debris_disk`, GEN.1071) cleared every disk field of a star before it asked whether the star could hold a debris disk. A protostar or pre-main-sequence star that kept a belt (its belts survive the 10 Myr planet cut) lost its class 0, I or II disk (`pre_main_sequence.draw_disk`, GEN.242) this way. The step now leaves any disk that is not a debris disk alone; a debris disk left by an earlier pass is still cleared when its belt is gone, and a main-sequence star is unchanged.

generation-output: changed
