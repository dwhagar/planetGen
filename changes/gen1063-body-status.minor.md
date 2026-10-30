**Headline:** A body now has a status worked out from its mass, composition, orbit and neighbours (planet, dwarf planet, minor body or moon), kept with the numbers behind it; nothing in the generator uses it yet.

### Added
- **Body status (GEN.1063).** A new `physics/body_status` module gives a body its IAU word and the numbers: how round it is (the compaction taper shifted by composition, 1e20 kg of ice and 1e21 kg of rock are half way to round), Soter's mu against the mass of the others in its zone, Margot's Pi and Stern and Levison's Lambda. Earth, Mars and Jupiter come out planets, Ceres, Pluto and Eris dwarf planets and a Trojan a minor body, as test targets and never as inputs. The registry has new rows (BR22 to BR24). The belt wiring (GEN.1069 and later) calls it on a later generation train.

generation-output: unchanged
