**Headline:** The brown dwarfs, neutron stars and black holes of a system's orbit tree now give light to the planets beside them, and a page's sky names them when they are bright enough.

### Changed

- The dark members of the orbit tree are sources of the radiation field (GEN.1441) like any star (`physics/dark_light.py`, `generation/system.py`): a brown dwarf gives the light its cooling gives at the system's age (`brown_dwarf.state`); a neutron star gives its thermal glow from a cooling curve of its age since the collapse; a black hole gives the light of an accretion disk fed by the winds of the stars of its tree (the Bondi-Hoyle-Lyttleton capture, turned to light at 10 percent and never past the Eddington limit), and nothing when no wind reaches it. The temperature, the light and the sky text of a planet near one count it.
- The wind law that sizes a star's bubble now lives in one place, `physics/stellar_wind.py`, and the heliosphere and the black hole's feed both read it; the heliosphere numbers are unchanged.
- A dark member gives bolometric light only: it adds no XUV, flare or cumulative XUV term, and it does not move the truncation of a disk (gravity is outside this item). Systems with no dark member that glows give the numbers they gave before.
- Stand-ins DL1 to DL7 and a new SKY3 are registered in section 6.259 of the gap plan. The effect ledger does not store the dark sources yet, and a stripped helium star and a Roche-lobe overflow feed are not modelled; both are filed as follow-ups.

generation-output: changed
