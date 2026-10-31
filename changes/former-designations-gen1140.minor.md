**Headline:** A renamed or renumbered object now shows what it used to be called: a "Formerly" line on its page, a `former_designations` list in its API answer, and a search by an old name finds it.

### Added
- **Former designations (GEN.1140, part d).** The sector, system and phenomenon pages print a "Formerly" line (each old name or ID with why it stopped and the generator epoch). The answers of a sector, a system (with each star, planet, moon, belt and comet in it) and a phenomenon carry `former_designations`, newest first. The name search of the Sectors, Systems, Stars, Planets and Moons panels also finds an object by a name it used to have (from the start of the old name), and a reference written with a retired ID (`nebula:<ID>`) names the object that held it. Added fields only, so the API version stays as it is.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (TODO closes GEN.1140 when the branch merges)
