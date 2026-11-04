**Headline:** A planet's or moon's page and map read the words of its class from one table when they are shown, so a change to the wording reaches every world at once.

### Changed
- **GEN.1212.** `planets` and `moons` lose the `description` column (migration 0200), and generation no longer copies the sentence of a body's class into it; `tuning.PLANET_CLASSES` keeps the composition text and the type of each class but no description. The words are read at view time from `physics/data/planet_class_text.json` (`physics/planet_class_text.py`, UX.252). `GET /api/systems/<id>` gives a planet or moon no `description`, a breaking change that raises the API version to 13; the class words are in the system's description document.

generation-output: changed
api-version: changed

TODO items resolved: GEN.1212.
