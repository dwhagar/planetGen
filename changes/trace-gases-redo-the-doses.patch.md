**Headline:** The ultraviolet index and the radiation doses of a world now follow its air after the trace gases are added to it.

### Changed
- **GEN.524.** `generation/life.py` gives the trace gases of table m to the air and, when they changed it, makes the doses and the stored ultraviolet index again (`radiation.generate_dose` on the stream `_rebuild` takes), so they read the air the body stores; before, they read the air one step earlier, which differed by up to a part in a thousand.

generation-output: changed

TODO items advanced: GEN.524.
