**Headline:** A frozen world's page now says how much high-pressure ice lies under its water, including none.

### Changed
- **GEN.781: the high-pressure ice thickness is stored for every body whose water is ice or an ocean frozen over.** The water column already worked it out, but a column with none was stored as empty, so the sheet read "not stored" for the thousands of frozen worlds that have none. A body in the state `ice` or `ice-covered ocean` now stores 0 km where no ice VI or VII lies under its water, and the thickness where it does; an open ocean, a dry body and a vapour world still store it only when there is some, as the sheet rule says. The fill audit no longer excuses the column.

generation-output: changed
