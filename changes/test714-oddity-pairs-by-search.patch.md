**Headline:** Two oddity tests and the brown dwarf test now find their case by searching instead of naming the draw that had it.

### Changed
- **The dark-side air loss test (W19) and the second climate branch test (W20) search draw keys (TEST.714).** Each named one key (`w19-M3V-2`, `w20-G2V-10`) for a planet in a rare state: a locked M3V planet whose thin air the dark side takes, a G2V planet whose climate has a second branch. The W20 key was already replaced twice when the draws moved. Each now takes the first key of a bounded range that has the case and keeps its assertions, and names how many keys it tried when none does.
- **The brown dwarf system helper in the dark light test keeps its six seeds and searches 200 more** if none of them draws a brown dwarf any more.

generation-affects: none
