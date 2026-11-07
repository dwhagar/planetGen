**Headline:** The bring-to-epoch help no longer shows an internal item code, and two tests that were red on main are fixed.

### Changed
- **bring-to-epoch help (TEST.684).** The text for `--no-wide-binaries` named an internal item code; it now says only what the option does. The tests that keep item codes out of shown text stay strict.
- **Tests (TEST.684).** The ordered-pair fuzz check no longer reads the white dwarf helium envelope table (each pair is a white dwarf mass and its helium mass, not a low and a high) as a range.

generation-output: unchanged

TODO items resolved: TEST.684.
