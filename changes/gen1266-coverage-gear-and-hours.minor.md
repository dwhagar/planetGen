**Headline:** The gear a world asks for and the hours outside are now those of its best site, with a note naming the worst conditions elsewhere.

### Changed

- The required suit and body protection are chosen for the best site (the place livable the largest share of the year) at its worst hour, and the heavier gear that the worst place and hour would need is recommended with where and when (GEN.1266). The Protection text adds a sentence saying so, and "livable only in the terminator belt" for a world that turns once an orbit.
- The hours a person can be outside in a day are worked from the best site's own coldest and hottest hour through the year; a locked world uses the real range of its terminator belt. The API gains a `temperature_coverage` object under the equipment block (livable share, area and time, the best site, the coldest and hottest place) and the property sheet gains the matching rows.
- Additive: no API version bump.
