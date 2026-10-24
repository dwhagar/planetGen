**Headline:** Every phenomenon page now says what the object is in words, and links to the phenomenon before and after it.

### Changed
- A stored nebula, quasar, rogue planet, interstellar comet, asteroid field, supernova remnant, black hole or neutron star page has an "About" paragraph, built from the stored row by the same sentences the generator's description uses (UX.185). A row too bare to describe shows no paragraph.
- Each phenomenon page has Previous and Next links to the free phenomenon of its type before and after it. The API answer for a phenomenon carries them as `previous` and `next`, and a supernova remnant's centre as `compact_remnant_name`.
- The neighbour links wrap on a narrow screen.
