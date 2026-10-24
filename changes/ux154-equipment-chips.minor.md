**Headline:** The planet lists and the system page now show each piece of survival equipment as its own chip, solid when required and outlined when only recommended.

### Changed
- **Equipment chips (UX.154).** On a system page each planet and moon row shows one chip per item of its needs list instead of one chip with the items joined by plus signs. A required item is a solid chip, a recommended one an outlined chip (with "recommended" read aloud for screen readers), each in the colour of the PHI-6 domain the item comes from. The Equipment column of the Search page's planet and moon panels draws the same chips; a search row's API entry gains `equipment_items` (`item`, `severity`, `domain`). A world that needs nothing still shows no equipment chip.
