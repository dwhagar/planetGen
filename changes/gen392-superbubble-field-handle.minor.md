**Headline:** The galaxy shape now carries a superbubble handle: a bubble's effect on the young stars (0.6 inside a full-age cavity, a derived shell of 13.1 times) is built and bounded, ready for the density terms to read.

### Added

- GEN.392: `GalaxyShape.bubbles` holds a `BubbleField` (`galaxy/bubble_field.py`), the galaxy's seed, with a box cut of the bubbles, the young-star factor at a point, a certified maximum over a box and a bubble-free density read. A bubble's age comes from its radius by Weaver's law at one driving power, capped at the 40 Myr the driving power lasts; the shell is 0.06 of the radius thick with star formation at its density to the power 1.5. Stand-ins registered as section 6.339 (SB1 to SB5); design in `docs/design/galaxy/superbubble-density.md` section 8. Nothing reads the factor yet, so the density and every generated galaxy are unchanged until GEN.393.

generation-output: unchanged
