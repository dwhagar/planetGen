**Headline:** A wet or water-rich world no longer brings nitrogen with all of its water-linked carbon, only with the carbon up to a dry world's amount and a drawn share of the rest.

### Changed
- **GEN.1584 B: wet-class nitrogen (stand-in NB5, L3).** `rocky_inventory` counted nitrogen from all the carbon, and a wet or water-rich class's carbon is tied to its water, so a water-rich body carried several times the bulk silicate Earth's nitrogen. Nitrogen is now counted from the carbon up to the dry class's log-middle (0.056 wt% CO2) plus a drawn share (`nitrogen_with_wet_carbon`, 0.03 to 0.5, log) of what is above it, for the classes b2, b3 and b4; the dry class is untouched. The share is a design call (the Earth took about 2 percent carbonaceous material and its nitrogen mostly from the dry reservoir, [R2923], Johansen and others 2023), hand-chosen against the Earth analog; this slice stands alone on its own commit so that it can be dropped. Before and after in the PR.

generation-output: changed
generation-affects: fill
