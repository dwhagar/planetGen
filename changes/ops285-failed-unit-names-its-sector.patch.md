**Headline:** A unit that fails for good now names the sector, or the id range, it was working on (OPS.285).

### Changed

- The text of a unit's failure says where the unit was working: `(on sector 12/3/0)` for a sector unit (shell, band, wedge, as the seed spells an address) and `(on ids 1,200 to 2,400)` for a batch unit. It is in the held unit's error, the run's error, the `error.log` line and the attempt's `error_text` in the ledger, so a failure names its own address after its stage ended and the stage's item file was deleted. The unit row keeps its arguments as before.
