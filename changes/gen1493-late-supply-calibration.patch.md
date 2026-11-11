**Headline:** The late supply a rocky planet takes from the belts of its own system is built and calibrated, but nothing reads it yet, so nothing generated changes.

### Changed

- New `physics/late_supply.py` (GEN.1493 slice 0): the share of a belt's early mass that lands on a planet, with the other planets, the giants above all, as the sinks (capture cross-section with focusing), reusing `influx.reach`. Stand-in rows LS1 to LS4 (section 6.254), burn-down lines for LS1 to LS3, a wiring-register line for the module.
- A census of 581 generated rocky planets (`tests/audit/late_supply_census.py`) puts the fitted inward share at 0.105 (0.1 kept) and shows that planets of systems with a giant take about 20 times less, that 31 percent take nothing (no belt outside them), and that the spread is much wider than the present law's. These go to slice 1.
- Nothing is wired, stored or drawn; output is unchanged.

generation-output: unchanged
generation-affects: none
