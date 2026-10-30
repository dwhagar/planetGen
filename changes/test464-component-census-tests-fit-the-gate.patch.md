**Headline:** The component pass test that took over 900 s on a four-core machine now takes about 35 s, and the register names the open item for two more unread things.

### Fixed

- **The same-seed component test draws a smaller Cloud (TEST.464).** The test drew a satellite galaxy of 2e8 solar masses (138,000 stars, 150 s a draw) five times; it now draws one of 2e7 (10,000 stars), so the same-seed, order-free and per-component-stream checks run in about 35 s instead of over 900 s. The row test draws 5e7, the least at which the gas-rich Cloud places young massive stars. Nothing the generator draws changes.
- **The whole-sphere component census fills 30 sectors of each region, not 100 (TEST.464).** The rates are dialed up, so thirty show every kind the disk makes; the fixture no longer runs past the 900 s setup limit.
- **The wiring register points `noaa_proton_events` at GEN.805 and `ground_albedo` at GEN.684,** the open items that will read them, now that GEN.747 and GEN.462 are closed.
