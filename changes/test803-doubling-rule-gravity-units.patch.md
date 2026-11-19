**Headline:** The radiation sweep's doubling rule now gives the table the world's gravity in m/s2 instead of in g.

### Fixed
- TEST.803: `test_doubling_the_infrared_active_gases_never_cools_the_ground` passed `body.gravity` (in g) to the greenhouse table, which takes m/s2, so it asked about a world a tenth as heavy as the one it sweeps. On Train 50's draws it hit a 12 bar carbon dioxide air at 1.19 g and read 353.4 K falling to 336.6 K under twice the infrared-active gases; with the gravity in m/s2 the same air rises from 299 K to 327 K. The rule and its tolerance are unchanged.

generation-output: unchanged
generation-affects: none
