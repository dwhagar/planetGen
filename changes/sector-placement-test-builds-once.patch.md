**Headline:** The sector placement stress test builds its twelve large systems once instead of up to 180 times, and finishes in under three minutes.

### Changed

- `test_bughunt_sector_placement.py` built a new large B3V system for each of up to 12 placements in each of 15 seeds (up to 180 systems at 13 to 20 s each) and four large O5V systems (1 to 32 s each), so the file ran for the better part of an hour and hung the partial gate (PERF.241, TEST.417). It builds 12 large B3V systems and one large O5V system once for the module (`large_systems`) and each seed fills its own sector with the 12 in an order of its own. The law under test reads each system's Hill sphere, not what is in the system, so the same pairs of sizes meet in different places and orders; the check on every placed pair, the full-sector `ValueError` and the explicit-position bypass are unchanged.

generation-output: unchanged
