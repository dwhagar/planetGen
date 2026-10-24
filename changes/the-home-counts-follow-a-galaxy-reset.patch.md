**Headline:** After a galaxy reset the Home page's counts drop to zero at once instead of showing the old galaxy's figures.

### Fixed
- **Home showed the old counts after a reset (UX.161, Boss's live observation).** A stored count (`countcache`, PERF.64) made before the reset was served once more while the next one was made, and the page cache (PERF.104) kept that answer; with the galaxy wiped, nothing moved the content stamp again to replace it. The stored counts are now keyed by the sector deletion epoch, which a reset moves up, so no count made before a reset is served after it, and `planetgen reset` tells the page cache at once (`dataversion.record(..., "everything")`), as a generation run does when it finishes. Test in `test_reset_counts.py`, with the caches on and off.
TODO items resolved: UX.161. Added: none.
