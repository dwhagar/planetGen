**Headline:** A search that takes a while to prepare, such as all M class planets, now finishes and shows its results instead of preparing over and over.

### Fixed
- **A slow search no longer loops on "preparing" (UX.234).** A page that is slow to build is built by a queued job and shown when the job is done. When the page cache would not keep the finished answer (data changed while it was built, as it does during a fill, or Redis had no room for the cache), the page found nothing, queued the same job again and showed "preparing" again, for ever. A finished job now keeps its answer for the page that waited for it (five minutes), and the page shows it.
- **A search waits for all its answers at once (UX.234).** The Search page needs eight answers (the facets and one per result list). It asked for them one after another, a new "preparing" screen for each. It now starts all eight together and waits for the slowest.
TODO items resolved: UX.234 (TODO closes it after reading this). Added: none.
