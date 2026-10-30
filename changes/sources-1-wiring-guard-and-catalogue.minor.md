**Headline:** A test now fails when something the app builds or downloads is never read, and setup downloads the Harris cluster catalog, the NOAA solar proton events and the Baraffe 2015 isochrones.

### Added
- **The wiring guard (GEN.754).** `scripts/wiring_guard.py` finds every module, setup dataset, dataset reader, data file and stand-in registry row that nothing in the app reads, and `src/tests/test_wiring_guard.py` fails unless `docs/design/research/wiring-register.md` names it with the open item that will wire it (or says why it is read by design). ExoMol and HITRAN, the 1,453 MB of the 1,806 MB setup downloads that no generation step reads yet, are named there against GEN.850.
- **Three setup downloads and their readers (GEN.297, GEN.747, GEN.936, GEN.948).** `harris_clusters` (157 globular clusters, 49 kB), `noaa_proton_events` (265 solar proton events, 1976 to 2017, 45 kB) and `baraffe_isochrones` (30 ages of low-mass stars and brown dwarfs, 65 kB) join the catalogue as registry rows SC4 to SC6. Nothing in generation reads them yet, so generated output does not change; `./update.sh` fetches them.
- **The burn-down register (GEN.946).** `docs/design/research/burn-down-register.md` gives each of the 73 hand-tuned stand-in rows (tweak level L3 or L4) a disposition, `no data` with the reason or `replace` with the open item that will replace it, and `test_burn_down_register.py` fails for a new hand-tuned row without one. `stand-in-values.md` shows it, and the wiring status, as two new columns.
- **Stand-in rows AC1 to AC10 (GEN.724)** for the accreting compact binaries block of `tuning.py`, with a test that no constant of the block lacks a row.

### Changed
- **The millisecond pulsar registry row (GEN.724)** now says a recycled pulsar's neutron star companion is 0.2 percent, as the code has kept, not 5 percent.
