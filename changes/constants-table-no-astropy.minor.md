### Changed
- The physical constants are a table in `physics/constants.py` (CODATA 2018 and IAU 2015 values) instead of being read from astropy, so astropy is no longer installed on a server: each process starts about 280 ms sooner and uses about 29 MB less, and the install is 53 MB smaller. No stored value changes. astropy stays in the `test` extra, where `test_astropy_constants.py` holds the table against astropy's own values (GEN.211).

### Removed
- The astropy and pyerfa run-time dependencies (`setup.py`, `install-python-deps.sh`, both lock files).
