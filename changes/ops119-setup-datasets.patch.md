**Headline:** The test machine setup now downloads the datasets it needs and says so when it cannot.

### Fixed
- **Setup fetched nothing and still ended well (OPS.119).** The cloud environment setup ran the dataset download with `python -I`, which also drops `PYTHONPATH`, so `planetgen` could not be imported, the step only warned, and setup finished as if it had worked. It now runs with `-P` (no current directory, `PYTHONPATH` kept), checks the result with `datasets status`, and ends with a SETUP INCOMPLETE line when the data is not all there. The test preflight also fetches every dataset the math check reads (isochrones, spectra, the helium-hydrogen equation of state, HYG and Mamajek, besides the white-dwarf tables), where it fetched only the white-dwarf data and the math check then stopped a fresh machine's run with 16 DataMissing failures.

generation-output: unchanged

TODO items resolved: OPS.119.
