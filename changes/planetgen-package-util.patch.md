### Changed
- **The code starts moving into the `planetgen` package (OPS.24, step 1 of 14).** `log`, `appconfig` and `serialization` are now `planetgen.util.log`, `planetgen.util.appconfig` and `planetgen.util.serialization`. The version file is now `src/planetgen/_version.py`. Every caller moved with them, with no compatibility stubs left behind.
