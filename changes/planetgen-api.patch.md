### Changed
- **The JSON API moves into `planetgen.api` (OPS.24, step 11 of 14).** `src/html/api/` is now `src/planetgen/api/`, with the same module names (`app`, `routes`, `auth`, `config`, ...). The Apache example drops its deny rule for the old folder, since the package sits outside the DocumentRoot. Every caller moved too.
