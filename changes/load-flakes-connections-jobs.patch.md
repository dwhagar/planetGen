### Fixed

- `GET /api/databases` counts each schema's sectors and systems over the one connection it lists them with, instead of opening a connection pool per schema and keeping it; under the parallel test suite that ran MariaDB out of connections (TEST.93).
- A Generate-page job in its first seconds shows as starting, not interrupted, while its runner process hasn't yet started running the job (TEST.92).
