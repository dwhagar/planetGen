**Headline:** A running generation job now links to the systems it has made so far.

### Changed
- **"Systems made this run" link on the job panel (ADM.150).** It opens the Systems list limited to systems created since the run started (`/systems?since=<Unix seconds>`; the API's `GET /api/systems` takes the same `since`). The list is built by the existing slow-page job (claimed under the memory gate, shared by everyone who asks), so an admin waits on the loading bar and the site is never held up.
