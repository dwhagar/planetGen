**Headline:** A running generation job now links to the systems it has made so far.

### Changed
- **New admin page "Generated so far" (ADM.150), `/admin/generated`.** It lists the sectors, systems, planets and moons a generation run has made, newest first, one tab per kind, from a start time (`?since=`, by default the start of the latest galaxy run). The job panel links to it as "Made so far this run", because the Galaxy Map is closed while a run goes. The API has `GET /api/generated` for it; it reads only the rows made since by the key, so it is quick on any galaxy and is never cached.
