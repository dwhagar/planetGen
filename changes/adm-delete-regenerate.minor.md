### Added
- **Delete and Regenerate buttons on sectors, systems, planets, moons,
  asteroid belts and phenomena (ADM.8).** An admin sees an Edit panel on
  the sector, system and phenomenon pages; each button asks for
  confirmation first. Regenerating rolls the object again in the same
  place and keeps its name; after a body changes, the rest of the system
  is re-checked and re-spaced until it is stable, and the page says what
  moved and anything still unstable. Deleting a sector removes everything
  in it and leaves its place in the galaxy free to generate again. New
  API routes: `DELETE`/`POST .../regenerate` on `/api/planets`,
  `/api/moons`, `/api/belts` and `/api/phenomena/<type>`, plus
  `DELETE /api/sectors/<id>/contents` and `POST /api/sectors/<id>/regenerate`.
