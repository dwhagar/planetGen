### Fixed

- System page: an asteroid belt's row in the object list shows only its
  density, its range ("2.1 AU to 3.3 AU") and its three largest minerals,
  no longer the belt's distance twice (UX.19). No row (planets, moons,
  belts) shows the zone any more; the zone stays in each body's own
  description.
- API: `GET /api/systems/<id>` gives each belt its `composition`
  (`{component, concentration}`, largest share first).
