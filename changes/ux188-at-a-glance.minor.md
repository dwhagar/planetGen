**Headline:** Planet and moon pages now open with an "At a glance" table of what a captain needs first, and leave out properties that hold nothing.

### Changed
- A planet's or moon's page has an "At a glance" table just under the render: its habitability symbols and life scores on a line of their own, then class, surface temperature, pressure, gravity, radiation, hours outside, day length, distance and year (UX.188). The full property sheet follows, as before.
- Object pages no longer list a property whose value is none or 0, except where the absence matters (no atmosphere, no magnetic field, no surface water, a quiet star, the habitability, life and rating rows). A line under the sheet says how many are not listed.
- `GET /api/objects/<ref>/sheet` rows carry `absent` and `show_absent`; the API still returns every row.
