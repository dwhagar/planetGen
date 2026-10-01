### Added
- **Generate a Galaxy Map block.** `generate.py galaxy --block M.I.S.SLAB [--block-layer J]` fills one drill-down block (or one of its layers), and the admin Generate page has a matching form.
- **Neighborhood radius in light-years.** The Generate page's single-sector neighborhood takes a radius in ly (13 ly and up), converted to parsecs.
- **JSON answers from the Generate page.** A post with `Accept: application/json` returns the started job's id and status URL, so the Galaxy Map can start jobs without leaving the map.
