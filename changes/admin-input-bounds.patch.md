### Security
- **Admin generation inputs have upper bounds.** The Generate page, the
  one-off system page, the API and `generate.py` now reject a radius over
  200 pc (about 652 ly for the API's generate-neighborhood `radius_ly`,
  down from 1,000,000 ly), a ring or highest ring over 100,000, a ring
  `--limit` over 628,322 (ring 100,000's slot count), and more than 500
  orbital slots (also in a `--system-file`). All four share the constants
  in `src/stellarObjects/generationLimits.py`, and the page's number
  inputs carry them as `max`.
