**Headline:** The sector map's gravity layer now says "not computed yet" when the galaxy has no shape, instead of failing with a "data unavailable" error, and a page whose data call is refused with a conflict now shows what the refusal said.

### Changed
- **MAP.405.** `GET /sector/<id>/gravity` answers 409 with a JSON message ("The gravity layer is not computed yet: the galaxy has no shape yet.") when the API says the galaxy has no shape, in place of a 502.
- An API conflict (409) that no page catches itself now shows a "Not possible now" page with the message the API gave (it was a 502 "Data unavailable" page).

TODO items advanced: MAP.405.
