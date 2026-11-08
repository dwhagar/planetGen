**Headline:** On a small phone screen the Galaxy Map's info card no longer scrolls itself when it is expanded: its three values scroll in their own box and the "more" region keeps room.

### Changed

- UX.292: with a long name and long values in a window only 640 px high (390 by 640), an expanded star-system card had a first view of 328 px against a card limit of 295 px, so the "more" region got no height and the card itself scrolled. An expanded card now keeps at least 4 rem for its "more" region, and the three values of the first view shrink into a box of their own that scrolls, so the card never scrolls (the MAP.347 contract). Cards that fit are drawn as before.

generation-output: unchanged
