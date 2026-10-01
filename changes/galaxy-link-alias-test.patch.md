### Fixed
- **Test suite green again.** The route fuzz test that tries ids like `0001` or `١` on every page now expects the "Show on Galaxy Map" links (`/sector/<id>/galaxy`, `/system/<id>/galaxy`) to answer with their redirect to the Galaxy Map instead of failing on it. The pages themselves were already correct.
