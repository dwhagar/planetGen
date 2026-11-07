**Headline:** The API link audit knows the belt member route and the site update route, so it is green again.

### Changed

- `src/tests/test_api_link_ids.py` (TEST.662): `api.object_belt_member_route` (GEN.1070) and `update.get_update` each get a reason in the list of endpoints the audit does not fetch (a belt member is named by the ref and rank in its path, and the update state carries no linked object).

generation-output: unchanged

TODO items resolved: TEST.662.
