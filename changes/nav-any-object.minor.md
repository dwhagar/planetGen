### Added
- NAV takes any object: `/nav` and `GET /api/nav` accept a star, planet, moon, asteroid belt or comet as either end, written as an object reference (`planet:40`). A course to or from a body gets legs: out of its system to the heliopause, between the systems, into the destination body (all in the System Local Frame inside a system). Two objects in one system are one in-system leg, using the same warp and fold tables with a note (NAV.16).

### Changed
- `GET /api/nav` takes `from`/`to` as object references; the `from_kind`/`to_kind`/`from_type`/`to_type` parameters are gone (`from=nebula:2` replaces `from=2&from_kind=phenomenon&from_type=nebula`).
