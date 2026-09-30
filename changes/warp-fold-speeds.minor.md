### Changed

- Warp travel times follow Boss's warp curve, which matches warp^(10/3) up to about warp 9 and then climbs toward warp 10 (warp 9.995 is about 12,200 c). The NAV warp table covers warp 1, 2, 4, 8, 9, 9.5, 9.9 and 9.995.

### Added

- Dimensional fold travel times (6F^4 / (10 - F) times c) at fold 4 to 8.5 on the NAV page and in `/api/nav`'s `fold_times`.
