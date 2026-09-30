### Changed

- NAV courses read "bearing mark mark" (for example "045 mark 330") on Boss's nested reference frames: bearing 000 points at the sector's center for a course inside one sector and at the galactic core between sectors, and the mark is the elevation mod 360 (270 is straight down). The NAV page shows the course and its frame in place of the Azimuth and Altitude rows, and the NAV Map's compass arrow marks bearing 000.
- `/api/nav`'s `direct` carries `bearing_deg`, `mark_deg`, `elevation_deg` and `frame` instead of `azimuth_deg` and `altitude_deg`.
