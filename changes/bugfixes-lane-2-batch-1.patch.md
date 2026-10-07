### Fixed

- Admin scripts reject a `--mysql-port` outside 1 to 65535 with a usage error before connecting, instead of a connection error later (OPS.6).
- Whole numbers stay in plain digits up to 999,999 and go scientific from 7 digits; numbers shown with decimals still go scientific from 5 whole digits (UX.36).
- The System Map side panel shows a planet's or moon's surface pressure under its surface temperature (MAP.117).
