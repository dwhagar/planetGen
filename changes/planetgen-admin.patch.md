### Changed
- **The admin modules move into `planetgen.admin` (OPS.24, step 7 of 14).** `adminAuth`, `loginThrottle`, `totp`, `qrcodegen`, `activitylog` and `adminEdits` are now `planetgen.admin.auth`, `throttle`, `totp`, `qrcode`, `activity_log` and `edits`. The common-password list and its license moved with them, and every caller moved too.
