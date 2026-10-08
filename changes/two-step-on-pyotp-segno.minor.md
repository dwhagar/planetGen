### Changed

- Two-step sign-in codes are checked with `pyotp` and the setup QR code is
  drawn by `segno` (SEC.29). Secrets, the 30-second step, the one-step
  clock window and the no-reuse rule are unchanged, so enrolled admins keep
  working. The two hand-written modules (`totp.py`, `qrcode.py`, with a
  copy of Nayuki's QR generator) are deleted.
- The setup page's `otpauth://` link now reads `planetGen:<name>` where it
  read `planetGen%3A<name>` (the same label, written the usual way).
