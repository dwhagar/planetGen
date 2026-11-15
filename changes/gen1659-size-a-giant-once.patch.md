**Headline:** A giant is sized once, on its final orbit, where it was sized twice.

### Changed
- `_settle_body` no longer sizes a giant before the tides of its orbit and again after them; the second call, on the final orbit, already gave the final radius (GEN.1659). A close-in giant whose first call moved the radius could feed the tides a different radius, so the generated output may change for a few giants.

generation-output: changed
generation-affects: fill
