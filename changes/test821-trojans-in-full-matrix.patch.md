**Headline:** The full-matrix test no longer mistakes a trojan on its giant's orbit for an overlap.

### Fixed
- **`test_full_matrix` checks a trojan (GEN.1150, GEN.1151) against its giant and leaves it out of the spacing of the bodies (TEST.821).** The fuzz order test has done so since TEST.818; the check is now one helper (`tests/trojan_support.assert_trojans_on_their_giants`) both call. Every other spacing check stays as strict as it was. No generator change.
