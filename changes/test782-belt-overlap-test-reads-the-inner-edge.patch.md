**Headline:** The belt overlap test now checks the belt's inner edge, the way the spacing pass measures it.

### Fixed
- TEST.782: `test_validate_system_belt_overlap_correction_lands_exactly_past_the_belt` expected the pushed belt's centre to land at the outer edge of the belt inside it plus the separation. Since GEN.1592 a belt is measured and pushed by its inner edge (`shift_belt` moves the centre and both edges together), so the inner edge lands exactly there and the centre sits half a width further out (2.45 AU for a 0.4 AU belt past 2.25). The test now asserts the inner edge, the kept width and the centre; the generator is unchanged.

generation-output: none
