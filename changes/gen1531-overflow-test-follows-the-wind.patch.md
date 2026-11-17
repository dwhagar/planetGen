**Headline:** The black hole overflow test no longer assumes a giant's wind at 5 AU is small against its overflow.

### Fixed

- GEN.1531: `test_dark_light` compared a giant's wind feed at 5 AU against a hundredth of its overflow, which holds only for the wind law of before Train 41's GEN.1625. It now checks that the far feed is under the near feed less the overflow. Test only.

generation-output: unchanged
