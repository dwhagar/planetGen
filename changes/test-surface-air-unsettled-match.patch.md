**Headline:** The surface-air storage test now skips an unsettled air by its gases as well as its pressure.

### Changed

tests/test_surface_air.py: the test that reads the stored density, scale height and greenhouse multiplier of every generated body skipped any body whose pressure equalled an unsettled air's pressure to a part in a billion. A Train 45 sample (seeds 0 to 9 after GEN.1520) held one air that did not settle at exactly 10 kPa and three bodies of 0.003 Earth masses at 8 K with exactly 10 kPa of another mixture, so four bodies were skipped against a bound of two. The skip now also requires the stored gases that hold 1 percent of the pressure or more to match the unsettled air's to a part in a thousand. The bound (max(1, bodies // 20)) is unchanged. Test only; no generation output changes.

generation-output: unchanged
