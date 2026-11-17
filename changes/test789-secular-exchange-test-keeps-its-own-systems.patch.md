**Headline:** The secular exchange generation test no longer depends on which test ran before it.

### Fixed

- TEST.789: `test_secular_exchange_generation` failed once under load with an inclination of 1.697 degrees against the stored 0.997. One test of the file moves every planet's plane by 0.7 degrees on the module-shared systems to check how the spin follows, so a test that ran after it on the same worker (the order differs under xdist) saw the moved orbits. That test now builds its own systems (`scratch_systems`); the assertions and bounds are unchanged. Test only.

generation-output: unchanged
