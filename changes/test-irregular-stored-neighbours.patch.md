**Headline:** A test of saving and reading back a planet's irregular moons no longer assumes its neighbouring planets hold none.

### Changed

- The stored-irregular-moons test compares each neighbouring planet with the irregular set the generator gave it before saving (none, or an empty set), instead of requiring none. Since GEN.1406 the sample system's outer planets hold an empty set and the test failed on main. Test only.

generation-output: unchanged
