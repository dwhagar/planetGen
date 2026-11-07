**Headline:** The body-marker test picks its not-a-class letter from the class list instead of a fixed letter that has become a class.

### Changed

- **Body-marker test follows the class list (TEST.674).** `tests/js/bodymarkers.test.mjs` used R as the letter that is not a planet class, and R is a class now, so the test failed. It takes the first capital letter that `PLANET_CLASS_LETTERS` leaves out, so adding a class cannot break it again (it asserts a letter is left). Test only.

generation-output: unchanged
