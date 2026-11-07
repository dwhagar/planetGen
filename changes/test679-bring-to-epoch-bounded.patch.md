**Headline:** The bring-to-epoch end-to-end test builds one system a cell and makes three cells again, so it takes seconds instead of more than seventeen minutes.

### Changed

- **`test_a_sector_made_under_an_older_epoch_is_really_generated_again` bounded (TEST.679).** It was not stuck: a system takes about fifteen seconds to build on the tree of train 15, a profile-sized sector holds a dozen, the first sector throws 78 rogue planets into sectors of their own (GEN.1111), and the test marked all 79 sectors behind and made every one again, each filling with stars whose rogues landed in more cells. The hang guard killed it after 17 minutes alone. Every cell the test makes (the first and each making again) now holds one system, and the sector and two of its landing cells are behind; the rest are at the running epoch. It still deletes the sector, makes it again from its address, checks it holds the same systems and carries the running epoch, and that a second run finds nothing to do. The file runs in two minutes. Test only.

generation-output: unchanged
