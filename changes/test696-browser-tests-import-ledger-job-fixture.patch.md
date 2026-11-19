**Headline:** The ten browser and menu test files that use `page_targets` import the `sample_ledger_job` fixture it now needs.

### Fixed

- `page_targets` (in `test_web_a11y.py`) gained a `sample_ledger_job` argument with OPS.271, but the ten test files that import `page_targets` by name did not import that fixture, so every test using it ended in `fixture 'sample_ledger_job' not found` (370 errors in `test_web_browser_controls.py` alone, TEST.696). Each file now imports it beside `sample_job`. No product code changed.

generation-output: unchanged
