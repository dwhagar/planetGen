**Headline:** Four tests red on main pass again; the Gaia nearby-star query qualifies its columns.

### Fixed
- **Main reds (TEST.740, part).** The Gaia nearby-star ADQL query names its table `g` and qualifies every column, so `dec` is no longer an unquoted reserved word (same rows, same bytes). The warm-up test reads the capitalised log line the command writes. The Gaia sources test reads the VizieR fallback's current name. A planet architecture docstring cites its source by ID. The unused `run_common.write_text_output` is removed and the wiring register drops the roundness row RN1, which is read now.

generation-output: unchanged
