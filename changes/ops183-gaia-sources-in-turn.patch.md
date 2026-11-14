**Headline:** The Gaia 25 parsec sample is fetched at setup from four sources in turn, so the ESA archive's 2,000-row cap or a full queue no longer leaves it missing.

### Changed
- **`datasets/fetch.py` takes mirrors for a piece (OPS.183).** When a piece's own address fails or serves the wrong bytes, `fetch` tries each `Mirror` in order and logs which one served and, for the Gaia file, how many sources it holds; it is missing only when all fail, with every reason. A mirror either gives the same bytes (checked against the sha256) or converts another catalogue's answer, which `accept` judges and a `.source.json` record keeps.
- **`datasets/gaia_nearby.py` lists four:** the ESA archive's asynchronous query (with the optional `PLANETGEN_ESA_USER` and `PLANETGEN_ESA_PASSWORD` account, the password masked by SEC.48), the Heidelberg ARI mirror (byte for byte the same file), VizieR's Gaia DR3 table and the Gaia Catalogue of Nearby Stars (I/360), the last two kept only when they hold about 5,300 well-measured stars.

generation-output: unchanged
generation-affects: none
