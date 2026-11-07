**Headline:** The locked nltk is 3.10.4, and the API refuses "nan" and "inf" in its numbers with a 400 instead of a server error.

### Changed

- `requirements.lock` and `requirements-server.lock` (SEC.40, audit finding S-L1): nltk 3.10.3 to 3.10.4 (advisory PYSEC-2026-3740, not reachable: nltk only reads the words corpus in `names/wordlists.py`), with the hashes of the 3.10.4 wheel and source archive; its dependencies are unchanged. `pip-audit` finds no known vulnerability in either lock, and the name tests pass on 3.10.4.
- `util/checks.py` `finite_float` (SEC.41, finding S-L2): the one place a request's text becomes a float, refusing "nan", "inf" and "-inf" with a `ValueError`. `GET /api/galaxy/made` (`since`, `until`), `GET /api/galaxy/cloud-volume` and the web `/galaxy/cloud-volume` use it, so they answer the usual 400 where "nan" used to reach MySQL as a parameter and answer an anonymous 500 (and the web route answered 502). The other routes that read floats already refused non-finite numbers.

generation-output: unchanged
