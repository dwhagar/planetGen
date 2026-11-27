**Headline:** The static URL tests expect the version the way the URL carries it, with the plus quoted

### Fixed

- TEST.849: `test_page_shell.py::test_static_url_appends_package_version_and_the_static_fingerprint` and `test_static_url_keeps_subdirectories` were red on main. `static_url` has quoted the version in the query since the fingerprint of the static files (`quote(token, safe='')`); the tests compared with the bare `__version__`, which held no plus until OPS.64 (PR 7315, 8.0.4653) made it `8.0.4653+db...gen...api...`. A raw plus in a query reads as a space, so the encoded `%2B` is the right form and `static_url` stays as it is; both tests now expect `quote(__version__, safe='')`. Guard files run on the change: `test_page_shell.py`, `test_fmt_bad_numbers.py`, `test_fuzz_text_rendering.py` (959 passed). Tests only.

generation-output: unchanged
