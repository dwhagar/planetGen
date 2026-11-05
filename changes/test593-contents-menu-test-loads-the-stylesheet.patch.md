**Headline:** The contents menu test now loads the site stylesheet, so it checks the one-row menu the pages really show and passes.

### Fixed
- `test_the_contents_menu_is_one_row_when_it_fits_and_the_grouped_list_when_it_does_not` (TEST.593) built its page with no stylesheet. The one-row layout is the site's own CSS (`.app-main.toc-strip`), so without it every contents list wrapped inside the nav and counted as fitting, and the 30-section page got the row. With `style.css` loaded the page does what the test says: 3 links fit one row (974 of 974 px), 60 links get the grouped menu. No page change was needed; the test was the broken part.
