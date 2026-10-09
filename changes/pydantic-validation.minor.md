### Changed
- Every API request body and the Generate page's form fields are checked by Pydantic models with the old limits, and a bad request lists every problem field at once (an `errors` list in the API response, one joined message on the page) instead of stopping at the first. The hand-written checks are deleted (ADM.21).
