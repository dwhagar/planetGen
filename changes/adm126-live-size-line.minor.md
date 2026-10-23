**Headline:** The Generate form shows how big the run will be while you fill it in.

### Added
- **A live size line under the Add sectors and Check forms (ADM.126).** A moment after a box changes, the line says how many sectors, star systems and stars the run would make, how much it adds to the database and how long it takes, or why the values cannot be sized yet. It is the estimate the Review step shows (`planetgen galaxy --estimate-only`, cached by the form's values); nothing is started or written. A command typed in the admin palette shows the same figures on its row.
- **`GET /admin/generate/task/<key>/size`**, admins only, returns that estimate as JSON for the form's fields.
