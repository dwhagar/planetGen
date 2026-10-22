**Headline:** The database gets the tables for comments and notices on edit requests, ready for the 10.0 roles.

### Changed
- **`request_comments` and `request_notifications` (DB.39, schema v126).** The conversation on an edit request (who wrote what, in which role) and the notices owed to accounts about it, each with its plan. A request is a plan row (DB.36), so the requester, approver and decision are not repeated. The migration only makes the tables; nothing reads them yet.

TODO items resolved / Added: DB.39 / Added: none
