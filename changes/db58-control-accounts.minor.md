**Headline:** The control database gets the account columns, the single-Owner rule and the token tables that user accounts will use.

### Changed
- **Control schema v15 (DB.58): roles, email and status on accounts, `account_tokens` and `invite_uses`.** Each account has a `role` (user, moderator, admin or owner), an optional unique `email` and a `status` (active or disabled). Exactly one account can be the owner, and MySQL refuses a second one; the account the install makes is the owner, and in an existing control database the lowest account id becomes the owner and the others stay admins. `account_tokens` (hashed one-time tokens with a purpose and an expiry) and `invite_uses` (who registered through an invite link) are made empty. Nothing reads or writes the new columns or tables yet.

TODO items resolved / Added: DB.58 / Added: none
