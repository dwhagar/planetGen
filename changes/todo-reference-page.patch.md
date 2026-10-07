### Added
- **A searchable TODO reference page**: `docs/plan/todo-reference.html`
  explains every TODO ID ever issued (full text, status, phase, thread,
  prerequisites and what each unblocks) with search and filters.
  `python scripts/build_todo_docs.py` rebuilds it and the tier plan from
  `docs/TODO.md` and the plan files, and a test fails when either page is
  out of date.
