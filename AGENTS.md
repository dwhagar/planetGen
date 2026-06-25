# Repository Agent Rules (`AGENTS.md`)

## 1. Execution Architecture and Loop Mitigation
- **Deterministic Failure:** If a file write, regex patch, or refactoring tool fails to match or execute 2 consecutive times, terminate the current execution loop immediately and output the active state to the user. Do not attempt semantic variations of the same edit.
- **State Verification:** Before executing any file modification, parse the local file buffer to ensure the line anchors and structural offsets match the intended target exactly.
- **Decomposition Threshold:** Tasks requiring modifications across more than 2 distinct files or exceeding 50 lines of code change must be broken into discrete, sequential phases. Execute and verify Phase $N$ completely before initiating Phase $N+1$.

## 2. Python Code Generation Guardrails
- **Surgical Block Edits:** Modify files via narrow, targeted block replacements or specific line edits. Do not rewrite, regenerate, or output unaffected adjacent functions or classes.
- **AST Integrity:** Ensure all Python blocks maintain absolute syntactic correctness. All indentation levels (4-space standard), block colon terminators, and closing parentheses, brackets, or braces must be structurally balanced before committing a write tool action.
- **No Structural Truncation:** Do not use ellipsis markers, placeholders, or code summaries (e.g., `# ... rest of function unchanged`) within code generation blocks. Output the complete implementation block targeted for modification.

## 3. Strict Ordered Phase Protocol
- **Phase 1: Logic Implementation:** Execute all algorithmic, control flow, and structural modifications first. Do not touch, add, or alter any comments, docstrings, or type hints during this phase.
- **Phase 2: Validation:** Verify code correctness against the active file environment.
- **Phase 3: Documentation and Robust Commenting:** Update or generate necessary docstrings (following Google or Sphinx style conventions based on file context) and inline comments. 

## 4. Documentation Stability and Comment Quality
- **Comment Idempotency:** Once a docstring or inline comment is established and sufficient, it is a frozen state. Do not modify, reword, or reformat it in subsequent iterations unless the underlying functional logic it describes has explicitly changed.
- **Technical Rigor:** Ensure inline comments explain the underlying rationale, mathematical properties, and state mutations rather than restating the Python syntax. Avoid high-level prose or vague summaries.
- **Format Consistency:** Adhere strictly to the pre-existing comment style and syntax patterns found within the target file. Do not toggle between different phrasing styles or equivalent comment structures.

## 5. Architectural Planning Protocol (Design-First)
- **Pre-Execution Gate:** When a prompt asks for structural options, architectural variations, or an execution plan, you are strictly prohibited from invoking file-write or patch tools.
- **Option Generation Requirements:** Present exactly 2 or 3 distinct technical approaches to solve the problem. For each approach, provide:
  1. A concise, low-level technical summary of the implementation strategy.
  2. The algorithmic, memory, or complexity trade-offs involved ($O(N)$ implications, dependency impacts, etc.).
  3. A brief, pseudo-code structural block illustrating the pattern.
- **Explicit Standby State:** Conclude the response by explicitly asking which option to execute. Enter a non-interactive standby state until the user provides the selection.

## 6. Monolithic Module Dependency Guardrails
- **Namespace Anchoring:** When utilizing whole-module imports (e.g., `import config`), do not perform autonomous repository searches to discover internal sub-properties.
- **Explicit Property Restraint:** Restrict property utilization strictly to the specific attributes explicitly defined within the active user prompt's `[DEPENDENCY MAP]`. Treat these explicit attributes as static string tokens to eliminate line-matching ambiguity and prevent AST verification loops in the host IDE process space.