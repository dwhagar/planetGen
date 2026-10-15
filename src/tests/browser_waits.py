# tests/browser_waits.py

"""
A longer wait for every browser test's page (TEST.126, TEST.138, TEST.144).

Playwright gives each wait (`goto`, `click`, `wait_for_selector`, `wait_for_function`, a screenshot) 30 s. The
Galaxy and Sector Map tests draw WebGL in software (SwiftShader) on the same cores the other workers use, so on a
busy machine one honest wait can pass 30 s and the test fails with a `TimeoutError` that passes alone. Every
context a test makes through `Browser.new_context` gets `DEFAULT_TIMEOUT_MS` as its default instead. It is a ceiling,
not a delay: a wait still returns the moment its condition holds, an assertion is unchanged, and a condition that
never holds still fails (after 2 minutes, well inside the 900 s a test is given). A wait that passes its own
`timeout=` keeps it.
"""

DEFAULT_TIMEOUT_MS = 120_000


def install(browser_class=None):
    """Makes the contexts of `browser_class` (Playwright's sync `Browser` by default) use the longer wait.

    Idempotent. Does nothing when Playwright is not installed.

    Returns:
        bool: True when the class is patched (now or before).
    """
    if browser_class is None:
        try:
            from playwright.sync_api import Browser as browser_class
        except ImportError:
            return False
    if getattr(browser_class.new_context, "_longer_waits", False):
        return True
    original = browser_class.new_context

    def new_context(self, *args, **kwargs):
        context = original(self, *args, **kwargs)
        context.set_default_timeout(DEFAULT_TIMEOUT_MS)
        return context

    new_context._longer_waits = True
    browser_class.new_context = new_context
    return True
