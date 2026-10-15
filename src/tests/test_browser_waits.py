# tests/test_browser_waits.py

"""The longer default wait for browser contexts (tests/browser_waits.py, TEST.126)."""

from tests import browser_waits


class _Context:
    def __init__(self):
        self.timeout = None

    def set_default_timeout(self, ms):
        self.timeout = ms


class _Browser:
    def new_context(self, *args, **kwargs):
        context = _Context()
        context.args, context.kwargs = args, kwargs
        return context


def test_a_new_context_gets_the_longer_default_wait_and_keeps_its_arguments():
    assert browser_waits.install(_Browser)
    context = _Browser().new_context(viewport={"width": 1}, color_scheme="dark")
    assert context.timeout == browser_waits.DEFAULT_TIMEOUT_MS
    assert context.kwargs == {"viewport": {"width": 1}, "color_scheme": "dark"}


def test_installing_twice_wraps_once():
    class Other(_Browser):
        pass

    browser_waits.install(Other)
    first = Other.new_context
    assert browser_waits.install(Other) and Other.new_context is first
