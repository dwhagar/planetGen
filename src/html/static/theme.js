// html/static/theme.js
//
// Light / dark / system theme switch. Loaded as a plain blocking script in
// the <head>, before style.css (html/lib/page.py's `head_html`), so the
// saved choice is applied to <html data-theme="..."> before the first
// paint and the page never flashes the OS theme first. A separate file
// because the site's Content-Security-Policy allows no inline script.
//
// "system" means no data-theme attribute at all: style.css then follows
// prefers-color-scheme. localStorage can be missing or throw (private
// windows, blocked site data), so every access is guarded; without it the
// switch still works for the current page.
//
// The map canvases read their colours from the CSS tokens when they start,
// so they pick up a change on the next page view.

(function () {
  "use strict";

  var KEY = "planetgen-theme";
  var ORDER = ["system", "light", "dark"];
  var LABELS = { system: "System", light: "Light", dark: "Dark" };
  var root = document.documentElement;

  function saved() {
    try {
      var value = window.localStorage.getItem(KEY);
      return ORDER.indexOf(value) >= 0 ? value : "system";
    } catch (e) {
      return "system";
    }
  }

  function store(value) {
    try {
      if (value === "system") {
        window.localStorage.removeItem(KEY);
      } else {
        window.localStorage.setItem(KEY, value);
      }
    } catch (e) {
      // Not remembered; still applied to this page.
    }
  }

  // Browser UI colour: the page ships one theme-color per OS scheme
  // (media="(prefers-color-scheme: ...)"). An explicit choice points both
  // at that theme's colour; "system" restores each one's own.
  function applyThemeColor(theme) {
    var metas = document.querySelectorAll('meta[name="theme-color"]');
    var chosen = null;
    for (var i = 0; i < metas.length; i++) {
      var meta = metas[i];
      if (!meta.hasAttribute("data-own-content")) {
        meta.setAttribute("data-own-content", meta.getAttribute("content"));
      }
      if (theme !== "system" && (meta.getAttribute("media") || "").indexOf(theme) >= 0) {
        chosen = meta.getAttribute("data-own-content");
      }
    }
    for (var j = 0; j < metas.length; j++) {
      metas[j].setAttribute("content", chosen || metas[j].getAttribute("data-own-content"));
    }
  }

  function apply(theme) {
    if (theme === "system") {
      root.removeAttribute("data-theme");
    } else {
      root.setAttribute("data-theme", theme);
    }
  }

  var current = saved();
  apply(current);

  function wire() {
    applyThemeColor(current);
    var buttons = document.querySelectorAll("[data-theme-toggle]");
    function label() {
      for (var i = 0; i < buttons.length; i++) {
        buttons[i].textContent = "Theme: " + LABELS[current];
        buttons[i].setAttribute("title", "Colour theme: " + LABELS[current] + ". Click to change.");
      }
    }
    for (var i = 0; i < buttons.length; i++) {
      buttons[i].hidden = false;
      buttons[i].addEventListener("click", function () {
        current = ORDER[(ORDER.indexOf(current) + 1) % ORDER.length];
        apply(current);
        applyThemeColor(current);
        store(current);
        label();
      });
    }
    label();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", wire);
  } else {
    wire();
  }
})();
