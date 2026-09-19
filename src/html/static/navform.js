// html/static/navform.js
//
// Loaded on every page (html/lib/page.py's `render`) so an element
// carrying `data-nav-target`/`data-nav-params` navigates by POSTing a
// throwaway hidden form instead of a plain `<a href>` -- the one place in
// html/ that still needs this after `lib/page.py`'s `post_link` moved
// every *plain* HTML link to a real, no-JS-required POST `<form>` (a
// `<form>` can't be nested inside an SVG shape the way galaxy.py's/
// sector.py's/nav.py's map markers are, so those stay real `<a>` tags,
// intercepted here instead). `data-nav-params` is a JSON object string
// (build it with `json.dumps` server-side, HTML-attribute-escaped the
// same way any other data-* value is).
//
// This is the one part of html/'s navigation that now requires
// JavaScript -- every plotted sector/system/phenomenon also has a plain,
// no-JS-required row in the table below its own map (see sector.py/
// galaxy.py/nav.py), so a map marker is never the *only* way to reach
// something.

(function () {
  "use strict";

  function submitNav(action, params) {
    var form = document.createElement("form");
    form.method = "post";
    form.action = action;
    form.style.display = "none";
    Object.keys(params || {}).forEach(function (key) {
      var input = document.createElement("input");
      input.type = "hidden";
      input.name = key;
      input.value = params[key];
      form.appendChild(input);
    });
    document.body.appendChild(form);
    form.submit();
  }

  document.addEventListener("click", function (event) {
    var el = event.target.closest("[data-nav-target]");
    if (!el) {
      return;
    }
    event.preventDefault();
    var params = {};
    try {
      params = JSON.parse(el.dataset.navParams || "{}");
    } catch (err) {
      params = {};
    }
    submitNav(el.dataset.navTarget, params);
  });

  window.planetgenSubmitNav = submitNav;
})();
