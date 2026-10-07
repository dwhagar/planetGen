// html/static/localtime.js
//
// Shows every server-rendered time in the viewer's own time zone. The
// server writes each one as <time datetime="...Z" data-local-time> with a
// UTC fallback text (planetgen/web/lib/fmt.py `utc_time_html`); this rewrites the
// text with Intl.DateTimeFormat in the browser's zone, adding the zone's
// short name, and keeps the UTC time in the title. Without script the
// page still reads correctly, labelled UTC. Other scripts that add times
// later (generatejobs.js) call window.planetgenLocalizeTimes(root).

(function () {
  "use strict";

  var formatter = null;
  try {
    formatter = new Intl.DateTimeFormat(undefined, {
      year: "numeric", month: "short", day: "numeric",
      hour: "2-digit", minute: "2-digit", timeZoneName: "short",
    });
  } catch (e) {
    formatter = null;
  }

  function localize(root) {
    if (!formatter) {
      return;
    }
    var nodes = (root || document).querySelectorAll("time[data-local-time]");
    for (var i = 0; i < nodes.length; i++) {
      var node = nodes[i];
      var when = new Date(node.getAttribute("datetime"));
      if (isNaN(when.getTime())) {
        continue;
      }
      if (!node.hasAttribute("title")) {
        node.setAttribute("title", node.textContent);
      }
      node.textContent = formatter.format(when);
    }
  }

  window.planetgenLocalizeTimes = localize;

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", function () { localize(document); });
  } else {
    localize(document);
  }
})();
