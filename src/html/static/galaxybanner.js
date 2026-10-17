// html/static/galaxybanner.js
//
// Keeps the galaxy banner (UX.3) current. base.html draws the banner when the page is made, from what the
// server knew then; this asks /galaxy-status every 30 seconds, and again when a hidden tab comes back, so the
// banner appears when a background run starts and goes when it ends, without a reload. The server keeps one
// answer for a few seconds, so many open tabs cost it little. Without script the banner still shows as drawn.

(function () {
  "use strict";

  var banner = document.querySelector(".galaxy-banner[data-status-url]");
  if (!banner || typeof fetch !== "function") {
    return;
  }
  var url = banner.getAttribute("data-status-url");
  var seconds = parseInt(banner.getAttribute("data-poll-seconds"), 10) || 30;
  var timer = null;

  function show(status) {
    if (status && status.active && status.text) {
      banner.textContent = status.text;
      banner.hidden = false;
    } else {
      banner.hidden = true;
      banner.textContent = "";
    }
  }

  function poll() {
    fetch(url, { credentials: "same-origin", cache: "no-store", headers: { Accept: "application/json" } })
      .then(function (response) {
        return response.ok ? response.json() : null;
      })
      .then(function (status) {
        if (status) {
          show(status);
        }
      })
      .catch(function () {
        // leave the banner as it is; the next poll tries again
      });
  }

  function schedule() {
    window.clearInterval(timer);
    timer = null;
    if (!document.hidden) {
      timer = window.setInterval(poll, seconds * 1000);
    }
  }

  document.addEventListener("visibilitychange", function () {
    if (!document.hidden) {
      poll();
    }
    schedule();
  });
  schedule();
})();
