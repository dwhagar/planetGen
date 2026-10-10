// NAV.8: a link to /system/<id>#planet-12 (or #star-, #moon-, #belt-, #comet-) opens that body's row on the
// system page, highlights it, scrolls to it and selects the body on the System Map. Planets, moons, belts and
// comets are anchors on the system page, not pages of their own.

var TARGET_CLASS = "body-targeted";

// The map's click handlers go on once it marks itself ready (`data-ready="true"`).
function whenMapReady(callback) {
  var root = document.getElementById("sysmap-root");
  if (!root || root.dataset.ready === "true") {
    callback();
    return;
  }
  var observer = new MutationObserver(function () {
    if (root.dataset.ready === "true") {
      observer.disconnect();
      callback();
    }
  });
  observer.observe(root, { attributes: true, attributeFilter: ["data-ready"] });
}

function selectOnMap(row) {
  var match = /^(planet|moon|star|belt|comet)-(\d+)$/.exec(row.id);
  if (!match) return;
  var markers = document.querySelectorAll(".sysmap-body[data-kind='" + match[1] + "'][data-id='" + match[2] + "']");
  var marker = markers[0];
  if (!marker) return;
  // A moon sits in its planet's scene: enter that scene first by clicking the planet's marker.
  var scene = marker.closest("svg[data-scene]");
  if (scene && scene.classList.contains("sysmap-hidden")) {
    var sceneId = scene.getAttribute("data-scene");
    var opener = document.querySelector(".sysmap-body[data-scene='" + sceneId + "']");
    if (opener) opener.dispatchEvent(new MouseEvent("click", { bubbles: true }));
  }
  marker.dispatchEvent(new MouseEvent("click", { bubbles: true }));
}

function target() {
  var previous = document.querySelector("." + TARGET_CLASS);
  if (previous) previous.classList.remove(TARGET_CLASS);
  var id = decodeURIComponent(window.location.hash.slice(1));
  if (!id) return;
  var row = document.getElementById(id);
  if (!row || !row.querySelector(":scope > details.body-row")) return;
  // Open the row and every row or moon group it sits in.
  for (var node = row; node; node = node.parentElement) {
    if (node.tagName === "DETAILS") node.open = true;
  }
  var own = row.querySelector(":scope > details.body-row");
  if (own) own.open = true;
  row.classList.add(TARGET_CLASS);
  row.scrollIntoView({ block: "center" });
  whenMapReady(function () { selectOnMap(row); });
}

window.addEventListener("hashchange", target);
target();
