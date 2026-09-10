// html/static/systemmap.js
//
// Click-for-info and drill-into-moons behavior for the "System Map" panel
// built by `lib/systemmap.py`. Unlike `sectormap.js` (which has to resolve
// clicks by geometry through a rotated 3D `preserve-3d` stack -- see that
// file's own comment for why), this map is flat, static, fixed-size SVG
// with no rotation/scroll/zoom, so a plain event-target lookup is all
// clicking needs.
//
// Every planet/moon/belt/star is one `<g data-kind="..." data-*="...">` --
// clicking (or Enter/Space on a focused one) either fills the info side
// panel from its `data-*` attributes, or -- for a planet with moons,
// marked with `data-scene="planet-<id>"` -- swaps which `<svg data-scene>`
// is visible so that planet takes the star's place with its own moons
// arranged around it. Built with plain DOM calls (never innerHTML with
// unescaped content), same as `sectormap.js`, since every data-* value is
// still database content.

(function () {
  "use strict";

  function addField(dl, label, value) {
    if (!value && value !== 0) {
      return;
    }
    var dt = document.createElement("dt");
    dt.textContent = label;
    var dd = document.createElement("dd");
    dd.textContent = value;
    dl.appendChild(dt);
    dl.appendChild(dd);
  }

  function classField(el) {
    var cls = el.dataset.class;
    if (!cls) {
      return "";
    }
    return el.dataset.classdesc ? "Class " + cls + " -- " + el.dataset.classdesc : "Class " + cls;
  }

  function showInfo(el) {
    var panel = document.getElementById("sysmap-info");
    if (!panel) {
      return;
    }
    panel.textContent = "";

    var heading = document.createElement("h3");
    heading.textContent = el.dataset.name || "Unknown";
    panel.appendChild(heading);

    var dl = document.createElement("dl");
    var kind = el.dataset.kind;
    if (kind === "star") {
      addField(dl, "Role", el.dataset.role);
      addField(dl, "Star type", el.dataset.type);
      addField(dl, "Temperature", el.dataset.temp);
      addField(dl, "Mass", el.dataset.mass);
      addField(dl, "Radius", el.dataset.radius);
      addField(dl, "Luminosity", el.dataset.lum);
    } else if (kind === "belt") {
      addField(dl, "Density", el.dataset.density);
      addField(dl, "Distance", el.dataset.distance);
      addField(dl, "Composition", el.dataset.composition);
    } else {
      addField(dl, "Class", classField(el));
      addField(dl, "Type", el.dataset.bodytype);
      addField(dl, "Zone", el.dataset.zone);
      addField(dl, "Distance", el.dataset.distance);
      addField(dl, "Period", el.dataset.period);
      addField(dl, "Gravity", el.dataset.gravity);
      addField(dl, "Life Chemistry", el.dataset.life);
      if (kind === "moon") {
        addField(dl, "Orbits", el.dataset.parent);
      } else if (el.dataset.moons) {
        addField(dl, "Moons", el.dataset.moons);
      }
    }
    panel.appendChild(dl);

    if (el.dataset.scene) {
      var hint = document.createElement("p");
      hint.className = "hint";
      hint.textContent = "Click again to view its moon system.";
      panel.appendChild(hint);
    }
  }

  function resetInfo(panel) {
    panel.textContent = "";
    var hint = document.createElement("p");
    hint.className = "hint";
    hint.textContent = "Click a star, planet, moon, or asteroid belt for details.";
    panel.appendChild(hint);
  }

  function initSystemMap(root) {
    var info = document.getElementById("sysmap-info");
    var crumb = document.getElementById("sysmap-crumb");
    if (!info) {
      return;
    }
    var scenes = Array.prototype.slice.call(root.querySelectorAll(".sysmap-svg"));

    function sceneFocusName(sceneEl) {
      var self = sceneEl.querySelector('[data-self="true"]');
      return self ? self.dataset.name : "";
    }

    function showScene(sceneId) {
      var active = null;
      scenes.forEach(function (scene) {
        var isActive = scene.dataset.scene === sceneId;
        scene.classList.toggle("sysmap-hidden", !isActive);
        if (isActive) {
          active = scene;
        }
      });
      if (!active) {
        return;
      }

      crumb.textContent = "";
      if (sceneId !== "system") {
        var backBtn = document.createElement("button");
        backBtn.type = "button";
        backBtn.className = "starmap-btn sysmap-back-btn";
        backBtn.textContent = "← Back to system";
        backBtn.addEventListener("click", function () {
          showScene("system");
        });
        crumb.appendChild(backBtn);
        var name = sceneFocusName(active);
        if (name) {
          var current = document.createElement("span");
          current.className = "sysmap-crumb-current hint";
          current.textContent = "Moons of " + name;
          crumb.appendChild(current);
        }
      }

      var self = active.querySelector('[data-self="true"]');
      if (self) {
        showInfo(self);
      } else {
        resetInfo(info);
      }
    }

    root.addEventListener("click", function (event) {
      var el = event.target.closest("[data-kind]");
      if (!el) {
        return;
      }
      if (el.dataset.scene) {
        showScene(el.dataset.scene);
      } else {
        showInfo(el);
      }
    });

    root.addEventListener("keydown", function (event) {
      if (event.key !== "Enter" && event.key !== " ") {
        return;
      }
      var el = event.target.closest("[data-kind]");
      if (!el) {
        return;
      }
      event.preventDefault();
      if (el.dataset.scene) {
        showScene(el.dataset.scene);
      } else {
        showInfo(el);
      }
    });

    showScene("system");
  }

  var rootEl = document.getElementById("sysmap-root");
  if (rootEl) {
    initSystemMap(rootEl);
  }
})();
