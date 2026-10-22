// static/generateplace.js
//
// The Add sectors form (web/templates/generate_task.html, ADM.108): shows the one choice that is made, the way to name
// a place or the grid range, and hides the rest. Without script every choice is shown and the radios still decide
// what is posted (`web/generate_page.where`).

function show(root, attribute, value) {
  for (const panel of root.querySelectorAll(`[${attribute}]`)) {
    panel.hidden = panel.getAttribute(attribute) !== value;
  }
}

function checked(root, name) {
  const radio = root.querySelector(`input[name="${name}"]:checked`);
  return radio ? radio.value : null;
}

export function init(root) {
  const update = () => {
    show(root, "data-what-panel", checked(root, "what"));
    show(root, "data-place-panel", checked(root, "place"));
    show(root, "data-grid-panel", checked(root, "mode"));
    for (const shared of root.querySelectorAll("[data-place-shared]")) {
      shared.hidden = !shared.getAttribute("data-place-shared").split(" ").includes(checked(root, "place"));
    }
  };
  root.addEventListener("change", update);
  root.classList.add("tabs-on");
  update();
}

document.querySelectorAll("[data-add-sectors]").forEach(init);
