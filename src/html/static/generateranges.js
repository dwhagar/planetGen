// static/generateranges.js
//
// The admin Generate page's sliders (web/templates/generate.html, GEN.183): shows a range input's value in
// the <output> its data-range-output names, as the slider moves. Without JavaScript the output keeps the
// value the server drew.

for (const input of document.querySelectorAll("input[type=range][data-range-output]")) {
  const output = document.getElementById(input.dataset.rangeOutput);
  if (!output) continue;
  const show = () => { output.textContent = input.value; };
  input.addEventListener("input", show);
  show();
}
