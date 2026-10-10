// static/generatestarmix.js
//
// The Generate page's star mix (TODO ADM.45): single stars, close binaries and
// wide pairs are shares of systems that must total exactly 100%. This keeps a
// running total next to them and says in words which way to move; the server
// refuses a set that is off (planetgen.generation.prevalence.star_mix_prevalences).

export const TOLERANCE = 0.01;

// The sentence for shares totalling `total`.
export function totalText(total) {
  if (!Number.isFinite(total)) return "Every share needs a number.";
  const rounded = Math.round(total * 1e4) / 1e4;
  if (Math.abs(total - 100) <= TOLERANCE) return "Total " + rounded + "%. That is right.";
  const verb = total > 100 ? "Lower" : "Raise";
  return "Total " + rounded + "%, not 100%. " + verb + " one of them by " + Math.round(Math.abs(total - 100) * 1e4) / 1e4 + " points.";
}

export function init() {
  return Array.from(document.querySelectorAll("[data-star-mix]")).map((box) => {
    const inputs = Array.from(box.querySelectorAll("input"));
    const out = box.querySelector("[data-star-mix-total]");
    const update = () => {
      const total = inputs.reduce((sum, input) => sum + (input.value === "" ? NaN : Number(input.value)), 0);
      out.textContent = totalText(total);
      out.classList.toggle("error", !(Math.abs(total - 100) <= TOLERANCE));
    };
    inputs.forEach((input) => input.addEventListener("input", update));
    update();
    return update;
  });
}

init();
