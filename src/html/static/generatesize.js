// static/generatesize.js
//
// The live "how big" line under a Generate form (ADM.126): the fields of the form go to the server a moment after the
// last change, and the line shows how many sectors, star systems and stars the run would make, what it adds to the
// database and how long it takes, or why the values cannot be sized yet. The server answers with the estimate the
// Review step shows (`planetgen galaxy --estimate-only`, cached by the form's values); nothing is started or written.

const WAIT_MS = 400;
const NOT_SIZED = "Not sized yet: ";

function show(note, text, isError) {
  note.textContent = text;
  note.classList.toggle("error", Boolean(isError));
  note.hidden = !text;
}

function draw(note, answer) {
  if (answer.state === "ready") {
    show(note, answer.refused && answer.refusal ? `${answer.line} ${answer.refusal}` : answer.line, answer.refused);
  } else if (answer.state === "invalid") {
    show(note, NOT_SIZED + answer.message, false);
  } else {
    show(note, "", false);
  }
}

for (const note of document.querySelectorAll("[data-size-report]")) {
  const form = note.closest("form");
  let timer = null;
  let sent = 0;
  let pending = null;
  const ask = async () => {
    const query = new URLSearchParams();
    for (const [name, value] of new FormData(form)) {
      if (name !== "csrf_token" && typeof value === "string") query.append(name, value);
    }
    const mine = ++sent;
    if (pending) pending.abort();
    pending = new AbortController();
    note.setAttribute("aria-busy", "true");
    try {
      const response = await fetch(`${note.dataset.url}?${query}`, {
        credentials: "same-origin", headers: { Accept: "application/json" }, signal: pending.signal,
      });
      const answer = response.ok ? await response.json() : { state: "none" };
      if (mine === sent) draw(note, answer);
    } catch (error) {
      if (mine === sent && error.name !== "AbortError") show(note, "", false);
    } finally {
      if (mine === sent) note.removeAttribute("aria-busy");
    }
  };
  const later = () => { clearTimeout(timer); timer = setTimeout(ask, WAIT_MS); };
  form.addEventListener("input", later);
  form.addEventListener("change", later);
  note.hidden = true;
  ask();
}
