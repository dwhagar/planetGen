// static/generatesize.js
//
// The live "how big" line under a Generate form (ADM.126): the fields of the form go to the server a moment after the
// last change, and the line shows how many sectors, star systems and stars the run would make, what it adds to the
// database and how long it takes, or why the values cannot be sized yet. The server answers with the estimate the
// Review step shows (`planetgen galaxy --estimate-only`, cached by the form's values); nothing is started or written.
// A big galaxy takes minutes: the server then answers "working" at once and the line asks again every few seconds.

const WAIT_MS = 400;
const NOT_SIZED = "Not sized yet: ";
const WORKING = "Working out how big this would be. ";
const RETRY_MS = 3000;
const MAX_RETRIES = 400;   // the server works the estimate out in the work queue (ADM.141); about twenty minutes

function show(note, text, isError) {
  note.textContent = text;
  note.classList.toggle("error", Boolean(isError));
  note.hidden = !text;
}

function draw(note, answer) {
  if (answer.state === "ready") {
    show(note, answer.refused && answer.refusal ? `${answer.line} ${answer.refusal}` : answer.line, answer.refused);
  } else if (answer.state === "working") {
    show(note, answer.message || WORKING, false);
  } else if (answer.state === "invalid") {
    show(note, NOT_SIZED + answer.message, false);
  } else {
    show(note, "", false);
  }
}

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { fieldsReady } = await import(`./fieldsready.js${VERSION_QUERY}`);
await fieldsReady();

for (const note of document.querySelectorAll("[data-size-report]")) {
  const form = note.closest("form");
  let timer = null;
  let sent = 0;
  let pending = null;
  let retries = 0;
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
      if (mine === sent) {
        draw(note, answer);
        // The estimate is being worked out in the background: ask again until it is there.
        if (answer.state === "working" && retries++ < MAX_RETRIES) {
          clearTimeout(timer);
          timer = setTimeout(ask, (Number(answer.retry_seconds) || RETRY_MS / 1000) * 1000);
        } else {
          retries = 0;
        }
      }
    } catch (error) {
      if (mine === sent && error.name !== "AbortError") show(note, "", false);
    } finally {
      if (mine === sent) note.removeAttribute("aria-busy");
    }
  };
  const later = () => { clearTimeout(timer); retries = 0; timer = setTimeout(ask, WAIT_MS); };
  form.addEventListener("input", later);
  form.addEventListener("change", later);
  form.addEventListener("sl-change", later);
  note.hidden = true;
  ask();
}
