// static/joblog.js
//
// The job output box (ADM.22): a terminal (Xterm.js, vendor/xterm/) that
// takes the job's raw output as it streams, so carriage-return progress
// redraws and colour codes show as they would in a console. If the terminal
// can't start, it falls back to writing the text into the page's <pre>.
//
//   const log = await createJobLog(preElement);
//   log.write(text);   // more output
//
// Imported by generatejobs.js with this module's own `?v=` query, as
// systemmap.js explains.

const VERSION_QUERY = new URL(import.meta.url).search;

// The text a plain <pre> shows: each line as its last carriage return left it.
function collapse(text) {
  return text.split("\n").map((line) => line.replace(/\r$/, "").split("\r").pop()).join("\n");
}

function plainLog(pre) {
  // The stream replays the end of the log, so the server-rendered text is replaced, not added to.
  let buffer = "";
  pre.textContent = "";
  return {
    terminal: false,
    write(text) {
      buffer += text.replace(/\r\n/g, "\n");
      const atBottom = pre.scrollHeight - pre.scrollTop - pre.clientHeight < 24;
      pre.textContent = collapse(buffer);
      if (atBottom) pre.scrollTop = pre.scrollHeight;
    },
  };
}

export async function createJobLog(pre) {
  // The terminal starts empty and is filled by the stream, which replays the end of the log.
  try {
    const [{ Terminal }, { FitAddon }] = await Promise.all([
      import(`./vendor/xterm/xterm.mjs${VERSION_QUERY}`),
      import(`./vendor/xterm/addon-fit.mjs${VERSION_QUERY}`),
    ]);
    const css = document.createElement("link");
    css.rel = "stylesheet";
    css.href = new URL(`./vendor/xterm/xterm.css${VERSION_QUERY}`, import.meta.url).href;
    document.head.appendChild(css);

    const host = document.createElement("div");
    host.className = "job-terminal";
    host.setAttribute("role", "log");
    host.setAttribute("aria-label", "Job output");
    pre.after(host);
    pre.hidden = true;

    const style = getComputedStyle(document.documentElement);
    const color = (name, fallback) => style.getPropertyValue(name).trim() || fallback;
    const terminal = new Terminal({
      convertEol: true,
      disableStdin: true,
      scrollback: 20000,
      fontSize: 13,
      fontFamily: color("--font-mono", "monospace"),
      cursorInactiveStyle: "none",
      theme: { background: color("--bg", "#14151e"), foreground: color("--text", "#e7e9f3") },
    });
    const fit = new FitAddon();
    terminal.loadAddon(fit);
    terminal.open(host);
    fit.fit();
    new ResizeObserver(() => fit.fit()).observe(host);
    return { terminal: true, write: (text) => terminal.write(text) };
  } catch (err) {
    const stray = document.querySelector(".job-terminal");
    if (stray) stray.remove();
    pre.hidden = false;
    return plainLog(pre);
  }
}
