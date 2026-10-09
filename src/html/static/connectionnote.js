// html/static/connectionnote.js
//
// When to tell the admin a live stream was lost (ADM.40). A job's log
// stream (`/admin/generate/jobs/<id>/stream`) ends on purpose every
// STREAM_SECONDS to stay under the web server's request timeout, and the
// browser's EventSource reports every such end as an "error" before it
// reconnects with the last event id, so nothing is lost. Saying "lost
// connection" at each of those was noise. The note now appears only when
// the stream has stayed down for `graceMs`, and clears when it is back.

/**
 * @param {(text: string) => void} show Shows (or, with "", clears) the note.
 * @param {number} graceMs How long a stream may stay down before the note shows.
 * @param {{setTimeout: Function, clearTimeout: Function}} timers For tests.
 * @returns {{lost: () => void, restored: () => void}}
 */
export function quietReconnect(show, graceMs, timers = globalThis) {
  let pending = null;
  const cancel = () => {
    if (pending !== null) {
      timers.clearTimeout(pending);
      pending = null;
    }
  };
  return {
    lost() {
      if (pending === null) {
        pending = timers.setTimeout(() => {
          pending = null;
          show("Connection lost; reconnecting...");
        }, graceMs);
      }
    },
    restored() {
      cancel();
      show("");
    },
  };
}
