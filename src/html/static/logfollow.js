// html/static/logfollow.js
//
// Admin log boxes that follow the end of the log: a `<pre data-follow-end>` (the Update page's update.log, a queue
// run's log) is scrolled to its last line when the page loads and whenever its text changes. The reader who scrolls
// up pauses the following, so the lines they are reading stay put; scrolling back to the end, or pressing End, picks
// it up again. The Update page reloads itself while an update runs, so the pause (and the place read) is kept in
// sessionStorage across the reload. A separate file because the site's Content-Security-Policy allows no inline
// script. Without script the box shows the log from its top, as before.

const NEAR_END_PX = 24;

function atEnd(pre) {
  return pre.scrollHeight - pre.scrollTop - pre.clientHeight < NEAR_END_PX;
}

function keyFor(pre, index) {
  return `logfollow:${location.pathname}:${pre.id || index}`;
}

function remember(key, value) {
  try {
    if (value === null) {
      sessionStorage.removeItem(key);
    } else {
      sessionStorage.setItem(key, String(value));
    }
  } catch (err) {
    // Storage blocked: the reload then follows the end, which is the default.
  }
}

function recall(key) {
  try {
    const text = sessionStorage.getItem(key);
    return text === null ? null : Number(text);
  } catch (err) {
    return null;
  }
}

function follow(pre, index) {
  const key = keyFor(pre, index);
  const saved = recall(key);
  let following = saved === null || Number.isNaN(saved);
  const toEnd = () => { pre.scrollTop = pre.scrollHeight; };
  if (following) {
    toEnd();
  } else {
    pre.scrollTop = saved;
  }
  pre.dataset.following = following ? "true" : "false";
  pre.addEventListener("scroll", () => {
    following = atEnd(pre);
    pre.dataset.following = following ? "true" : "false";
    remember(key, following ? null : pre.scrollTop);
  });
  new MutationObserver(() => {
    if (following) {
      toEnd();
    }
  }).observe(pre, { childList: true, characterData: true, subtree: true });
  // Fonts and the layout settle after load and can lengthen the box.
  window.addEventListener("load", () => {
    if (following) {
      toEnd();
    }
  });
}

document.querySelectorAll("pre[data-follow-end]").forEach(follow);
