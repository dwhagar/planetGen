// html/static/objecttree.js
//
// The summary tree of object pages, the system page and the sector page (UX.172, web/lib/objecttree.py).
//
// The server draws a tree as nested lists: a row is the kind's mark, the name (a link to the object's page) and a summary
// line, and a row with children also has a toggle button and a list of them. Without this script every level is just
// listed and each name is a link. The script adds what script is for:
//
// - the toggle buttons: a row not drawn open (`data-open`) is closed until its toggle is pressed;
// - a row whose children are not on the page (its list holds a "pending" link to its own page) loads them when it is
//   opened, from the route in its `data-tree-url` (`GET /object/<ref>/tree`: `{html, next}`), and a level longer than a
//   page ends in a "Show more" button that reads the next page;
// - the keys of the WAI-ARIA tree pattern: Up and Down move between the visible rows, Right opens a row or goes to its
//   first child, Left closes it or goes to its parent, Home and End go to the first and last row. One row is in the
//   Tab order at a time (the names are the rows; the toggles are left to the mouse and these keys).
//
// A failed load says so in the row and tries again the next time the row is opened. Nothing here writes markup of its
// own except those messages and the "Show more" row: the rows are the server's.

const ROW = "a.tree-name, [data-tree-more]";

function reachable(row) {
  // Not inside the list of a closed row.
  return !row.closest("ul.tree-children[hidden]");
}

function rows(root) {
  return Array.from(root.querySelectorAll(ROW)).filter(reachable);
}

function levelOf(row) {
  let level = 1;
  for (let node = row.parentElement; node; node = node.parentElement) {
    if (node.classList && node.classList.contains("tree-children")) {
      level += 1;
    }
    if (node.hasAttribute && node.hasAttribute("data-tree-root")) {
      break;
    }
  }
  return level;
}

function toggleOf(node) {
  return node.querySelector(":scope > .tree-row > .tree-toggle");
}

function listOf(node) {
  return node.querySelector(":scope > ul.tree-children");
}

function isOpen(node) {
  return node.hasAttribute("data-open");
}

function setOpen(node, open) {
  node.toggleAttribute("data-open", open);
  node.setAttribute("aria-expanded", open ? "true" : "false");
  const toggle = toggleOf(node);
  if (toggle) {
    toggle.setAttribute("aria-expanded", open ? "true" : "false");
  }
  const list = listOf(node);
  if (list) {
    list.hidden = !open;
  }
}

function enhance(root, scope) {
  const container = scope || root;
  const nodes = Array.from(container.querySelectorAll("li.tree-node"));
  if (container.matches && container.matches("li.tree-node")) {
    nodes.push(container);
  }
  for (const node of nodes) {
    node.setAttribute("role", "treeitem");
    const row = node.querySelector(":scope > .tree-row a.tree-name");
    if (row) {
      node.setAttribute("aria-level", String(levelOf(row)));
    }
    const toggle = toggleOf(node);
    if (toggle) {
      toggle.hidden = false;
      toggle.tabIndex = -1;
      if (!node.dataset.enhanced) {
        node.dataset.enhanced = "1";
        setOpen(node, isOpen(node));
      }
    }
  }
  for (const group of container.querySelectorAll("ul.tree-children")) {
    group.setAttribute("role", "group");
  }
  for (const pending of container.querySelectorAll("li.tree-more, li.tree-pending")) {
    pending.setAttribute("role", "none");
  }
  for (const row of container.querySelectorAll(ROW)) {
    row.tabIndex = -1;
  }
  if (!root.querySelector(ROW + "[tabindex='0']")) {
    const first = root.querySelector(ROW);
    if (first) {
      first.tabIndex = 0;
    }
  }
}

function focusRow(root, row) {
  for (const other of root.querySelectorAll(ROW + "[tabindex='0']")) {
    other.tabIndex = -1;
  }
  row.tabIndex = 0;
  row.focus();
}

function treeUrl(element) {
  const holder = element.closest("[data-tree-url]");
  return holder ? holder.dataset.treeUrl : null;
}

function fetchLevel(url, after) {
  const address = after ? url + (url.includes("?") ? "&" : "?") + "after=" + encodeURIComponent(after) : url;
  return fetch(address, { credentials: "same-origin", headers: { Accept: "application/json" } }).then(function (response) {
    if (!response.ok) {
      throw new Error("HTTP " + response.status);
    }
    return response.json();
  });
}

function htmlNodes(html) {
  const template = document.createElement("template");
  template.innerHTML = html;
  return template.content;
}

function moreRow(next) {
  const item = document.createElement("li");
  item.className = "tree-more";
  const button = document.createElement("button");
  button.type = "button";
  button.className = "btn btn-small";
  button.dataset.treeMore = "";
  button.dataset.next = next;
  button.textContent = "Show more";
  item.appendChild(button);
  return item;
}

function say(list, text) {
  const note = document.createElement("li");
  note.className = "tree-pending tree-error";
  note.setAttribute("role", "none");
  note.textContent = text;
  list.appendChild(note);
}

function loadInto(root, list, url, after, before) {
  return fetchLevel(url, after).then(function (answer) {
    const fragment = htmlNodes(answer.html);
    const added = Array.from(fragment.children);
    list.insertBefore(fragment, before || null);
    if (answer.next) {
      list.appendChild(moreRow(answer.next));
    }
    for (const node of added) {
      enhance(root, node);
    }
    enhance(root, list);
  });
}

function openRow(root, node) {
  setOpen(node, true);
  const list = listOf(node);
  const pending = list && list.querySelector(":scope > li.tree-pending");
  if (!pending || node.dataset.treeLoading) {
    return;
  }
  node.dataset.treeLoading = "1";
  loadInto(root, list, node.dataset.treeUrl, null, null).then(function () {
    pending.remove();
  }, function () {
    pending.remove();
    say(list, "Could not load this level. Close the row and open it again to retry.");
  }).then(function () {
    delete node.dataset.treeLoading;
  });
}

function toggleRow(root, node) {
  if (isOpen(node)) {
    setOpen(node, false);
  } else {
    openRow(root, node);
  }
}

function showMore(root, button) {
  const item = button.closest("li");
  const list = item.parentElement;
  const url = treeUrl(list);
  if (!url || button.disabled) {
    return;
  }
  button.disabled = true;
  const focused = document.activeElement === button;
  const after = button.dataset.next;
  loadInto(root, list, url, after, item).then(function () {
    const first = item.previousElementSibling && item.previousElementSibling.querySelector(ROW);
    item.remove();
    if (focused && first) {
      focusRow(root, first);
    }
  }, function () {
    button.disabled = false;
    say(list, "Could not load more. Press Show more to try again.");
  });
}

function keydown(root, event) {
  const row = event.target.closest(ROW);
  if (!row || event.altKey || event.ctrlKey || event.metaKey) {
    return;
  }
  const visible = rows(root);
  const index = visible.indexOf(row);
  const node = row.closest("li.tree-node");
  const branch = node && toggleOf(node) ? node : null;
  let target = null;
  switch (event.key) {
    case "ArrowDown":
      target = visible[index + 1];
      break;
    case "ArrowUp":
      target = visible[index - 1];
      break;
    case "Home":
      target = visible[0];
      break;
    case "End":
      target = visible[visible.length - 1];
      break;
    case "ArrowRight":
      if (branch && !isOpen(branch)) {
        openRow(root, branch);
        event.preventDefault();
        return;
      }
      target = branch ? listOf(branch).querySelector(ROW) : null;
      break;
    case "ArrowLeft":
      if (branch && isOpen(branch)) {
        setOpen(branch, false);
        event.preventDefault();
        return;
      }
      target = (function () {
        const parent = (node || row.closest("li")).parentElement.closest("li.tree-node");
        return parent ? parent.querySelector(":scope > .tree-row a.tree-name") : null;
      })();
      break;
    default:
      return;
  }
  event.preventDefault();
  if (target) {
    focusRow(root, target);
  }
}

function setup(root) {
  enhance(root);
  root.setAttribute("role", "tree");
  root.addEventListener("click", function (event) {
    const button = event.target.closest("[data-tree-more]");
    if (button) {
      showMore(root, button);
      return;
    }
    const toggle = event.target.closest(".tree-toggle");
    if (toggle) {
      toggleRow(root, toggle.closest("li.tree-node"));
    }
  });
  root.addEventListener("keydown", function (event) {
    keydown(root, event);
  });
  root.addEventListener("focusin", function (event) {
    const row = event.target.closest(ROW);
    if (row && row.tabIndex !== 0) {
      for (const other of root.querySelectorAll(ROW + "[tabindex='0']")) {
        other.tabIndex = -1;
      }
      row.tabIndex = 0;
    }
  });
}

for (const root of document.querySelectorAll("[data-tree-root]")) {
  setup(root);
}
