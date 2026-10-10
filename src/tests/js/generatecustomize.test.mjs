// tests/js/generatecustomize.test.mjs -- static/generatecustomize.js, the
// Generate page's Customize window (ADM.28): settings groups become tabs,
// one shows at a time, the inputs stay in the form, and a rejected value
// opens the window on its tab.

import { test } from "node:test";
import assert from "node:assert/strict";

import { FakeEvent, h, installDom, STATIC_URL } from "./fakedom.mjs";

installDom("http://localhost/");

function build() {
  const input = h("input", { name: "prevalence_comets" });
  const second = h("input", { name: "directive_systems" });
  const groups = [
    h("fieldset", { "data-customize-tab": "Prevalence" }, [h("legend", {}, "Prevalence"), input]),
    h("fieldset", { "data-customize-tab": "Override" }, [h("legend", {}, "Override"), second]),
  ];
  const submit = h("button", { type: "submit" }, "Generate");
  const actions = h("div", { class: "search-actions" }, [submit]);
  const form = h("form", { "data-customize": "" }, [...groups, actions]);
  document.body.appendChild(form);
  return { form, groups, input, second, actions, submit };
}

const mod = await import(new URL("generatecustomize.js", STATIC_URL).href);

test("a form's settings groups become tabs in a window opened by a Customize button", () => {
  const { form, groups, actions } = build();
  const built = mod.enhance(form);
  assert.ok(built);
  assert.equal(built.button.textContent, "Customize");
  assert.equal(actions.children[0], built.button, "the button leads the form's own actions");
  assert.deepEqual(built.dialog.querySelectorAll(".customize-tab").map((tab) => tab.textContent), ["Prevalence", "Override"]);
  assert.ok(groups.every((group) => built.dialog.contains(group)), "the groups moved into the window");
  assert.ok(form.contains(built.dialog), "the window stays inside the form so its inputs still submit");
});

test("only the selected tab's group shows, and clicking a tab switches", () => {
  const { form, groups } = build();
  const built = mod.enhance(form);
  assert.deepEqual(groups.map((group) => group.hidden), [false, true]);
  built.dialog.querySelectorAll(".customize-tab")[1].dispatchEvent(new FakeEvent("click"));
  assert.deepEqual(groups.map((group) => group.hidden), [true, false]);
});

test("the Customize button opens the window and Done closes it", () => {
  const { form } = build();
  const built = mod.enhance(form);
  built.button.dispatchEvent(new FakeEvent("click"));
  assert.ok(built.dialog.hasAttribute("open"));
  built.dialog.querySelectorAll("button").find((b) => b.textContent === "Done").dispatchEvent(new FakeEvent("click"));
  assert.ok(!built.dialog.hasAttribute("open"));
});

test("a rejected value opens the window on the tab that holds it", () => {
  const { form, groups, second } = build();
  const built = mod.enhance(form);
  second.dispatchEvent(new FakeEvent("invalid", { bubbles: false }));
  assert.deepEqual(groups.map((group) => group.hidden), [true, false]);
  assert.ok(built.dialog.hasAttribute("open"));
});

test("a form with no settings groups is left alone", () => {
  const form = h("form", { "data-customize": "" }, [h("div", { class: "search-actions" })]);
  assert.equal(mod.enhance(form), null);
});
