/* Streamlit glue: receives args on every Python rerun, diffs them against
 * what was last applied, and forwards only the changes to the long-lived
 * ViewerController. Events (atom clicks, pose selections) go back to Python
 * via setComponentValue with a strictly increasing nonce.
 */

import "molstar/build/viewer/molstar.css";
import "./style.css";

import { Streamlit, type RenderData } from "streamlit-component-lib";

import { PoseBar } from "./poseBar";
import { fingerprint } from "./synthetic";
import type { ComponentArgs, DepictionArg, ViewerEvent, ViewerEventBody } from "./types";
import { ViewerController } from "./viewer";

const DEFAULT_HEIGHT = 620;

const viewerEl = document.getElementById("viewer") as HTMLElement;
const poseBarEl = document.getElementById("poseBar") as HTMLElement;

let nonce = Date.now();
function emit(event: ViewerEventBody): void {
  const payload: ViewerEvent = { ...event, nonce: ++nonce };
  Streamlit.setComponentValue(payload);
}

const viewer = new ViewerController((atom) => emit({ type: "click", atom }));
const poseBar = new PoseBar(poseBarEl, {
  onSelect: (idx, light) => {
    viewer.enqueue(() => viewer.selectPose(idx, light));
    if (!light) emit({ type: "pose", index: idx });
  },
});

// Last-applied fingerprint per arg; deleting a key forces a re-apply.
const applied = new Map<string, string>();
function changed(key: string, value: unknown): boolean {
  const fp = fingerprint(value);
  if (applied.get(key) === fp) return false;
  applied.set(key, fp);
  return true;
}

let height = DEFAULT_HEIGHT;
function updateFrameHeight(): void {
  Streamlit.setFrameHeight(height + poseBar.height);
}

// -- 2D depiction overlay -----------------------------------------------------

let depictionUrl: string | null = null;

// Mounted inside .msp-viewport (not #viewer, which also holds Mol*'s side
// panels) and re-attached on every render, since Mol* can remount the
// viewport on layout changes.
function depictionEl(): HTMLElement {
  let el = document.getElementById("ligandDepiction");
  if (!el) {
    el = document.createElement("div");
    el.id = "ligandDepiction";
    el.className = "hidden";
    const name = document.createElement("div");
    name.className = "depname";
    const img = document.createElement("img");
    img.alt = "2D structure of the docked ligand";
    el.append(name, img);
  }
  const viewport = viewer.viewportElement(viewerEl);
  if (el.parentElement !== viewport) viewport.appendChild(el);
  return el;
}

// Rendered through an <img> blob URL, never innerHTML'd, so an SVG can't
// execute anything even if its escaping ever had a gap.
function showDepiction(depiction: DepictionArg | null): void {
  const el = depictionEl();
  if (depictionUrl) { URL.revokeObjectURL(depictionUrl); depictionUrl = null; }
  if (!depiction) { el.classList.add("hidden"); return; }
  depictionUrl = URL.createObjectURL(new Blob([depiction.svg], { type: "image/svg+xml" }));
  (el.querySelector("img") as HTMLImageElement).src = depictionUrl;
  (el.querySelector(".depname") as HTMLElement).textContent = depiction.name;
  el.classList.remove("hidden");
}

// -- render -------------------------------------------------------------------

function normalize(raw: Partial<ComponentArgs>): ComponentArgs {
  return {
    receptor: raw.receptor ?? null,
    het_ligands: raw.het_ligands ?? [],
    style: raw.style ?? null,
    box: raw.box ?? null,
    picks: raw.picks ?? [],
    constraints: raw.constraints ?? [],
    poses: raw.poses ?? null,
    depiction: raw.depiction ?? null,
    focus: raw.focus ?? null,
    height: raw.height ?? DEFAULT_HEIGHT,
  };
}

const ready = viewer.init(viewerEl);
let lastFocusNonce: number | null = null;

async function onRender(event: Event): Promise<void> {
  const args = normalize((event as CustomEvent<RenderData>).detail.args);

  if (args.height !== height || !viewerEl.style.height) {
    height = args.height;
    viewerEl.style.height = `${height}px`;
    updateFrameHeight();
  }
  await ready;

  // A new receptor wipes every overlay and pose, so force them all to be
  // re-applied against it.
  if (changed("receptor", args.receptor?.id ?? null)) {
    for (const key of ["style", "box", "picks", "constraints", "poses"]) applied.delete(key);
    poseBar.hide();
    viewer.enqueue(() => viewer.setReceptor(args.receptor));
  }
  if (changed("het_ligands", args.het_ligands)) viewer.setHetLigands(args.het_ligands);
  if (args.style && changed("style", args.style)) {
    const style = args.style;
    viewer.enqueue(() => viewer.setStyle(style));
    // hide_ligand changes which receptor copy the contacts merge in
    if (viewer.poseCount && poseBar.active >= 0) {
      viewer.enqueue(() => viewer.showInteractions(poseBar.active));
    }
  }
  if (changed("box", args.box)) viewer.enqueue(() => viewer.setBox(args.box));
  if (changed("picks", args.picks)) viewer.enqueue(() => viewer.setPicks(args.picks));
  if (changed("constraints", args.constraints)) viewer.enqueue(() => viewer.setConstraints(args.constraints));

  if (changed("poses", args.poses?.key ?? null)) {
    const poses = args.poses;
    poseBar.hide();
    await viewer.enqueue(() => viewer.setPoses(poses));
    if (poses && viewer.poseCount) {
      poseBar.show(poses.name, poses.scores, viewer.poseCount);
      poseBar.select(0, false);
    }
    updateFrameHeight();
  }

  if (changed("depiction", args.depiction)) showDepiction(args.depiction);
  else if (args.depiction) depictionEl();  // re-attach after a possible viewport remount

  if (args.focus && args.focus.nonce !== lastFocusNonce) {
    lastFocusNonce = args.focus.nonce;
    const focus = args.focus;
    viewer.enqueue(async () => viewer.focus(focus));
  }
}

Streamlit.events.addEventListener(Streamlit.RENDER_EVENT, (e) => {
  onRender(e).catch((err) => console.warn("render failed", err));
});
Streamlit.setComponentReady();
Streamlit.setFrameHeight(DEFAULT_HEIGHT);

// Read-only debug hook for the browser console.
(window as any).__dockingViewer = { viewer, poseBar };
