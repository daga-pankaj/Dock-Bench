/* Pose bar under the viewer: ⏮ ▶/❚❚ ⏭, one button per pose, "all".
 *
 * Owns the active pose index and playback timer. `onSelect(idx, light)` is
 * called for every change; `light` is true for playback steps (no contacts,
 * no event to Python), and the pose playback stopped on is re-selected with
 * light=false.
 */

export const POSE_PLAY_INTERVAL_MS = 1200;

export interface PoseBarHandlers {
  onSelect: (idx: number, light: boolean) => void;
}

// Next active pose for a step of `delta`, wrapping; from "all" (-1) it lands
// on the first (forward) or last (back) pose.
export function stepIndex(active: number, delta: number, n: number): number {
  if (active < 0) return delta > 0 ? 0 : n - 1;
  return (active + delta + n) % n;
}

export class PoseBar {
  active = 0;
  private count = 0;
  private timer: ReturnType<typeof setInterval> | null = null;
  private playBtn: HTMLButtonElement | null = null;

  constructor(private el: HTMLElement, private handlers: PoseBarHandlers) {}

  get height(): number {
    return this.el.classList.contains("hidden") ? 0 : this.el.offsetHeight + 10;
  }

  hide(): void {
    this.stop(false);
    this.count = 0;
    this.el.innerHTML = "";
    this.el.classList.add("hidden");
  }

  show(name: string, scores: (number | null)[], count: number): void {
    this.stop(false);
    this.count = count;
    this.el.innerHTML = "";
    const title = document.createElement("span");
    title.className = "title";
    title.textContent = `${name} - poses:`;
    this.el.appendChild(title);

    const stepBtn = (text: string, tip: string, onclick: () => void) => {
      const b = document.createElement("button");
      b.className = "pose-play";
      b.textContent = text;
      b.title = tip;
      b.onclick = onclick;
      this.el.appendChild(b);
      return b;
    };
    stepBtn("⏮", "previous pose", () => { this.stop(false); this.step(-1, false); });
    this.playBtn = stepBtn("▶", "play through the docked poses", () => this.togglePlayback());
    stepBtn("⏭", "next pose", () => { this.stop(false); this.step(1, false); });

    for (let i = 0; i < count; i++) {
      const b = document.createElement("button");
      const score = scores[i];
      b.className = "pose";
      b.textContent = `${i + 1}${typeof score === "number" ? `  ${score.toFixed(2)}` : ""}`;
      b.onclick = () => { this.stop(false); this.select(i, false); };
      this.el.appendChild(b);
    }
    const all = document.createElement("button");
    all.className = "pose";
    all.textContent = "all";
    all.onclick = () => { this.stop(false); this.select(-1, false); };
    this.el.appendChild(all);
    this.el.classList.remove("hidden");
  }

  select(idx: number, light: boolean): void {
    this.active = idx;
    this.el.querySelectorAll(".pose").forEach((b, i) => {
      b.classList.toggle("active", i === (idx === -1 ? this.count : idx));
    });
    this.handlers.onSelect(idx, light);
  }

  private step(delta: number, light: boolean): void {
    if (!this.count) return;
    this.select(stepIndex(this.active, delta, this.count), light);
  }

  private togglePlayback(): void {
    if (this.timer) { this.stop(true); return; }
    if (this.count < 2) return;
    const step = () => this.step(1, true);
    step();
    this.timer = setInterval(step, POSE_PLAY_INTERVAL_MS);
    if (this.playBtn) {
      this.playBtn.textContent = "❚❚";
      this.playBtn.title = "pause";
      this.playBtn.classList.add("active");
    }
  }

  /** `refresh` re-selects the pose playback stopped on, with contacts. */
  stop(refresh: boolean): void {
    if (!this.timer) return;
    clearInterval(this.timer);
    this.timer = null;
    if (this.playBtn) {
      this.playBtn.textContent = "▶";
      this.playBtn.title = "play through the docked poses";
      this.playBtn.classList.remove("active");
    }
    if (refresh) this.select(this.active, false);
  }
}
