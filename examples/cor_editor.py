"""Manual placement of centres of rotation (CoR) on Hologic whole-body DXA scans.

Usage
-----
    python examples/cor_editor.py SCAN1.dcm SCAN2.dcm ... --cor-dir data_private/cor

For each scan ``<stem>.dcm`` the markers are saved to ``<cor-dir>/cor_manual_<stem>.json``
as ``{marker: [row, col]}`` (whole-body grid, row 0 = patient right, col 0 = top of
head).  The stem is the file name as given, never a DICOM tag value.  The JSON files are
read by ``hologic-dxa segment-profiles --cor``.

Opening a scan reloads its saved markers and selects the first missing one; after each
click the next *missing* marker is selected, so only new markers need placing.

Mouse and keys
--------------
    left click  place the selected marker      right click / Delete  erase it
    Tab / Shift+Tab  next / previous marker    Enter  save
    Left / Right arrow  previous / next scan

Image orientation: the left edge of the image is the patient's RIGHT (letter R).
EXPERIMENTAL: display and marker placement only, no measurement is made here.
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
import tkinter as tk
from pathlib import Path
from tkinter import ttk

import matplotlib
import numpy as np
import pydicom

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from hologic_dxa.segments import (
    MARKER_NAMES,
    build_segments,
    load_cor_json,
    load_whole_body_mu,
)

GRID_ROWS, GRID_COLS = 106, 150

MIDLINE_COLOR = "#55aaff"
HEAD_COLOR = "#44ee77"
LEFT_COLORS = {
    "glenohumeral": "#ff4444", "elbow": "#dd2222", "wrist": "#ff7777", "hand_end": "#ffaaaa",
    "hip": "#cc44ff", "knee": "#aa22dd", "ankle": "#881199", "foot_end": "#bb66dd",
}  # fmt: skip
RIGHT_COLORS = {
    "glenohumeral": "#ffaa33", "elbow": "#dd8822", "wrist": "#ffcc66", "hand_end": "#ffe0a0",
    "hip": "#ff44cc", "knee": "#dd22aa", "ankle": "#991188", "foot_end": "#dd66bb",
}  # fmt: skip

ENHANCEMENTS = [
    "tag1002_ch1", "tag1002_ch1_raw", "tag1002_ch1_log",
    "tag1002_ch2", "tag1002_ch2_raw", "tag1002_ch2_log",
    "clahe_sq", "clahe", "dual_sub", "gamma", "log", "raw",
]  # fmt: skip


def marker_color(name: str) -> str:
    if name in ("head_top", "atlanto_occipital"):
        return HEAD_COLOR
    if name.endswith("_L"):
        return LEFT_COLORS.get(name[:-2], "#ffffff")
    if name.endswith("_R"):
        return RIGHT_COLORS.get(name[:-2], "#ffffff")
    return MIDLINE_COLOR


def parse_hd_channels(ds: pydicom.Dataset) -> dict[str, np.ndarray]:
    """High-resolution float64 images of tag (0023,1002): ch1/ch2/ch3, shape (rows, cols).

    Rows run head to foot, columns are lateral (unlike the 106x150 grid).
    """
    if (0x0023, 0x1002) not in ds:
        return {}
    data = bytes(ds[0x0023, 0x1002].value)
    types = {370: "ch1", 372: "ch2", 374: "ch3"}
    channels: dict[str, np.ndarray] = {}
    pos = 0
    while pos + 6 <= len(data):
        rtype = struct.unpack_from("<H", data, pos)[0]
        length = struct.unpack_from("<I", data, pos + 2)[0]
        if length < 6 or pos + length > len(data):
            break
        payload = data[pos + 6 : pos + length]
        if rtype in types and len(payload) >= 40:
            rows = struct.unpack_from("<I", payload, 0)[0]
            cols = struct.unpack_from("<I", payload, 4)[0]
            if len(payload) >= 40 + rows * cols * 8:
                raw = np.frombuffer(payload[40 : 40 + rows * cols * 8], dtype="<f8")
                channels[types[rtype]] = raw.reshape(rows, cols).copy()
        if rtype == 0:
            break
        pos += length
    return channels


def _normalise(a: np.ndarray, lo_pct: float = 0, hi_pct: float = 99) -> np.ndarray:
    body = a[a > 0]
    if body.size == 0:
        return a
    lo, hi = np.percentile(body, lo_pct), np.percentile(body, hi_pct)
    return np.clip((a - lo) / (hi - lo + 1e-12), 0, 1)


def _clahe(img: np.ndarray, clip: float) -> np.ndarray:
    try:
        from skimage.exposure import equalize_adapthist

        return equalize_adapthist(img.astype(np.float32), clip_limit=clip, nbins=256)
    except ImportError:
        return img


def enhance(
    mode: str,
    mu_h: np.ndarray,
    mu_l: np.ndarray,
    hd: dict[str, np.ndarray],
    row_size_m: float,
    col_size_m: float,
) -> tuple[np.ndarray, bool]:
    """Return (image, is_hd): grid modes are (lateral, head-to-foot), HD is transposed."""
    if mode.startswith("tag1002") and ("ch1" if "ch1" in mode else "ch2") in hd:
        img = hd["ch1" if "ch1" in mode else "ch2"].astype(np.float64)
        body = img != 0.0
        if body.any():
            lo, hi = np.percentile(img[body], 1), np.percentile(img[body], 99)
            img = np.clip((img - lo) / (hi - lo + 1e-12), 0, 1)
            img[~body] = 0.0
        if mode.endswith("_raw"):
            return img, True
        if mode.endswith("_log"):
            return np.log1p(img * 10) / np.log1p(10), True
        return _clahe(img, 0.02), True

    if mode in ("clahe_sq", "tag1002_ch1", "tag1002_ch2"):
        from scipy.ndimage import zoom

        up = zoom(mu_h, (4.0, 4.0 * col_size_m / row_size_m), order=3)
        return _clahe(_normalise(up, 0, 99), 0.02), False
    if mode == "clahe":
        return _clahe(_normalise(mu_h), 0.03), False
    if mode == "dual_sub":
        return _normalise(np.clip(mu_h - 0.4 * mu_l, 0, None)), False
    if mode == "gamma":
        return np.power(_normalise(mu_h), 0.5), False
    if mode == "log":
        return _normalise(np.log1p(mu_h), 1, 99), False
    return _normalise(mu_h), False


class CorEditor(tk.Tk):
    def __init__(self, scans: list[Path], cor_dir: Path) -> None:
        super().__init__()
        self.title("DXA CoR editor")
        self.configure(bg="#1e1e2e")
        self.scans, self.cor_dir, self.scan_idx = scans, cor_dir, 0
        self.mu_h = self.mu_l = None
        self.row_size_m = self.col_size_m = 0.0
        self.hd: dict[str, np.ndarray] = {}
        self.joints: dict[str, tuple[float, float]] = {}
        self.sel = MARKER_NAMES[0]
        self.plo, self.phi = tk.DoubleVar(value=0.0), tk.DoubleVar(value=99.0)
        self.cmap = tk.StringVar(value="hot")
        self.mode = tk.StringVar(value="tag1002_ch1")
        self.status = tk.StringVar(value="Ready")
        self._build_ui()
        self._load_scan(0)
        self.bind("<Key>", self._on_key)

    def _build_ui(self) -> None:
        bg, bg2, fg, acc = "#1e1e2e", "#181825", "#cdd6f4", "#cba6f7"
        top = tk.Frame(self, bg=bg)
        top.pack(fill="x", padx=8, pady=4)
        tk.Button(top, text="<", command=lambda: self._go(-1)).pack(side="left")
        self.lbl_scan = tk.Label(top, bg=bg, fg=acc, font=("Segoe UI", 12, "bold"))
        self.lbl_scan.pack(side="left", padx=8)
        tk.Button(top, text=">", command=lambda: self._go(1)).pack(side="left")
        tk.Button(top, text="Save", bg="#a6e3a1", command=self._save).pack(side="right")

        body = tk.Frame(self, bg=bg)
        body.pack(fill="both", expand=True, padx=8)
        self.fig = Figure(figsize=(5.5, 12), facecolor=bg2)
        self.ax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
        self.canvas = FigureCanvasTkAgg(self.fig, master=body)
        self.canvas.get_tk_widget().pack(side="left", fill="both", expand=True)
        self.canvas.mpl_connect("button_press_event", self._on_click)

        right = tk.Frame(body, bg=bg2, width=230)
        right.pack(side="right", fill="y", padx=(6, 0))
        right.pack_propagate(False)
        self.lb = tk.Listbox(
            right, bg="#313244", fg=fg, exportselection=False, activestyle="none",
            font=("Consolas", 9), selectbackground=acc, selectforeground=bg,
        )  # fmt: skip
        self.lb.pack(fill="both", expand=True, padx=6, pady=6)
        self.lb.bind("<<ListboxSelect>>", self._on_lb)
        self.lbl_info = tk.Label(right, bg=bg2, fg=fg, font=("Consolas", 9), justify="left")
        self.lbl_info.pack(anchor="w", padx=8)
        for label, var in (("Colormap", self.cmap), ("Processing", self.mode)):
            row = tk.Frame(right, bg=bg2)
            row.pack(fill="x", padx=8, pady=2)
            tk.Label(row, text=label, bg=bg2, fg=fg).pack(side="left")
            is_cmap = var is self.cmap
            values = ["hot", "bone", "inferno", "plasma", "gray"] if is_cmap else ENHANCEMENTS
            cb = ttk.Combobox(row, textvariable=var, values=values, state="readonly", width=14)
            cb.pack(side="right")
            cb.bind("<<ComboboxSelected>>", lambda _e: self._draw())
        for label, var, lo, hi in (("Floor %", self.plo, 0, 60), ("Ceiling %", self.phi, 40, 100)):
            tk.Label(right, text=label, bg=bg2, fg=fg).pack(anchor="w", padx=8)
            tk.Scale(
                right, variable=var, from_=lo, to=hi, orient="horizontal", bg=bg2, fg=fg,
                highlightthickness=0, command=lambda _v: self._draw(),
            ).pack(fill="x", padx=8)  # fmt: skip
        tk.Label(self, textvariable=self.status, bg=bg2, fg="#6c7086", anchor="w").pack(fill="x")

    # ── data ────────────────────────────────────────────────────────────────
    def _cor_path(self) -> Path:
        return self.cor_dir / f"cor_manual_{self.scans[self.scan_idx].stem}.json"

    def _load_scan(self, idx: int) -> None:
        path = self.scans[idx]
        self.lbl_scan.config(text=f"[{idx + 1}/{len(self.scans)}]  {path.stem}")
        try:
            ds = pydicom.dcmread(str(path), force=True)
            mu = load_whole_body_mu(ds)
        except ValueError as exc:
            self.status.set(str(exc))
            return
        self.mu_h, self.mu_l = mu.mu_h, mu.mu_l
        self.row_size_m, self.col_size_m = mu.row_size_m, mu.col_size_m
        self.hd = parse_hd_channels(ds)
        self.joints = {}
        if self._cor_path().exists():
            self.joints = load_cor_json(self._cor_path(), (GRID_ROWS, GRID_COLS))
        missing = [m for m in MARKER_NAMES if m not in self.joints]
        self.status.set(
            f"{len(self.joints)} markers loaded, {len(missing)} to place"
            + (f" (HD {next(iter(self.hd.values())).shape})" if self.hd else "")
        )
        self._select(missing[0] if missing else MARKER_NAMES[0])

    def _go(self, step: int) -> None:
        if 0 <= self.scan_idx + step < len(self.scans):
            self.scan_idx += step
            self._load_scan(self.scan_idx)

    def _save(self) -> None:
        self.cor_dir.mkdir(parents=True, exist_ok=True)
        data = {k: list(v) for k, v in self.joints.items()}
        self._cor_path().write_text(json.dumps(data, indent=2), encoding="utf-8")
        self.status.set(f"Saved {self._cor_path().name}")

    # ── selection ───────────────────────────────────────────────────────────
    def _select(self, name: str) -> None:
        self.sel = name
        self.lb.delete(0, "end")
        for i, m in enumerate(MARKER_NAMES):
            self.lb.insert("end", f"{'x' if m in self.joints else '.'} {m}")
            self.lb.itemconfig(i, fg=marker_color(m))
        i = MARKER_NAMES.index(name)
        self.lb.selection_set(i)
        self.lb.see(i)
        if name in self.joints:
            r, c = self.joints[name]
            self.lbl_info.config(text=f"{name}\nrow {r:.1f}  col {c:.1f}")
        else:
            self.lbl_info.config(text=f"{name}\nnot placed")
        self._draw()

    def _next_missing(self, after: str) -> str:
        n, i = len(MARKER_NAMES), MARKER_NAMES.index(after)
        for k in range(1, n + 1):
            cand = MARKER_NAMES[(i + k) % n]
            if cand not in self.joints:
                return cand
        return MARKER_NAMES[(i + 1) % n]

    def _on_lb(self, _event: object) -> None:
        sel = self.lb.curselection()
        if sel:
            self._select(MARKER_NAMES[sel[0]])

    def _on_key(self, event: tk.Event) -> None:
        k, i = event.keysym, MARKER_NAMES.index(self.sel)
        if k == "Return":
            self._save()
        elif k in ("Tab", "ISO_Left_Tab"):
            step = -1 if k == "ISO_Left_Tab" else 1
            self._select(MARKER_NAMES[(i + step) % len(MARKER_NAMES)])
        elif k == "Left":
            self._go(-1)
        elif k == "Right":
            self._go(1)
        elif k == "Delete":
            self.joints.pop(self.sel, None)
            self._select(self.sel)

    def _on_click(self, event: object) -> None:
        if event.inaxes != self.ax or event.xdata is None:
            return
        if event.button == 3:
            self.joints.pop(self.sel, None)
            self._select(self.sel)
            return
        row = float(np.clip(event.xdata, 0, GRID_ROWS - 1e-3))
        col = float(np.clip(event.ydata, 0, GRID_COLS - 1e-3))
        self.joints[self.sel] = (row, col)
        self._select(self._next_missing(self.sel))

    # ── drawing ─────────────────────────────────────────────────────────────
    def _draw(self) -> None:
        ax = self.ax
        ax.cla()
        if self.mu_h is not None:
            img, is_hd = enhance(
                self.mode.get(), self.mu_h, self.mu_l, self.hd, self.row_size_m, self.col_size_m
            )
            shown = np.clip(img if is_hd else img.T, 0, None)
            body = shown[shown > 0]
            lo_p, hi_p = self.plo.get(), self.phi.get()
            vmin, vmax = (
                (np.percentile(body, lo_p), np.percentile(body, max(hi_p, lo_p + 1)))
                if body.size
                else (0.0, 1.0)
            )
            ax.imshow(
                shown, origin="upper", cmap=self.cmap.get(), aspect="auto",
                extent=[0, GRID_ROWS, GRID_COLS, 0], vmin=vmin, vmax=vmax,
            )  # fmt: skip
        for _definition, prox, dist in build_segments(self.joints):
            ax.plot([prox[0], dist[0]], [prox[1], dist[1]], "-", color="white", alpha=0.35, lw=1.2)
        for name, (row, col) in self.joints.items():
            is_sel = name == self.sel
            ax.plot(
                row, col, "o", color=marker_color(name), markersize=13 if is_sel else 7,
                markeredgecolor="white", markeredgewidth=2.5 if is_sel else 1.0, zorder=10,
            )  # fmt: skip
            ax.text(
                row + 2, col, name.replace("_", " "), color=marker_color(name), fontsize=6.5,
                va="center", alpha=1.0 if is_sel else 0.75,
            )  # fmt: skip
        label = dict(fontsize=9, fontweight="bold", va="center", color="white")
        ax.text(2, 75, "R", ha="left", **label)
        ax.text(104, 75, "L", ha="right", **label)
        ax.set_xlim(0, GRID_ROWS)
        ax.set_ylim(GRID_COLS, 0)
        ax.axis("off")
        self.canvas.draw_idle()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("scans", nargs="+", type=Path, help="Hologic whole-body .dcm files")
    parser.add_argument("--cor-dir", type=Path, default=Path("data_private/cor"))
    args = parser.parse_args()
    missing = [p for p in args.scans if not p.is_file()]
    if missing:
        sys.exit(f"File not found: {', '.join(str(p) for p in missing)}")
    app = CorEditor(args.scans, args.cor_dir)
    app.geometry("950x920")
    app.mainloop()


if __name__ == "__main__":
    main()
