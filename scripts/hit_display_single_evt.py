#!/usr/bin/env python3
"""
MuCol hit display — tracker, calorimeter, or combined.

Reads ROOT files produced by slcio_to_root.py. Pass --tracker and/or --calo
to select which systems to include; at least one is required.

  --tracker only  →  R–Z (tracker scale), X–Y (full), X–Y (inner r < 620 mm)
  --calo    only  →  R–Z (calo scale),    X–Y (full), X–Y (ECAL  r < 2200 mm)
  both            →  R–Z (full scale),    X–Y (full), X–Y (tracker r < 620 mm)

Calorimeter marker area scales as log(energy); tracker markers are fixed.

Usage:
    python hit_display_single_evt.py --tracker                            # default file, event 1
    python hit_display_single_evt.py --calo                               # default file, event 1
    python hit_display_single_evt.py --tracker --calo                     # combined, both defaults
    python hit_display_single_evt.py --tracker hits.root -e 3 -o out.png
    python hit_display_single_evt.py --tracker hits.root --calo calo_hits.root -e 2
"""

import os
import argparse

import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_TRACKER = "/lustre/cmsdata/fnardi/mucoll_tutorial/pgun/pgun_tuple.root"
DEFAULT_CALO    = "/lustre/cmsdata/fnardi/mucoll_tutorial/pgun/pgun_calo.root"

# Tracker: blue (VXD), green (IT), teal (OT)
TRK_LABELS = ["VXD Barrel", "VXD Endcap", "IT Barrel", "IT Endcap", "OT Barrel", "OT Endcap"]
TRK_COLORS = ["#1f77b4", "#6baed6", "#2ca02c", "#98df8a", "#17becf", "#9edae5"]

# Calorimeter: red/orange (ECAL), purple (HCAL)
CAL_LABELS = ["ECAL Barrel", "ECAL Endcap", "HCAL Barrel", "HCAL Endcap"]
CAL_COLORS = ["#d62728", "#ff7f0e", "#9467bd", "#c5b0d5"]

# Reference lines for the R–Z axis (r_mm, label)
_TRK_REFS = [(120, "VXD"), (600, "IT"), (1500, "OT")]
_CAL_REFS = [(1500, "ECAL in"), (2100, "HCAL in"), (3400, "HCAL out")]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _load(root_file, tree_name, branches, event_id):
    arr = uproot.open(root_file)[tree_name].arrays(branches, library="np")
    if event_id is not None:
        sel = arr["evt"] == event_id
        if not sel.any():
            available = sorted(np.unique(arr["evt"]).tolist())
            raise ValueError(
                f"Event {event_id} not found in {root_file}. Available: {available}"
            )
        arr = {k: v[sel] for k, v in arr.items()}
    return arr


# ---------------------------------------------------------------------------
# Hit display
# ---------------------------------------------------------------------------

def plot_hit_display(tracker_file, calo_file, out_path, event_id=None):
    trk = (
        _load(tracker_file, "hits", ["evt", "col_id", "x", "y", "z", "r"], event_id)
        if tracker_file else None
    )
    cal = (
        _load(calo_file, "calo_hits", ["evt", "col_id", "x", "y", "z", "r", "energy"], event_id)
        if calo_file else None
    )

    has_trk = trk is not None
    has_cal = cal is not None
    lp      = _layout_params(has_trk, has_cal)
    alpha   = 0.5 if event_id is not None else 0.2

    title_evt = f"event {event_id}" if event_id is not None else "all events"
    n_trk = len(trk["x"]) if has_trk else 0
    n_cal = len(cal["x"]) if has_cal else 0

    # --- figure ---------------------------------------------------------
    fig = plt.figure(figsize=(18, 9))
    fig.patch.set_facecolor("#111111")
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                            width_ratios=[2.2, 1], hspace=0.38, wspace=0.28)
    ax_rz   = fig.add_subplot(gs[:, 0])
    ax_xy   = fig.add_subplot(gs[0, 1])
    ax_zoom = fig.add_subplot(gs[1, 1])
    for ax in (ax_rz, ax_xy, ax_zoom):
        _dark_ax(ax)

    # --- tracker hits ---------------------------------------------------
    if has_trk:
        for cid, (label, color) in enumerate(zip(TRK_LABELS, TRK_COLORS)):
            mask = trk["col_id"] == cid
            kw   = dict(s=5, color=color, alpha=alpha, linewidths=0, rasterized=True)
            ax_rz.scatter(trk["z"][mask], trk["r"][mask],
                          label=f"{label}  ({mask.sum()})", **kw)
            ax_xy.scatter(trk["x"][mask], trk["y"][mask], **kw)
            zm = trk["r"][mask] < lp["zoom_r"]
            ax_zoom.scatter(trk["x"][mask][zm], trk["y"][mask][zm], **kw)

    # --- calorimeter hits -----------------------------------------------
    if has_cal:
        for cid, (label, color) in enumerate(zip(CAL_LABELS, CAL_COLORS)):
            mask   = cal["col_id"] == cid
            energy = cal["energy"][mask]
            sizes  = _energy_to_size(energy)
            kw     = dict(s=sizes, color=color, alpha=alpha, linewidths=0, rasterized=True)
            ax_rz.scatter(cal["z"][mask], cal["r"][mask],
                          label=f"{label}  ({mask.sum()})", **kw)
            ax_xy.scatter(cal["x"][mask], cal["y"][mask], **kw)
            # calo hits are outside the tracker zoom window — draw only in calo-only mode
            if not has_trk:
                zm = cal["r"][mask] < lp["zoom_r"]
                ax_zoom.scatter(cal["x"][mask][zm], cal["y"][mask][zm], **kw)

    # --- R–Z styling ----------------------------------------------------
    ax_rz.set_xlabel("z  [mm]", color="#cccccc", fontsize=10)
    ax_rz.set_ylabel("r  [mm]", color="#cccccc", fontsize=10)
    ax_rz.set_title("R–Z View", color="white", fontsize=12, fontweight="bold", pad=6)
    ax_rz.set_xlim(*lp["rz_xlim"])
    ax_rz.set_ylim(*lp["rz_ylim"])
    ax_rz.grid(True, alpha=0.12, color="#888888", linestyle="--")
    ax_rz.legend(markerscale=3, fontsize=7.5, loc="upper right",
                 facecolor="#1a1a1a", edgecolor="#444444",
                 labelcolor="#cccccc", framealpha=0.85)
    for r_val, lbl in lp["ref_lines"]:
        ax_rz.axhline(r_val, color="#2a2a2a", linestyle=":", linewidth=0.8)
        ax_rz.text(lp["rz_xlim"][1] - 50, r_val + 40, lbl,
                   color="#666666", fontsize=7.0, ha="right")

    # --- X–Y styling ----------------------------------------------------
    _xy_style(ax_xy,   "X–Y  (full detector)",  lp["xy_lim"])
    _xy_style(ax_zoom, lp["zoom_title"],         lp["zoom_r"])

    # --- title ----------------------------------------------------------
    parts = []
    if tracker_file:
        parts.append(f"trk: {os.path.basename(tracker_file)}")
    if calo_file:
        parts.append(f"calo: {os.path.basename(calo_file)}")
    hit_counts = f"{n_trk:,} trk" if has_trk else ""
    if has_cal:
        hit_counts += ("  +  " if hit_counts else "") + f"{n_cal:,} calo"
    fig.suptitle(
        f"MuCol Hit Display  ·  {'  ·  '.join(parts)}  ·  {title_evt}  ·  {hit_counts} hits",
        color="white", fontsize=11, fontweight="bold", y=1.005,
    )

    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"Saved: {out_path}")
    plt.close()


# ---------------------------------------------------------------------------
# Layout parameters (adapt to which systems are present)
# ---------------------------------------------------------------------------

def _layout_params(has_trk, has_cal):
    if has_trk and has_cal:
        return dict(
            rz_xlim=(-4000, 4000), rz_ylim=(-30, 4000),
            xy_lim=3600,
            zoom_r=620,  zoom_title="X–Y  (tracker,  r < 620 mm)",
            ref_lines=_TRK_REFS + _CAL_REFS,
        )
    elif has_trk:
        return dict(
            rz_xlim=(-2500, 2500), rz_ylim=(-30, 1620),
            xy_lim=1600,
            zoom_r=620,  zoom_title="X–Y  (inner,  r < 620 mm)",
            ref_lines=_TRK_REFS,
        )
    else:  # calo only
        return dict(
            rz_xlim=(-4000, 4000), rz_ylim=(-30, 4000),
            xy_lim=3600,
            zoom_r=2200, zoom_title="X–Y  (ECAL,  r < 2200 mm)",
            ref_lines=_CAL_REFS,
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _energy_to_size(energy):
    e_clip = np.clip(energy, 1e-6, None)
    lo, hi = np.log10(1e-6), np.log10(1.0)
    return np.clip(1.5 + 8.0 * (np.log10(e_clip) - lo) / (hi - lo), 1.0, 12.0)


def _dark_ax(ax):
    ax.set_facecolor("#111111")
    ax.tick_params(colors="#cccccc", labelsize=8)
    for spine in ax.spines.values():
        spine.set_edgecolor("#333333")


def _xy_style(ax, title, lim):
    ax.set_xlabel("x  [mm]", color="#cccccc", fontsize=9)
    ax.set_ylabel("y  [mm]", color="#cccccc", fontsize=9)
    ax.set_title(title, color="white", fontsize=9.5, fontweight="bold", pad=4)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.12, color="#888888", linestyle="--")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="MuCol hit display (tracker / calorimeter / combined)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--tracker", nargs="?", const=DEFAULT_TRACKER, default=None, metavar="FILE",
        help="Include tracker hits. FILE defaults to %(const)s",
    )
    parser.add_argument(
        "--calo", nargs="?", const=DEFAULT_CALO, default=None, metavar="FILE",
        help="Include calorimeter hits. FILE defaults to %(const)s",
    )
    parser.add_argument("-e", "--event", metavar="N", type=int, default=1,
                        help="Event ID to display (default: 1)")
    parser.add_argument("-o", "--output", metavar="FILE", default=None,
                        help="Output PNG (default: hit_display_single_evt.png next to this script)")
    args = parser.parse_args()

    if args.tracker is None and args.calo is None:
        parser.error("Specify at least one of --tracker [FILE] and/or --calo [FILE]")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    out = args.output or os.path.join(script_dir, "hit_display_single_evt.png")

    plot_hit_display(args.tracker, args.calo, out, event_id=args.event)


if __name__ == "__main__":
    main()
