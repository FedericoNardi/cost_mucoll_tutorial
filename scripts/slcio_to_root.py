#!/usr/bin/env python3
"""
Convert tracker and calorimeter hits from SLCIO to ROOT TTrees.

Writes a single ROOT file containing:
  'hits'      TTree  — tracker (VXD, IT, OT)
                       branches: evt, col_id, layer, side, x, y, z, r, edep, t
  'calo_hits' TTree  — calorimeter (ECal Barrel/Endcap, HCal Barrel/Endcap)
                       branches: evt, col_id, layer, x, y, z, r, energy

Collection types (sim vs. digitised) are auto-detected from the first event
independently for tracker and calorimeter. Use --tracker-only or --calo-only
to write only one of the two trees.

Usage:
    python slcio_to_root.py input.slcio -o detector.root
    python slcio_to_root.py a.slcio b.slcio -o detector.root -m 10
    python slcio_to_root.py input.slcio --list-collections
    python slcio_to_root.py input.slcio --tracker-only -o hits.root
    python slcio_to_root.py input.slcio --calo-only    -o calo_hits.root
"""

import argparse
import sys
import math
import numpy as np

from pyLCIO.drivers.Driver import Driver
from pyLCIO import EVENT, UTIL
from pyLCIO.io.EventLoop import EventLoop
import ROOT as R

# ---------------------------------------------------------------------------
# Known collection tables
# ---------------------------------------------------------------------------

TRK_SIM = [
    ("VXD Barrel", "VertexBarrelCollection",       0),
    ("VXD Endcap", "VertexEndcapCollection",       1),
    ("IT Barrel",  "InnerTrackerBarrelCollection", 2),
    ("IT Endcap",  "InnerTrackerEndcapCollection", 3),
    ("OT Barrel",  "OuterTrackerBarrelCollection", 4),
    ("OT Endcap",  "OuterTrackerEndcapCollection", 5),
]
TRK_DIGI = [
    ("VXD Barrel", "VXDTrackerHits",       0),
    ("VXD Endcap", "VXDEndcapTrackerHits", 1),
    ("IT Barrel",  "ITBarrelHits",         2),
    ("IT Endcap",  "ITEndcapHits",         3),
    ("OT Barrel",  "OTBarrelHits",         4),
    ("OT Endcap",  "OTEndcapHits",         5),
]

CAL_SIM = [
    ("ECAL Barrel", "ECalBarrelCollection", 0),
    ("ECAL Endcap", "ECalEndcapCollection", 1),
    ("HCAL Barrel", "HCalBarrelCollection", 2),
    ("HCAL Endcap", "HCalEndcapCollection", 3),
]
CAL_DIGI = [
    ("ECAL Barrel", "EcalBarrelCollectionDigi", 0),
    ("ECAL Endcap", "EcalEndcapCollectionDigi", 1),
    ("HCAL Barrel", "HcalBarrelCollectionDigi", 2),
    ("HCAL Endcap", "HcalEndcapCollectionDigi", 3),
]


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

class FullDetectorConverterDriver(Driver):
    """Fills a 'hits' TTree (tracker) and/or 'calo_hits' TTree (calorimeter)."""

    def __init__(self, output_path, do_tracker=True, do_calo=True, list_mode=False):
        Driver.__init__(self)
        self.output_path = output_path
        self.do_tracker  = do_tracker
        self.do_calo     = do_calo
        self.list_mode   = list_mode
        self._trk_cols   = None   # resolved on first event
        self._trk_sim    = None
        self._cal_cols   = None
        self._cal_sim    = None

    # ------------------------------------------------------------------
    def startOfData(self):
        self.trk_data = {}
        self.cal_data = {}

        if self.do_tracker:
            self.trk_tree = R.TTree("hits", "Tracker hit positions")
            for name in ("evt", "col_id", "layer", "side"):
                self.trk_data[name] = np.zeros(1, dtype=np.int32)
                self.trk_tree.Branch(name, self.trk_data[name], f"{name}/I")
            for name in ("x", "y", "z", "r", "edep", "t"):
                self.trk_data[name] = np.zeros(1, dtype=np.float32)
                self.trk_tree.Branch(name, self.trk_data[name], f"{name}/F")

        if self.do_calo:
            self.cal_tree = R.TTree("calo_hits", "Calorimeter hit positions and energies")
            for name in ("evt", "col_id", "layer"):
                self.cal_data[name] = np.zeros(1, dtype=np.int32)
                self.cal_tree.Branch(name, self.cal_data[name], f"{name}/I")
            for name in ("x", "y", "z", "r", "energy"):
                self.cal_data[name] = np.zeros(1, dtype=np.float32)
                self.cal_tree.Branch(name, self.cal_data[name], f"{name}/F")

    # ------------------------------------------------------------------
    def _resolve(self, available):
        """On the first event, resolve which collections are present."""
        if self.list_mode:
            print("\nCollections in this file:")
            for name in sorted(available):
                print(f"  {name}")
            sys.exit(0)

        if self.do_tracker:
            for candidates, is_sim in [(TRK_SIM, True), (TRK_DIGI, False)]:
                found = [(l, n, c) for l, n, c in candidates if n in available]
                if found:
                    self._trk_cols = found
                    self._trk_sim  = is_sim
                    kind = "SimTrackerHit" if is_sim else "TrackerHit (digi)"
                    print(f"\nTracker — {kind}:")
                    for label, name, _ in found:
                        print(f"  {label:12s}  {name}")
                    break
            if self._trk_cols is None:
                print("WARNING: no tracker collections found.")
                self._trk_cols = []

        if self.do_calo:
            for candidates, is_sim in [(CAL_SIM, True), (CAL_DIGI, False)]:
                found = [(l, n, c) for l, n, c in candidates if n in available]
                if found:
                    self._cal_cols = found
                    self._cal_sim  = is_sim
                    kind = "SimCalorimeterHit" if is_sim else "CalorimeterHit (digi)"
                    print(f"\nCalorimeter — {kind}:")
                    for label, name, _ in found:
                        print(f"  {label:14s}  {name}")
                    break
            if self._cal_cols is None:
                print("WARNING: no calorimeter collections found.")
                self._cal_cols = []

        print()

    # ------------------------------------------------------------------
    def processEvent(self, event):
        if self._trk_cols is None and self._cal_cols is None:
            self._resolve(set(event.getCollectionNames()))

        evt_num = event.getEventNumber()

        # --- tracker ---
        for _, col_name, col_id in (self._trk_cols or []):
            try:
                col = event.getCollection(col_name)
            except Exception:
                continue
            encoding = col.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder  = UTIL.BitField64(encoding)
            d        = self.trk_data
            for i in range(col.getNumberOfElements()):
                hit = col.getElementAt(i)
                if self._trk_sim:
                    pos  = hit.getPositionVec()
                    x, y, z = pos.X(), pos.Y(), pos.Z()
                    edep = hit.getEDep()
                    t    = hit.getTime()
                else:
                    pos  = hit.getPosition()
                    x, y, z = pos[0], pos[1], pos[2]
                    edep = 0.0
                    t    = 0.0
                cell_id = (int(hit.getCellID0()) & 0xFFFFFFFF) | (int(hit.getCellID1()) << 32)
                decoder.setValue(cell_id)
                d["evt"][0]    = evt_num
                d["col_id"][0] = col_id
                d["layer"][0]  = int(decoder["layer"].value()) if "layer" in encoding else -1
                d["side"][0]   = int(decoder["side"].value())  if "side"  in encoding else 0
                d["x"][0]      = x;  d["y"][0] = y;  d["z"][0] = z
                d["r"][0]      = math.sqrt(x*x + y*y)
                d["edep"][0]   = edep;  d["t"][0] = t
                self.trk_tree.Fill()

        # --- calorimeter ---
        for _, col_name, col_id in (self._cal_cols or []):
            try:
                col = event.getCollection(col_name)
            except Exception:
                continue
            encoding = col.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder  = UTIL.BitField64(encoding)
            d        = self.cal_data
            for i in range(col.getNumberOfElements()):
                hit = col.getElementAt(i)
                pos  = hit.getPosition()
                x, y, z = pos[0], pos[1], pos[2]
                cell_id = (int(hit.getCellID0()) & 0xFFFFFFFF) | (int(hit.getCellID1()) << 32)
                decoder.setValue(cell_id)
                d["evt"][0]    = evt_num
                d["col_id"][0] = col_id
                d["layer"][0]  = int(decoder["layer"].value()) if "layer" in encoding else -1
                d["x"][0]      = x;  d["y"][0] = y;  d["z"][0] = z
                d["r"][0]      = math.sqrt(x*x + y*y)
                d["energy"][0] = hit.getEnergy()
                self.cal_tree.Fill()

        n_trk = self.trk_tree.GetEntries() if self.do_tracker else 0
        n_cal = self.cal_tree.GetEntries() if self.do_calo    else 0
        print(f"  Event {evt_num:4d}  |  trk: {n_trk}  calo: {n_cal}")

    # ------------------------------------------------------------------
    def endOfData(self):
        out = R.TFile(self.output_path, "RECREATE")
        if self.do_tracker:
            self.trk_tree.Write()
        if self.do_calo:
            self.cal_tree.Write()
        out.Close()
        n_trk = self.trk_tree.GetEntries() if self.do_tracker else 0
        n_cal = self.cal_tree.GetEntries() if self.do_calo    else 0
        print(f"\nSaved → {self.output_path}  ({n_trk} tracker + {n_cal} calo hits)")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert SLCIO hits (tracker + calorimeter) to ROOT TTrees",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", metavar="input.slcio", nargs="+",
                        help="Input SLCIO file(s)")
    parser.add_argument("-o", "--output", metavar="OUT.root", default="detector.root",
                        help="Output ROOT file path")
    parser.add_argument("-m", "--max-events", metavar="N", type=int, default=-1,
                        help="Maximum number of events to process (-1 = all)")
    parser.add_argument("-s", "--skip-events", metavar="N", type=int, default=0,
                        help="Number of events to skip at the start")
    parser.add_argument("--list-collections", action="store_true",
                        help="Print all collections in the first event and exit")
    parser.add_argument("--tracker-only", action="store_true",
                        help="Write only the tracker 'hits' tree")
    parser.add_argument("--calo-only", action="store_true",
                        help="Write only the calorimeter 'calo_hits' tree")
    opts = parser.parse_args()

    if opts.tracker_only and opts.calo_only:
        parser.error("--tracker-only and --calo-only are mutually exclusive")

    loop    = EventLoop()
    n_total = 0
    for path in opts.input:
        loop.addFile(path)
        n_total += loop.reader.getNumberOfEvents()
        print(f"Added: {path}")

    n_process = n_total if opts.max_events < 0 else min(opts.max_events, n_total)
    print(f"Events in file(s): {n_total}  |  Will process: {n_process}")

    driver = FullDetectorConverterDriver(
        output_path=opts.output,
        do_tracker=not opts.calo_only,
        do_calo=not opts.tracker_only,
        list_mode=opts.list_collections,
    )
    loop.add(driver)

    if opts.skip_events:
        loop.skipEvents(opts.skip_events)

    loop.loop(n_process)
    loop.printStatistics()
    print("Done.")


if __name__ == "__main__":
    main()
