#!/usr/bin/env python3
"""
assemble.py: merge the per-tile outputs under tiles/<OSTILE>/ into
  tiles/tileset.json  - master tileset referencing each OS tile's tileset.json (external tilesets)
  tiles/summary.json  - city-wide statistics for the dashboard, computed from every tile's
                        buildings.csv (histogram, medians, tallest building, per-tile table)

usage: python3 tools/assemble.py [tiles]
"""
import csv
import glob
import json
import math
import os
import sys
import time

import numpy as np
from pyproj.transformer import TransformerGroup

TILES_DIR = sys.argv[1] if len(sys.argv) > 1 else "tiles"
# histogram bin edges (m); the dashboard's height classes (<6, 6-10, 10-15, 15-25, >25) are all edges
HIST_BINS = [0, 3, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 1e9]
NODATA = -9999.0

tr = TransformerGroup("EPSG:27700", "EPSG:4326", always_xy=True).transformers[0]


def read_csv(path):
    rows = {"lod": [], "height": [], "height_max": [], "footprint_area": [], "centre_e": [], "centre_n": [], "uid": []}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows["lod"].append(int(r["lod"]))
            rows["height"].append(float(r["height"]))
            rows["height_max"].append(float(r["height_max"]))
            rows["footprint_area"].append(float(r["footprint_area"]))
            rows["centre_e"].append(float(r["centre_e"]))
            rows["centre_n"].append(float(r["centre_n"]))
            rows["uid"].append(r["uid"])
    out = {k: (np.array(v) if k != "uid" else v) for k, v in rows.items()}
    return out


children, per_tile, regions = [], [], []
all_h, all_area = [], []
tot = {"buildings": 0, "buildings_lod2": 0, "triangles_lod1": 0, "triangles_lod2": 0,
       "bytes_lod1": 0, "bytes_lod2": 0, "glbs": 0}
tallest = {"height_max_m": -1.0}
for idx_path in sorted(glob.glob(os.path.join(TILES_DIR, "*", "index.json"))):
    d = json.load(open(idx_path))
    tile = d["tile"]
    tdir = os.path.dirname(idx_path)
    ts = json.load(open(os.path.join(tdir, "tileset.json")))
    region = ts["root"]["boundingVolume"]["region"]
    regions.append(region)
    children.append({"boundingVolume": {"region": region}, "geometricError": ts["root"]["geometricError"],
                     "content": {"uri": "%s/tileset.json" % tile},
                     "extras": {"tile": tile, "lod2": d["lod2"] is not None}})
    b = read_csv(os.path.join(tdir, "buildings.csv"))
    l1 = b["lod"] == 1
    l2 = b["lod"] == 2
    h = b["height"][l1]
    h = h[h != NODATA]
    a = b["footprint_area"][l1]
    a = a[a != NODATA]
    all_h.append(h)
    all_area.append(a)
    s1, s2 = d["lod1"], d["lod2"]
    tot["buildings"] += int(l1.sum())
    tot["triangles_lod1"] += s1["triangles"]
    tot["bytes_lod1"] += s1["bytes"]
    tot["glbs"] += s1["cells"]
    row = {"tile": tile, "buildings": int(l1.sum()), "lod2": s2 is not None,
           "median_height_m": round(float(np.median(h)), 1) if len(h) else None,
           "footprint_km2": round(float(a.sum()) / 1e6, 3),
           "mb": round((s1["bytes"] + (s2["bytes"] if s2 else 0)) / 1e6, 1),
           "centre_lon": round(math.degrees((region[0] + region[2]) / 2), 5),
           "centre_lat": round(math.degrees((region[1] + region[3]) / 2), 5),
           "region_deg": [round(math.degrees(v), 5) for v in region[:4]]}
    if s2:
        tot["buildings_lod2"] += int(l2.sum())
        tot["triangles_lod2"] += s2["triangles"]
        tot["bytes_lod2"] += s2["bytes"]
        tot["glbs"] += s2["cells"]
        hm = b["height_max"].copy()
        hm[~l2] = -1
        i = int(np.argmax(hm))
        lon, lat = tr.transform(b["centre_e"][i], b["centre_n"][i])
        row["tallest_m"] = round(float(hm[i]), 1)
        row["tallest_uid"] = b["uid"][i]
        if hm[i] > tallest["height_max_m"]:
            tallest = {"height_max_m": round(float(hm[i]), 1), "uid": b["uid"][i], "tile": tile,
                       "lon": round(float(lon), 6), "lat": round(float(lat), 6)}
    else:
        hh = b["height"].copy()
        hh[~l1] = -1
        i = int(np.argmax(hh))
        row["tallest_m"] = round(float(hh[i]), 1)
        row["tallest_uid"] = b["uid"][i]
    per_tile.append(row)

if not children:
    sys.exit("no tiles found under %s" % TILES_DIR)

H = np.concatenate(all_h)
A = np.concatenate(all_area)
hist = np.histogram(H, bins=HIST_BINS)[0].tolist()
root_region = [min(r[0] for r in regions), min(r[1] for r in regions), max(r[2] for r in regions),
               max(r[3] for r in regions), min(r[4] for r in regions), max(r[5] for r in regions)]
master = {"asset": {"version": "1.1", "generator": "ubdc assemble.py",
                    "extras": {"source": "UBDC Glasgow 3D building model, doi:10.5281/zenodo.15000747"}},
          "geometricError": 1000000.0,
          "root": {"boundingVolume": {"region": root_region}, "geometricError": 1000000.0, "refine": "ADD",
                   "children": children}}
with open(os.path.join(TILES_DIR, "tileset.json"), "w") as f:
    json.dump(master, f, separators=(",", ":"))

per_tile.sort(key=lambda r: -r["buildings"])
summary = {
    "built": time.strftime("%Y-%m-%d"),
    "tiles": len(per_tile), "tiles_with_lod2": sum(1 for r in per_tile if r["lod2"]),
    "totals": {"buildings": tot["buildings"], "buildings_lod2": tot["buildings_lod2"],
               "triangles_lod1": tot["triangles_lod1"], "triangles_lod2": tot["triangles_lod2"],
               "mb_lod1": round(tot["bytes_lod1"] / 1e6, 1), "mb_lod2": round(tot["bytes_lod2"] / 1e6, 1),
               "mb_total": round((tot["bytes_lod1"] + tot["bytes_lod2"]) / 1e6, 1), "glbs": tot["glbs"],
               "footprint_km2": round(float(A.sum()) / 1e6, 2)},
    "height": {"bins": HIST_BINS[:-1], "hist": hist, "n": int(len(H)),
               "median_m": round(float(np.median(H)), 1), "mean_m": round(float(H.mean()), 1),
               "p90_m": round(float(np.percentile(H, 90)), 1)},
    "tallest": tallest,
    "root_region_deg": [round(math.degrees(v), 5) for v in root_region[:4]],
    "focus": {"lon": round(sum(r["buildings"] * r["centre_lon"] for r in per_tile) / max(1, tot["buildings"]), 5),
              "lat": round(sum(r["buildings"] * r["centre_lat"] for r in per_tile) / max(1, tot["buildings"]), 5)},
    "per_tile": per_tile,
}
with open(os.path.join(TILES_DIR, "summary.json"), "w") as f:
    json.dump(summary, f, indent=1)
print("master tileset: %d tiles (%d with LoD2), %d GLBs, %.1f MB; buildings %d (LoD2 %d); median height %.1f m; tallest %s %.1f m" %
      (len(per_tile), summary["tiles_with_lod2"], tot["glbs"], summary["totals"]["mb_total"], tot["buildings"],
       tot["buildings_lod2"], summary["height"]["median_m"], tallest.get("uid"), tallest["height_max_m"]))
