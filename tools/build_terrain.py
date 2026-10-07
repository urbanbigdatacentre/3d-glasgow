#!/usr/bin/env python3
"""
build_terrain.py: Cesium quantized-mesh terrain tiles from the UBDC Glasgow 0.5 m DTM.

  1. mosaic   every DTM GeoTIFF (EPSG:27700) is read at --res metres (average resampling) into one
              array, cached as <dtm>/dtm_<res>m.npz so later runs start from here.
  2. surround outside the DTM and in its nodata holes, heights come from Copernicus GLO-30 (EGM2008
              orthometric, within ~1 m of ODN here) feathered over --feather metres, so the DTM edge
              does not show as a cliff. --no-surround uses flat sea level instead.
  3. tiles    for every level 0..--max-level, the tiles covering the area (plus both level-0 roots) are
              sampled on a lon/lat grid, converted to WGS84 ellipsoidal heights with PROJ's OSTN15 +
              OSGM15 grids (EPSG:7405 -> EPSG:4979), meshed with Delatin and written as
              {z}/{x}/{y}.terrain (quantized-mesh-1.0, TMS geodetic tiling). layer.json lists availability.

usage: tools/venv/bin/python tools/build_terrain.py --dtm ../DTM_5x5km --out tiles/terrain
"""
import argparse
import glob
import json
import math
import os
import sys
import time
from multiprocessing import get_context

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.merge import merge as rio_merge
from scipy.ndimage import distance_transform_edt, map_coordinates
from pyproj.transformer import TransformerGroup
from pydelatin import Delatin
from quantized_mesh_encoder import encode

COPERNICUS = "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_{lat}_00_{lon}_00_DEM/Copernicus_DSM_COG_10_{lat}_00_{lon}_00_DEM.tif"

G = {}   # globals shared with forked workers


# ----------------------------------------------------------------------------- data preparation
def build_mosaic(dtm_dir, res):
    cache = os.path.join(dtm_dir, "dtm_%gm.npz" % res)
    if os.path.exists(cache):
        d = np.load(cache)
        print("mosaic: cached %s %s" % (cache, d["arr"].shape))
        return d["arr"], tuple(float(v) for v in d["bounds"])
    files = sorted(glob.glob(os.path.join(dtm_dir, "*.tif")))
    if not files:
        sys.exit("no GeoTIFFs in %s" % dtm_dir)
    bounds = []
    for f in files:
        with rasterio.open(f) as ds:
            bounds.append(ds.bounds)
    x0 = math.floor(min(b.left for b in bounds) / res) * res
    y0 = math.floor(min(b.bottom for b in bounds) / res) * res
    x1 = math.ceil(max(b.right for b in bounds) / res) * res
    y1 = math.ceil(max(b.top for b in bounds) / res) * res
    W, H = int((x1 - x0) / res), int((y1 - y0) / res)
    arr = np.full((H, W), np.nan, dtype=np.float32)
    t0 = time.time()
    for f, b in zip(files, bounds):
        with rasterio.open(f) as ds:
            fac = res / ds.res[0]
            oh, ow = int(round(ds.height / fac)), int(round(ds.width / fac))
            data = ds.read(1, out_shape=(oh, ow), resampling=Resampling.average, masked=True)
            data = np.ma.filled(data.astype(np.float32), np.nan)
            data[(data < -1000) | (data > 5000)] = np.nan          # stray nodata values
        c0 = int(round((b.left - x0) / res))
        r0 = int(round((y1 - b.top) / res))
        arr[r0:r0 + oh, c0:c0 + ow] = data
        print("  read %s -> %dx%d (%.0fs)" % (os.path.basename(f), ow, oh, time.time() - t0), flush=True)
    np.savez_compressed(cache, arr=arr, bounds=np.array([x0, y0, x1, y1]))
    print("mosaic: %dx%d cells at %g m, valid %.1f%%, bounds E%d-%d N%d-%d" %
          (W, H, res, 100 * np.isfinite(arr).mean(), x0, x1, y0, y1))
    return arr, (x0, y0, x1, y1)


def fetch_surround(dtm_dir, lonlat_bounds, margin, res_deg=None, name="copernicus_surround"):
    """Copernicus GLO-30 mosaic (lon/lat grid) covering the area plus margin. Cached.
    res_deg: output cell size in degrees (None = native 30 m); coarse reads come from the COG overviews."""
    cache = os.path.join(dtm_dir, "%s.npz" % name)
    if os.path.exists(cache):
        d = np.load(cache)
        print("%s: cached %s" % (name, str(d["arr"].shape)))
        return d["arr"], tuple(float(v) for v in d["transform"])
    w, s, e, n = lonlat_bounds
    w, s, e, n = w - margin, s - margin, e + margin, n + margin
    urls = []
    for lat in range(math.floor(s), math.ceil(n)):
        for lon in range(math.floor(w), math.ceil(e)):
            urls.append(COPERNICUS.format(lat="N%02d" % lat if lat >= 0 else "S%02d" % -lat,
                                          lon="E%03d" % lon if lon >= 0 else "W%03d" % -lon))
    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    srcs = [rasterio.open(u) for u in urls]
    kw = {"res": (res_deg, res_deg)} if res_deg else {}
    arr, transform = rio_merge(srcs, bounds=(w, s, e, n), resampling=Resampling.average if res_deg else Resampling.bilinear, **kw)
    for ds in srcs:
        ds.close()
    arr = arr[0].astype(np.float32)
    arr[~np.isfinite(arr) | (arr < -1000)] = 0.0
    tr = (transform.c, transform.f, transform.a, transform.e)        # lon0, lat0(top), dlon, dlat(<0)
    np.savez_compressed(cache, arr=arr, transform=np.array(tr))
    print("%s: Copernicus GLO-30 %s from %d tiles" % (name, str(arr.shape), len(urls)))
    return arr, tr


def _sample_grid(grid, lon, lat):
    arr, (lon0, lat0, dlon, dlat) = grid
    cols = (lon - lon0) / dlon - 0.5
    rows = (lat - lat0) / dlat - 0.5
    return map_coordinates(arr, [rows, cols], order=1, mode="nearest")


def sample_surround(lon, lat):
    """Orthometric heights outside the DTM: native-resolution Copernicus near the city, the coarse
    wide-area grid beyond it (so the horizon has real relief instead of a plateau)."""
    near, wide = G["surround"], G.get("surround_wide")
    if near[0] is None:
        return np.zeros_like(lon)
    lon, lat = np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)
    if wide is None or wide[0] is None:
        return _sample_grid(near, lon, lat)
    h = _sample_grid(wide, lon, lat)
    arr, (lon0, lat0, dlon, dlat) = near
    inside = (lon >= lon0) & (lon <= lon0 + dlon * arr.shape[1]) & (lat <= lat0) & (lat >= lat0 + dlat * arr.shape[0])
    if inside.any():
        h[inside] = _sample_grid(near, lon[inside], lat[inside])
    return h


def blend_surround(arr, bounds, res, feather):
    """Fill nodata and feather the DTM edge with the surround (in place on a copy)."""
    x0, y0, x1, y1 = bounds
    H, W = arr.shape
    valid = np.isfinite(arr)
    if G["surround"][0] is None:
        out = np.where(valid, arr, 0.0).astype(np.float32)
        return out
    # surround heights on the mosaic grid (cell centres)
    cols = np.arange(W) * res + x0 + res / 2
    rows = y1 - (np.arange(H) * res + res / 2)
    E, N = np.meshgrid(cols, rows)
    lon, lat = G["to_lonlat"].transform(E.ravel(), N.ravel())
    sur = sample_surround(np.asarray(lon), np.asarray(lat)).reshape(H, W).astype(np.float32)
    dist = distance_transform_edt(valid) * res                          # metres to the nearest hole/edge
    w = np.clip(dist / feather, 0.0, 1.0).astype(np.float32)
    out = np.where(valid, w * np.nan_to_num(arr) + (1 - w) * sur, sur)
    return out.astype(np.float32)


# ----------------------------------------------------------------------------- tiling
def tile_bounds(z, x, y):
    size = 180.0 / (1 << z)
    return (-180 + x * size, -90 + y * size, -180 + (x + 1) * size, -90 + (y + 1) * size)


def tile_range(z, lonlat):
    w, s, e, n = lonlat
    size = 180.0 / (1 << z)
    x0, x1 = int(math.floor((w + 180) / size)), int(math.floor((e + 180) / size))
    y0, y1 = int(math.floor((s + 90) / size)), int(math.floor((n + 90) / size))
    return max(0, x0), max(0, y0), min((2 << z) - 1, x1), min((1 << z) - 1, y1)


def sample_heights(lon, lat):
    """ODN/orthometric heights at lon/lat arrays (DTM mosaic inside its bounds, surround elsewhere)."""
    arr, (x0, y0, x1, y1), res = G["mosaic"], G["bounds"], G["res"]
    E, N = G["to_bng"].transform(lon, lat)
    E, N = np.asarray(E), np.asarray(N)
    cols = (E - x0) / res - 0.5
    rows = (y1 - N) / res - 0.5
    inside = (E >= x0) & (E <= x1) & (N >= y0) & (N <= y1)
    h = sample_surround(lon, lat).astype(np.float64)
    if inside.any():
        h[inside] = map_coordinates(arr, [rows[inside], cols[inside]], order=1, mode="nearest")
    return h, E, N


def build_tile(job):
    z, x, y = job
    w, s, e, n = tile_bounds(z, x, y)
    res, max_level, base_err = G["res"], G["max_level"], G["base_err"]
    lat_c = (s + n) / 2
    width_m = (e - w) * 111320.0 * math.cos(math.radians(lat_c))
    height_m = (n - s) * 111320.0
    nx = int(min(257, max(9, math.ceil(width_m / res) + 1)))
    ny = int(min(257, max(9, math.ceil(height_m / res) + 1)))
    lons = np.linspace(w, e, nx)
    lats = np.linspace(n, s, ny)                                       # row 0 = north
    LON, LAT = np.meshgrid(lons, lats)
    h, E, N = sample_heights(LON.ravel(), LAT.ravel())
    _, _, hell = G["to_ell"].transform(E, N, h)                       # ODN -> WGS84 ellipsoidal
    heights = np.asarray(hell, dtype=np.float32).reshape(ny, nx)
    heights[~np.isfinite(heights)] = 0.0
    err = base_err * (2 ** (max_level - z))
    tin = Delatin(heights, max_error=err)
    v = np.asarray(tin.vertices, dtype=np.float64)
    tri = np.asarray(tin.triangles, dtype=np.uint32).reshape(-1, 3)
    pos = np.empty_like(v)
    pos[:, 0] = w + v[:, 0] / (nx - 1) * (e - w)
    pos[:, 1] = s + v[:, 1] / (ny - 1) * (n - s)      # pydelatin's y counts from the LAST row (south) upwards
    pos[:, 2] = v[:, 2]
    # Cesium back-face culls terrain, so make every triangle counter-clockwise in lon/lat.
    a, b, c = pos[tri[:, 0], :2], pos[tri[:, 1], :2], pos[tri[:, 2], :2]
    cw = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]) < 0
    tri[cw] = tri[cw][:, ::-1]
    path = os.path.join(G["out"], str(z), str(x), "%d.terrain" % y)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        encode(f, pos, tri, bounds=[w, s, e, n], sphere_method="naive")
    return z, len(v), os.path.getsize(path)


def init_worker(g):
    G.update(g)
    G["to_bng"] = TransformerGroup("EPSG:4326", "EPSG:27700", always_xy=True).transformers[0]
    G["to_ell"] = TransformerGroup("EPSG:7405", "EPSG:4979", always_xy=True).transformers[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dtm", required=True, help="folder with the 0.5 m DTM GeoTIFFs")
    ap.add_argument("--out", required=True, help="output folder (layer.json + {z}/{x}/{y}.terrain)")
    ap.add_argument("--res", type=float, default=5.0, help="working resolution in metres")
    ap.add_argument("--max-level", type=int, default=15)
    ap.add_argument("--base-error", type=float, default=0.3, help="mesh error (m) at the finest level; doubles per coarser level")
    ap.add_argument("--feather", type=float, default=400.0, help="blend width (m) between DTM and surround")
    ap.add_argument("--margin", type=float, default=0.35, help="native-resolution surround margin (degrees) around the DTM")
    ap.add_argument("--wide-margin", type=float, default=1.5, help="coarse (300 m) surround margin (degrees) for the horizon, levels <= 10")
    ap.add_argument("--no-surround", action="store_true")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--version", default="1.0.0")
    args = ap.parse_args()

    t0 = time.time()
    to_lonlat = TransformerGroup("EPSG:27700", "EPSG:4326", always_xy=True).transformers[0]
    to_ell = TransformerGroup("EPSG:7405", "EPSG:4979", always_xy=True)
    if not to_ell.best_available:
        sys.exit("PROJ grids missing: run  tools/venv/bin/pyproj sync --file uk_os_OSTN15_NTv2_OSGBtoETRS  and  --file uk_os_OSGM15_GB")
    print("heights: %s" % to_ell.transformers[0].description)

    arr, bounds = build_mosaic(args.dtm, args.res)
    x0, y0, x1, y1 = bounds
    lons, lats = to_lonlat.transform([x0, x1, x0, x1], [y0, y0, y1, y1])
    dtm_ll = (min(lons), min(lats), max(lons), max(lats))
    G["to_lonlat"] = to_lonlat
    G["surround"] = (None, (0, 0, 1, -1))
    G["surround_wide"] = None
    if not args.no_surround:
        try:
            G["surround"] = fetch_surround(args.dtm, dtm_ll, args.margin)
            G["surround_wide"] = fetch_surround(args.dtm, dtm_ll, args.wide_margin, res_deg=0.003, name="copernicus_wide")
        except Exception as e:
            print("surround unavailable (%s); using flat sea level outside the DTM" % e)
    arr = blend_surround(arr, bounds, args.res, args.feather)
    G.update({"mosaic": arr, "bounds": bounds, "res": args.res, "out": args.out,
              "max_level": args.max_level, "base_err": args.base_error})

    # tile list: both roots, then the area (wide margin up to level 12, tight above)
    jobs, available = [], []
    far = (dtm_ll[0] - args.wide_margin, dtm_ll[1] - args.wide_margin, dtm_ll[2] + args.wide_margin, dtm_ll[3] + args.wide_margin)
    near = (dtm_ll[0] - args.margin, dtm_ll[1] - args.margin, dtm_ll[2] + args.margin, dtm_ll[3] + args.margin)
    tight = (dtm_ll[0] - 0.01, dtm_ll[1] - 0.01, dtm_ll[2] + 0.01, dtm_ll[3] + 0.01)
    for z in range(0, args.max_level + 1):
        xa, ya, xb, yb = tile_range(z, far if z <= 10 else near if z <= 12 else tight)
        if z == 0:
            xa, xb = 0, 1
        available.append([{"startX": xa, "startY": ya, "endX": xb, "endY": yb}])
        for x in range(xa, xb + 1):
            for y in range(ya, yb + 1):
                jobs.append((z, x, y))
    print("tiles to build: %d (levels 0-%d)" % (len(jobs), args.max_level), flush=True)

    os.makedirs(args.out, exist_ok=True)
    stats = {}
    ctx = get_context("fork")
    with ctx.Pool(args.workers, initializer=init_worker, initargs=(G,)) as pool:
        for i, (z, nv, nb) in enumerate(pool.imap_unordered(build_tile, jobs, chunksize=8)):
            s = stats.setdefault(z, [0, 0, 0])
            s[0] += 1
            s[1] += nv
            s[2] += nb
            if (i + 1) % 500 == 0:
                print("  %d / %d tiles" % (i + 1, len(jobs)), flush=True)
    total = sum(s[2] for s in stats.values())
    for z in sorted(stats):
        s = stats[z]
        print("  level %2d: %4d tiles, %6.0f vertices/tile, %6.1f KB/tile" % (z, s[0], s[1] / s[0], s[2] / s[0] / 1e3))
    layer = {"tilejson": "2.1.0", "name": "UBDC Glasgow DTM", "version": args.version, "format": "quantized-mesh-1.0",
             "scheme": "tms", "projection": "EPSG:4326", "bounds": [-180, -90, 180, 90],
             "tiles": ["{z}/{x}/{y}.terrain?v={version}"], "minzoom": 0, "maxzoom": args.max_level,
             "attribution": "Terrain: UBDC Glasgow 0.5 m LiDAR DTM (doi:10.5281/zenodo.13273124, CC BY 4.0); surroundings: Copernicus DEM GLO-30",
             "available": available,
             "extensions": [],
             "description": "Built by tools/build_terrain.py: %g m working resolution, Delatin error %g m at level %d, heights ODN -> WGS84 ellipsoid via OSGM15" % (args.res, args.base_error, args.max_level)}
    with open(os.path.join(args.out, "layer.json"), "w") as f:
        json.dump(layer, f, indent=1)
    print("done: %d tiles, %.1f MB in %.0fs -> %s" % (len(jobs), total / 1e6, time.time() - t0, args.out))


if __name__ == "__main__":
    main()
