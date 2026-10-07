#!/usr/bin/env python3
"""
build_tileset.py: one OS 5 km tile -> a 3D Tiles 1.1 tileset with two levels of detail.

  root (refine REPLACE)
   |- LoD1 cells  (625 m, coarse, geometricError ~50 m)   lod1/<ix>_<iy>.glb
       |- LoD2 cells (312.5 m, fine, geometricError 0)     lod2/<ix>_<iy>.glb

Far away Cesium shows the LoD1 prisms; within a couple of km it swaps in the LoD2 roofs.
Tiles without LoD2 data (the 19 outer tiles) are LoD1-only.

Each GLB is glTF 2.0 with
  - positions quantised to int16 (KHR_mesh_quantization, node scale dequantises)
  - vertex + index streams meshopt-compressed (EXT_meshopt_compression)
  - per-vertex feature id -> building (EXT_mesh_features)
  - a building property table (EXT_structural_metadata)
  - no normals: the viewer shades from screen-space derivatives

Input: folders of CityJSON 1.0.x files, one Building per file, EPSG:27700, ODN heights.
Dependencies: numpy, pyproj, tools/meshopt.py (ctypes binding to vendor/libmeshoptimizer.dylib).

usage:
  python3 build_tileset.py --tile NS56SE --lod1 <folder> [--lod2 <folder>] --out tiles/NS56SE
"""
import argparse
import csv
import glob
import json
import math
import os
import struct
import sys
import time
from multiprocessing import Pool

import numpy as np
from pyproj.transformer import TransformerGroup

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor", "py"))
import meshopt  # noqa: E402
import mapbox_earcut  # noqa: E402
from shapely.geometry import Polygon, MultiPolygon  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
NODATA = -9999.0

# unified property set: (name, CityJSON key in LoD1 files, key in LoD2 files)
PROPS = [
    ("height", "Building_Height", "MEAN"),          # identical definitions, verified
    ("height_max", None, "MAX"),
    ("height_min", None, "MIN"),
    ("height_pct90", None, "PCT90"),
    ("ground_z", "Ground_Z", "Ground_Z"),
    ("footprint_area", "Footprint_Area", "Footprint_Area"),
    ("footprint_length", "Footprint_Length", "Footprint_Length"),
]
HEIGHT_BINS = [0, 3, 6, 9, 12, 15, 20, 30, 50, 1e9]

GLB_MAGIC = 0x46546C67
CHUNK_JSON = 0x4E4F534A
CHUNK_BIN = 0x004E4942


# ----------------------------------------------------------------------------
# OS National Grid tile names -> origin of the 5 km tile
# ----------------------------------------------------------------------------
def os_tile_origin(name):
    name = name.upper()
    def idx(ch):
        i = ord(ch) - ord("A")
        return i - 1 if i > 8 else i          # the letter I is skipped
    l1, l2 = idx(name[0]), idx(name[1])
    e = (((l1 - 2) % 5) * 5 + (l2 % 5)) * 100000
    n = ((19 - (l1 // 5) * 5) - (l2 // 5)) * 100000
    e += int(name[2]) * 10000
    n += int(name[3]) * 10000
    quad = name[4:6]
    if quad in ("NE", "SE"):
        e += 5000
    if quad in ("NE", "NW"):
        n += 5000
    return e, n


# ----------------------------------------------------------------------------
# LoD1 prism rebuild: some tiles (e.g. NS56NE) have raster-traced outlines with a vertex
# every 0.5 m, so a plain prism carries ~1000 triangles. Recover the roof outline,
# simplify it with Douglas-Peucker and re-extrude.
# ----------------------------------------------------------------------------
def roof_rings(V, T, zmax):
    """Closed vertex-index loops of the roof outline (edges used once among roof triangles)."""
    z = V[:, 2]
    top = T[np.all(np.abs(z[T] - zmax) < 2e-3, axis=1)]
    if len(top) == 0:
        return None
    e = np.concatenate([top[:, [0, 1]], top[:, [1, 2]], top[:, [2, 0]]])
    e = np.sort(e, axis=1)
    uniq, cnt = np.unique(e, axis=0, return_counts=True)
    bedges = uniq[cnt == 1]
    adj = {}
    for a, b in bedges.tolist():
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    if any(len(v) != 2 for v in adj.values()):
        return None                      # not a clean outline
    rings, seen = [], set()
    for start in adj:
        if start in seen:
            continue
        ring, prev, cur = [start], None, start
        seen.add(start)
        while True:
            nb = adj[cur]
            nxt = nb[0] if nb[0] != prev else nb[1]
            if nxt == start:
                break
            if nxt in seen:
                return None
            ring.append(nxt)
            seen.add(nxt)
            prev, cur = cur, nxt
        if len(ring) >= 3:
            rings.append(ring)
    return rings or None


def rebuild_prism(V, T, tol):
    """Returns (V2, T2) for a simplified prism, or None to keep the original mesh."""
    z = V[:, 2]
    zmin, zmax = float(z.min()), float(z.max())
    if zmax - zmin < 0.05:
        return None
    distinct = np.unique(np.round(z, 2))
    if len(distinct) != 2:
        return None                      # not a flat-roofed prism
    rings = roof_rings(V, T, zmax)
    footprints = []                      # valid shapely polygons (with holes) of the roof
    if rings is not None:
        polys = [Polygon(V[r, :2]) for r in rings]
        polys = [p for p in polys if p.area > 1e-6]
        polys.sort(key=lambda p: p.area, reverse=True)
        outers = []                      # [shell, [holes...]]
        for p in polys:
            rp = p.representative_point()
            for o in outers:
                if o[0].contains(rp):
                    o[1].append(p)
                    break
            else:
                outers.append([p, []])
        for shell, holes in outers:
            poly = Polygon(shell.exterior.coords, [h.exterior.coords for h in holes])
            if not poly.is_valid:
                poly = poly.buffer(0)
            footprints.append(poly)
    else:
        # outline touches itself (raster-traced "pinch" vertices): union the roof triangles instead
        top = T[np.all(np.abs(z[T] - zmax) < 2e-3, axis=1)]
        u = unary_union([Polygon(V[t, :2]) for t in top]).buffer(0)
        footprints = list(u.geoms) if isinstance(u, MultiPolygon) else [u]
    parts = []
    for poly in footprints:
        if poly.is_empty:
            continue
        simp = poly.simplify(tol, preserve_topology=True)
        geoms = simp.geoms if isinstance(simp, MultiPolygon) else [simp]
        for g in geoms:
            if isinstance(g, Polygon) and not g.is_empty and g.area > 0.5:
                parts.append(g)
    if not parts:
        return None
    # guard: the simplified footprint must keep the source roof area (within 5 %)
    top = T[np.all(np.abs(z[T] - zmax) < 2e-3, axis=1)]
    P = V[top][:, :, :2]
    src_area = 0.5 * np.abs((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) -
                            (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1])).sum()
    new_area = sum(g.area for g in parts)
    if src_area <= 0 or abs(new_area - src_area) > 0.05 * src_area:
        return None
    verts2d, ring_ends, wall_edges = [], [], []
    for g in parts:
        base = len(verts2d)
        for ring in [g.exterior] + list(g.interiors):
            coords = np.asarray(ring.coords)[:-1]
            if len(coords) < 3:
                continue
            start = len(verts2d)
            verts2d.extend(coords.tolist())
            ring_ends.append(len(verts2d))
            n = len(coords)
            wall_edges.extend([(start + i, start + (i + 1) % n) for i in range(n)])
        # earcut per part so holes attach to the right shell
        pv = np.asarray(verts2d[base:], dtype=np.float64)
        pe = np.asarray([r - base for r in ring_ends if r > base], dtype=np.uint32)
        tri = mapbox_earcut.triangulate_float64(pv, pe).reshape(-1, 3) + base
        parts_tris = tri if base == 0 and len(parts) == 1 else tri
        if base == 0:
            roof = [parts_tris]
        else:
            roof.append(parts_tris)
    roof = np.concatenate(roof) if len(roof) > 1 else roof[0]
    n = len(verts2d)
    xy = np.asarray(verts2d, dtype=np.float64)
    V2 = np.concatenate([np.column_stack([xy, np.full(n, zmax)]), np.column_stack([xy, np.full(n, zmin)])])
    walls = []
    for i, j in wall_edges:
        walls.append([i, j, j + n])
        walls.append([i, j + n, i + n])
    T2 = np.concatenate([roof, roof[:, ::-1] + n, np.asarray(walls, dtype=np.int64)]).astype(np.int32)
    return V2, T2


# ----------------------------------------------------------------------------
# parsing (worker processes)
# ----------------------------------------------------------------------------
def parse_file(job):
    path, lod, outline_tol = job
    with open(path, "rb") as f:
        d = json.load(f)
    V = np.asarray(d["vertices"], dtype=np.float64)
    out = []
    for oid, obj in d["CityObjects"].items():
        tris = []
        for g in obj.get("geometry", []):
            gtype = g["type"]
            if gtype in ("MultiSurface", "CompositeSurface"):
                polys = g["boundaries"]
            elif gtype == "Solid":
                polys = [p for shell in g["boundaries"] for p in shell]
            else:
                raise ValueError("unsupported geometry type %s in %s" % (gtype, path))
            for poly in polys:
                ring = poly[0]
                n = len(ring)
                if n == 3:
                    tris.append(ring)
                elif n > 3:
                    for i in range(1, n - 1):
                        tris.append([ring[0], ring[i], ring[i + 1]])
        if not tris:
            continue
        T = np.asarray(tris, dtype=np.int64)
        used = np.unique(T)
        remap = np.full(len(V), -1, dtype=np.int64)
        remap[used] = np.arange(len(used))
        a = obj.get("attributes", {}) or {}
        props = {}
        for name, k1, k2 in PROPS:
            key = k1 if lod == 1 else k2
            v = a.get(key) if key else None
            props[name] = float(v) if v is not None else NODATA
        uid = str(a.get("U_ID") or oid)
        Vb, Tb = V[used], remap[T].astype(np.int32)
        src_tris = len(Tb)
        rebuilt = False
        if lod == 1 and outline_tol > 0:
            try:
                res = rebuild_prism(Vb, Tb, outline_tol)
            except Exception:
                res = None
            if res is not None:
                Vb, Tb = res
                rebuilt = True
        props["_src_triangles"] = src_tris
        props["_rebuilt"] = rebuilt
        out.append((uid, Vb, Tb, props))
    return out


# ----------------------------------------------------------------------------
# geodesy
# ----------------------------------------------------------------------------
def geodetic_to_ecef(lon_deg, lat_deg, h):
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    sl, cl = np.sin(lat), np.cos(lat)
    so, co = np.sin(lon), np.cos(lon)
    N = WGS84_A / np.sqrt(1.0 - WGS84_E2 * sl * sl)
    return np.stack([(N + h) * cl * co, (N + h) * cl * so, (N * (1.0 - WGS84_E2) + h) * sl], axis=1)


def enu_axes(lon_deg, lat_deg):
    lon, lat = math.radians(lon_deg), math.radians(lat_deg)
    east = np.array([-math.sin(lon), math.cos(lon), 0.0])
    north = np.array([-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat)])
    up = np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)])
    return np.stack([east, north, up], axis=1)


# ----------------------------------------------------------------------------
# GLB writer
# ----------------------------------------------------------------------------
def _pad(b, align, fill=b"\x00"):
    r = len(b) % align
    return b if r == 0 else b + fill * (align - r)


def _align(n, a):
    return (n + a - 1) // a * a


def build_schema():
    props = {"uid": {"name": "Building ID", "type": "STRING"},
             "lod": {"name": "Level of detail", "type": "SCALAR", "componentType": "UINT8"},
             "n_triangles": {"name": "Triangles", "type": "SCALAR", "componentType": "UINT32"}}
    for name, k1, k2 in PROPS:
        props[name] = {"name": name, "type": "SCALAR", "componentType": "FLOAT32", "noData": NODATA}
    return {"id": "ubdc_glasgow_3d", "name": "UBDC Glasgow 3D building model",
            "classes": {"building": {"name": "Building", "properties": props}}}


def write_glb(path, local, tris, fids, table, schema, name, compress=True):
    """
    local : (n,3) float positions in the glTF y-up local frame (metres)
    tris  : (m,3) int triangle indices
    fids  : (n,) int feature id per vertex (index into the property table)
    table : dict property name -> np.ndarray (numeric) | list[str]
    Returns the GLB size in bytes.
    """
    local = np.asarray(local, dtype=np.float64)
    n = len(local)
    n_feat = len(next(iter(table.values())))

    # --- quantise positions to int16; the node scale restores metres ------------
    amax = np.abs(local).max(axis=0)
    scale = np.maximum(amax / 32767.0, 1e-4)
    q = np.ascontiguousarray(np.clip(np.round(local / scale), -32767, 32767).astype(np.int16))
    rec = np.zeros((n, 12), dtype=np.uint8)            # 12-byte vertex: int16 xyz, pad, uint16 fid, pad
    rec[:, 0:6] = q.view(np.uint8)
    rec[:, 8:10] = np.ascontiguousarray(np.asarray(fids, dtype=np.uint16)).view(np.uint8).reshape(n, 2)
    idx = np.ascontiguousarray(tris, dtype=np.uint32).ravel()
    if compress:
        rec, idx = meshopt.optimize(rec, idx)
        n = len(rec)
    qpos = np.ascontiguousarray(rec[:, 0:6]).view(np.int16).reshape(n, 3)
    use32 = n > 65535
    isz = 4 if use32 else 2

    parts = []
    views = []
    offset = 0

    def add_raw(raw):
        nonlocal offset
        off = offset
        padded = _pad(raw, 8)
        parts.append(padded)
        offset += len(padded)
        return off

    if compress:
        ev = meshopt.encode_vertex_buffer(rec)
        ei = meshopt.encode_index_buffer(idx, n)
        off_v = add_raw(ev)
        off_i = add_raw(ei)
        fb_v_len = n * 12
        fb_i_off = _align(fb_v_len, 8)
        fb_i_len = int(idx.size) * isz
        fallback_len = fb_i_off + fb_i_len
        views.append({"buffer": 1, "byteOffset": 0, "byteLength": fb_v_len, "byteStride": 12, "target": 34962,
                      "extensions": {"EXT_meshopt_compression": {
                          "buffer": 0, "byteOffset": off_v, "byteLength": len(ev),
                          "byteStride": 12, "count": n, "mode": "ATTRIBUTES"}}})
        views.append({"buffer": 1, "byteOffset": fb_i_off, "byteLength": fb_i_len, "target": 34963,
                      "extensions": {"EXT_meshopt_compression": {
                          "buffer": 0, "byteOffset": off_i, "byteLength": len(ei),
                          "byteStride": isz, "count": int(idx.size), "mode": "TRIANGLES"}}})
    else:
        off_v = add_raw(rec.tobytes())
        views.append({"buffer": 0, "byteOffset": off_v, "byteLength": n * 12, "byteStride": 12, "target": 34962})
        ib = idx.astype(np.uint32 if use32 else np.uint16).tobytes()
        off_i = add_raw(ib)
        views.append({"buffer": 0, "byteOffset": off_i, "byteLength": len(ib), "target": 34963})

    accessors = [
        {"bufferView": 0, "byteOffset": 0, "componentType": 5122, "count": n, "type": "VEC3",
         "min": [int(v) for v in qpos.min(axis=0)], "max": [int(v) for v in qpos.max(axis=0)]},
        {"bufferView": 0, "byteOffset": 8, "componentType": 5123, "count": n, "type": "SCALAR"},
        {"bufferView": 1, "componentType": 5125 if use32 else 5123, "count": int(idx.size), "type": "SCALAR"},
    ]

    # --- property table (uncompressed, in buffer 0, 8-byte aligned, exact lengths) --
    pt_props = {}
    for pname, col in table.items():
        if isinstance(col, list):
            encoded = [s.encode("utf-8") for s in col]
            offs = np.zeros(len(encoded) + 1, dtype=np.uint32)
            offs[1:] = np.cumsum([len(e) for e in encoded])
            raw = b"".join(encoded)
            views.append({"buffer": 0, "byteOffset": add_raw(raw), "byteLength": len(raw)})
            bv_vals = len(views) - 1
            raw = offs.tobytes()
            views.append({"buffer": 0, "byteOffset": add_raw(raw), "byteLength": len(raw)})
            pt_props[pname] = {"values": bv_vals, "stringOffsets": len(views) - 1, "stringOffsetType": "UINT32"}
        else:
            raw = np.ascontiguousarray(col).tobytes()
            views.append({"buffer": 0, "byteOffset": add_raw(raw), "byteLength": len(raw)})
            pt_props[pname] = {"values": len(views) - 1}

    binary = b"".join(parts)
    buffers = [{"byteLength": len(binary)}]
    ext_used = ["EXT_mesh_features", "EXT_structural_metadata", "KHR_mesh_quantization"]
    ext_req = ["KHR_mesh_quantization"]
    if compress:
        buffers.append({"byteLength": fallback_len, "extensions": {"EXT_meshopt_compression": {"fallback": True}}})
        ext_used.append("EXT_meshopt_compression")
        ext_req.append("EXT_meshopt_compression")

    gltf = {
        "asset": {"version": "2.0", "generator": "ubdc build_tileset.py"},
        "extensionsUsed": ext_used,
        "extensionsRequired": ext_req,
        "scene": 0,
        "scenes": [{"nodes": [0]}],
        "nodes": [{"mesh": 0, "name": name, "scale": [float(s) for s in scale]}],
        "meshes": [{"name": name, "primitives": [{
            "attributes": {"POSITION": 0, "_FEATURE_ID_0": 1},
            "indices": 2, "mode": 4, "material": 0,
            "extensions": {"EXT_mesh_features": {"featureIds": [
                {"featureCount": n_feat, "attribute": 0, "propertyTable": 0}]}},
        }]}],
        "materials": [{"name": "building", "doubleSided": True,
                       "pbrMetallicRoughness": {"baseColorFactor": [1.0, 1.0, 1.0, 1.0],
                                                "metallicFactor": 0.0, "roughnessFactor": 1.0}}],
        "accessors": accessors,
        "bufferViews": views,
        "buffers": buffers,
        "extensions": {"EXT_structural_metadata": {
            "schema": schema,
            "propertyTables": [{"name": "buildings", "class": "building", "count": n_feat, "properties": pt_props}]}},
    }
    json_bytes = _pad(json.dumps(gltf, separators=(",", ":")).encode("utf-8"), 4, b" ")
    binary = _pad(binary, 4)
    total = 12 + 8 + len(json_bytes) + 8 + len(binary)
    with open(path, "wb") as f:
        f.write(struct.pack("<III", GLB_MAGIC, 2, total))
        f.write(struct.pack("<II", len(json_bytes), CHUNK_JSON))
        f.write(json_bytes)
        f.write(struct.pack("<II", len(binary), CHUNK_BIN))
        f.write(binary)
    return total


# ----------------------------------------------------------------------------
# one level of detail -> cells -> GLBs
# ----------------------------------------------------------------------------
def load_lod(folder, lod, workers, outline_tol=0.0):
    files = sorted(glob.glob(os.path.join(folder, "*.json")))
    if not files:
        sys.exit("no .json files in %s" % folder)
    t0 = time.time()
    with Pool(workers) as pool:
        buildings = [b for chunk in pool.imap_unordered(parse_file, [(f, lod, outline_tol) for f in files], chunksize=64)
                     for b in chunk]
    buildings.sort(key=lambda b: b[0])
    rebuilt = sum(1 for b in buildings if b[3]["_rebuilt"])
    src = sum(b[3]["_src_triangles"] for b in buildings)
    now = sum(len(b[2]) for b in buildings)
    print("  lod%d: parsed %d buildings from %d files in %.1fs; prisms rebuilt %d; triangles %d -> %d" %
          (lod, len(buildings), len(files), time.time() - t0, rebuilt, src, now))
    return buildings


def build_level(buildings, lod, cell, x0, y0, out_dir, transformer, flatten, compress, schema, tile_name,
                simplify_error=0.0):
    """Returns (cells dict (ix,iy) -> tile dict, stats dict, csv rows)."""
    nb = len(buildings)
    counts = np.array([len(b[1]) for b in buildings], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    V = np.concatenate([b[1] for b in buildings])
    owner = np.repeat(np.arange(nb), counts)
    zmin = np.full(nb, np.inf)
    np.minimum.at(zmin, owner, V[:, 2])
    if flatten:
        h = V[:, 2] - zmin[owner]                        # every building sits on the ellipsoid (no terrain)
        lon, lat = transformer.transform(V[:, 0], V[:, 1])
    else:                                                # true ODN heights -> WGS84 ellipsoidal (for terrain)
        lon, lat, h = transformer.transform(V[:, 0], V[:, 1], V[:, 2])
    lon, lat, h = np.asarray(lon), np.asarray(lat), np.asarray(h)
    P = geodetic_to_ecef(lon, lat, h)

    bmin = np.full((nb, 2), np.inf)
    bmax = np.full((nb, 2), -np.inf)
    np.minimum.at(bmin, owner, V[:, :2])
    np.maximum.at(bmax, owner, V[:, :2])
    centre = (bmin + bmax) / 2.0
    ix = np.floor((centre[:, 0] - x0) / cell).astype(int)
    iy = np.floor((centre[:, 1] - y0) / cell).astype(int)
    cells = {}
    for b in range(nb):
        cells.setdefault((int(ix[b]), int(iy[b])), []).append(b)

    os.makedirs(out_dir, exist_ok=True)
    tiles = {}
    total_bytes = 0
    total_tris = 0
    csv_rows = []
    t0 = time.time()
    for key, blist in sorted(cells.items()):
        cell_name = "%d_%d" % key
        sel = np.concatenate([np.arange(starts[b], starts[b] + counts[b]) for b in blist])
        lon_c, lat_c, h_c = lon[sel], lat[sel], h[sel]
        clon = float((lon_c.min() + lon_c.max()) / 2.0)
        clat = float((lat_c.min() + lat_c.max()) / 2.0)
        R = enu_axes(clon, clat)
        C = geodetic_to_ecef(np.array([clon]), np.array([clat]), np.array([0.0]))[0]
        enu = (P[sel] - C) @ R
        local = np.stack([enu[:, 0], enu[:, 2], -enu[:, 1]], axis=1)   # glTF y-up

        tri_parts, fid_parts = [], []
        off = 0
        for k, b in enumerate(blist):
            tri_parts.append(buildings[b][2] + off)
            fid_parts.append(np.full(counts[b], k, dtype=np.uint16))
            off += counts[b]
        tris = np.concatenate(tri_parts)
        fids = np.concatenate(fid_parts)
        if simplify_error > 0:
            # error-bounded decimation (metres); buildings are separate components so no cracks
            idx, _ = meshopt.simplify(local.astype(np.float32), tris.ravel(), simplify_error)
            used = np.unique(idx)
            remap = np.full(len(local), -1, dtype=np.int64)
            remap[used] = np.arange(len(used))
            local, fids, tris = local[used], fids[used], remap[idx].reshape(-1, 3)
            sel = sel[used]

        table = {"uid": [buildings[b][0] for b in blist]}
        for pname, _, _ in PROPS:
            table[pname] = np.array([buildings[b][3][pname] for b in blist], dtype=np.float32)
        table["lod"] = np.full(len(blist), lod, dtype=np.uint8)
        table["n_triangles"] = np.array([buildings[b][3]["_src_triangles"] for b in blist], dtype=np.uint32)

        glb_name = cell_name + ".glb"
        nbytes = write_glb(os.path.join(out_dir, glb_name), local, tris, fids, table, schema,
                           "%s_lod%d_%s" % (tile_name, lod, cell_name), compress)
        total_bytes += nbytes
        total_tris += len(tris)
        region = [math.radians(float(lon_c.min())), math.radians(float(lat_c.min())),
                  math.radians(float(lon_c.max())), math.radians(float(lat_c.max())),
                  float(h_c.min()), float(h_c.max())]
        transform = [float(R[0, 0]), float(R[1, 0]), float(R[2, 0]), 0.0,
                     float(R[0, 1]), float(R[1, 1]), float(R[2, 1]), 0.0,
                     float(R[0, 2]), float(R[1, 2]), float(R[2, 2]), 0.0,
                     float(C[0]), float(C[1]), float(C[2]), 1.0]
        tiles[key] = {"region": region, "transform": transform, "uri": "lod%d/%s" % (lod, glb_name),
                      "buildings": len(blist), "triangles": int(len(tris)), "bytes": nbytes}
        for b in blist:
            uid, Vb, Tb, props = buildings[b]
            csv_rows.append([uid, lod, cell_name] + [props[p] for p, _, _ in PROPS] +
                            [round(float(centre[b, 0]), 2), round(float(centre[b, 1]), 2), len(Tb), len(Vb),
                             round(float(zmin[b]), 3)])
    print("  lod%d: %d cells, %d GLBs, %.1f MB, %d triangles, %d vertices in %.1fs" %
          (lod, len(cells), len(tiles), total_bytes / 1e6, total_tris, len(V), time.time() - t0))

    heights = np.array([b[3]["height"] for b in buildings])
    heights = heights[heights != NODATA]
    areas = np.array([b[3]["footprint_area"] for b in buildings])
    areas = areas[areas != NODATA]
    hist = np.histogram(heights, bins=HEIGHT_BINS)[0].tolist() if len(heights) else [0] * (len(HEIGHT_BINS) - 1)
    tallest = max(buildings, key=lambda b: b[3]["height_max"] if lod == 2 else b[3]["height"])
    stats = {"buildings": nb, "cells": len(tiles), "triangles": int(total_tris), "vertices": int(len(V)),
             "bytes": int(total_bytes),
             "height_mean_m": float(heights.mean()) if len(heights) else None,
             "height_median_m": float(np.median(heights)) if len(heights) else None,
             "height_hist": hist,
             "footprint_area_total_m2": float(areas.sum()) if len(areas) else None,
             "tallest": {"uid": tallest[0], "height": tallest[3]["height_max"] if lod == 2 else tallest[3]["height"]},
             "bbox_bng": [float(V[:, 0].min()), float(V[:, 1].min()), float(V[:, 0].max()), float(V[:, 1].max())]}
    return tiles, stats, csv_rows


def relative_transform(parent, child):
    """3D Tiles composes transforms down the hierarchy, so a child's transform must be
    expressed relative to its parent's. Both are column-major 16-element lists."""
    Mp = np.array(parent, dtype=np.float64).reshape(4, 4).T
    Mc = np.array(child, dtype=np.float64).reshape(4, 4).T
    rel = np.linalg.inv(Mp) @ Mc
    return rel.T.ravel().tolist()


def union_region(regions):
    r = np.array(regions)
    return [float(r[:, 0].min()), float(r[:, 1].min()), float(r[:, 2].max()),
            float(r[:, 3].max()), float(r[:, 4].min()), float(r[:, 5].max())]


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", required=True, help="OS tile name, e.g. NS56SE")
    ap.add_argument("--lod1", required=True, help="folder with LoD1 CityJSON files")
    ap.add_argument("--lod2", default=None, help="folder with LoD2 CityJSON files (optional)")
    ap.add_argument("--out", required=True, help="output folder (tileset.json, lod1/, lod2/)")
    ap.add_argument("--cell", type=float, default=312.5, help="LoD2 cell size in metres")
    ap.add_argument("--coarse-cell", type=float, default=625.0, help="LoD1 cell size in metres (multiple of --cell)")
    ap.add_argument("--lod1-geometric-error", type=float, default=50.0,
                    help="geometricError of LoD1 tiles that have LoD2 children (controls when LoD2 loads)")
    ap.add_argument("--root-geometric-error", type=float, default=1000000.0,
                    help="root geometricError; huge so the root always refines to the LoD1 level "
                         "(the viewer's 'LoD1 only' mode relies on this)")
    ap.add_argument("--no-flatten", action="store_true", help="keep ODN heights instead of dropping buildings to height 0")
    ap.add_argument("--no-compress", action="store_true", help="skip meshopt compression (debugging)")
    ap.add_argument("--lod1-outline-tolerance", type=float, default=0.3,
                    help="Douglas-Peucker tolerance (m) for rebuilding LoD1 prisms from their roof outline; 0 = keep source mesh")
    ap.add_argument("--lod2-simplify", type=float, default=0.10,
                    help="error-bounded simplification tolerance (m) for LoD2 meshes; 0 = none")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = ap.parse_args()

    flatten = not args.no_flatten
    compress = not args.no_compress
    x0, y0 = os_tile_origin(args.tile)
    ratio = args.coarse_cell / args.cell
    if abs(ratio - round(ratio)) > 1e-9:
        sys.exit("--coarse-cell must be a multiple of --cell")
    ratio = int(round(ratio))
    if flatten:
        tg = TransformerGroup("EPSG:27700", "EPSG:4326", always_xy=True)
    else:
        tg = TransformerGroup("EPSG:7405", "EPSG:4979", always_xy=True)   # BNG + ODN height -> WGS84 3D
        if not tg.best_available:
            sys.exit("--no-flatten needs the OSTN15/OSGM15 grids: tools/venv/bin/pyproj sync --file uk_os_OSTN15_NTv2_OSGBtoETRS "
                     "&& tools/venv/bin/pyproj sync --file uk_os_OSGM15_GB")
    tr = tg.transformers[0]
    crs_note = "%s (accuracy %s m)" % (tr.description, tr.accuracy)
    schema = build_schema()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    print("[%s] origin E%d N%d, %s" % (args.tile, x0, y0, "flattened" if flatten else "ODN heights"))

    b1 = load_lod(args.lod1, 1, args.workers, args.lod1_outline_tolerance)
    lod1_tiles, s1, rows1 = build_level(b1, 1, args.coarse_cell, x0, y0, os.path.join(args.out, "lod1"),
                                        tr, flatten, compress, schema, args.tile)
    lod2_tiles, s2, rows2 = {}, None, []
    if args.lod2:
        b2 = load_lod(args.lod2, 2, args.workers)
        lod2_tiles, s2, rows2 = build_level(b2, 2, args.cell, x0, y0, os.path.join(args.out, "lod2"),
                                            tr, flatten, compress, schema, args.tile, args.lod2_simplify)

    # --- hierarchy --------------------------------------------------------------
    children_of = {}
    for key, t in lod2_tiles.items():
        parent = (key[0] // ratio, key[1] // ratio)
        children_of.setdefault(parent, []).append(t)
    lod1_nodes = []
    for key in sorted(set(lod1_tiles) | set(children_of)):
        t1 = lod1_tiles.get(key)
        kids = children_of.get(key, [])
        node = {"refine": "REPLACE"}
        regions = []
        if t1:
            regions.append(t1["region"])
            node["transform"] = t1["transform"]
            node["content"] = {"uri": t1["uri"]}
            node["extras"] = {"lod": 1, "buildings": t1["buildings"], "triangles": t1["triangles"]}
        if kids:
            regions += [k["region"] for k in kids]
            node["geometricError"] = args.lod1_geometric_error
            node["children"] = [{"boundingVolume": {"region": k["region"]}, "geometricError": 0.0,
                                 "transform": relative_transform(t1["transform"], k["transform"]) if t1 else k["transform"],
                                 "content": {"uri": k["uri"]},
                                 "extras": {"lod": 2, "buildings": k["buildings"], "triangles": k["triangles"]}}
                                for k in kids]
        else:
            node["geometricError"] = 0.0
        node["boundingVolume"] = {"region": union_region(regions)}
        lod1_nodes.append(node)

    root_region = union_region([n["boundingVolume"]["region"] for n in lod1_nodes])
    tileset = {
        "asset": {"version": "1.1", "generator": "ubdc build_tileset.py",
                  "extras": {"tile": args.tile, "source": "UBDC Glasgow 3D building model, doi:10.5281/zenodo.15000747",
                             "crs_note": "EPSG:27700 -> WGS84 via " + crs_note,
                             "flattened": flatten, "cell_m": args.cell, "coarse_cell_m": args.coarse_cell,
                             "lod1_outline_tolerance_m": args.lod1_outline_tolerance, "lod2_simplify_m": args.lod2_simplify}},
        "geometricError": args.root_geometric_error,
        "root": {"boundingVolume": {"region": root_region}, "geometricError": args.root_geometric_error,
                 "refine": "REPLACE", "children": lod1_nodes},
    }
    with open(os.path.join(args.out, "tileset.json"), "w") as f:
        json.dump(tileset, f, separators=(",", ":"))

    with open(os.path.join(args.out, "buildings.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["uid", "lod", "cell"] + [p for p, _, _ in PROPS] + ["centre_e", "centre_n", "n_triangles", "n_vertices", "zmin_odn"])
        w.writerows(rows1 + rows2)

    index = {"tile": args.tile, "origin_bng": [x0, y0], "root_region": root_region, "flattened": flatten,
             "compressed": compress, "crs_note": crs_note, "height_bins": HEIGHT_BINS[:-1],
             "lod1": s1, "lod2": s2, "built": time.strftime("%Y-%m-%d %H:%M")}
    with open(os.path.join(args.out, "index.json"), "w") as f:
        json.dump(index, f, indent=1)
    print("[%s] done in %.1fs -> %s" % (args.tile, time.time() - t0, args.out))


if __name__ == "__main__":
    main()
