# Glasgow 3D City Model dashboard

Interactive web dashboard for the UBDC Glasgow 3D building model
(Li, Zhao et al. 2025, [doi:10.5281/zenodo.15000747](https://doi.org/10.5281/zenodo.15000747), CC BY 4.0):
119,765 LoD1 buildings across 23 OS 5 km tiles, 38,165 of them with LoD2 roof detail in the city centre.

Static site, no backend: a CesiumJS viewer streaming 3D Tiles, plus KPIs, a height histogram and a
per-tile coverage table. Served by GitHub Pages at **https://urbanbigdatacentre.github.io/3d-glasgow/**
(every push to `main` redeploys; the `tiles/` folder is committed, so no build step runs on GitHub).

## Layout

```
index.html                 dashboard page (CesiumJS from CDN, no build step)
tiles/tileset.json         master tileset -> one external tileset per OS tile
tiles/summary.json         city-wide statistics consumed by the page
tiles/<OSTILE>/            tileset.json, lod1/*.glb (625 m cells), lod2/*.glb (312.5 m cells), buildings.csv, index.json
tools/build_tileset.py     CityJSON (one building per file, EPSG:27700) -> 3D Tiles 1.1 for one OS tile
tools/assemble.py          master tileset + summary.json from all tiles
tools/build_all.sh         rebuild everything (23 tiles) and assemble
tools/meshopt.py           ctypes binding to meshoptimizer (compression, simplification)
tools/build_meshopt.sh     builds tools/vendor/libmeshoptimizer.dylib from source (clang++, ~10 s)
```

The raw data (CityJSON / OBJ / FileGDB, 5.5 GB, on Zenodo) lives **outside** this folder and is not committed.

## How the tiles are made

Per OS tile, `build_tileset.py` produces a two-level tileset with `REPLACE` refinement:

* **LoD1** (625 m cells) is the coarse level shown from afar. Every LoD1 prism is rebuilt from its roof
  outline: Douglas-Peucker at 0.3 m, then re-extruded with earcut. This matters for tiles such as
  NS56NE whose source outlines are raster traces with a vertex every 0.5 m (1,000 triangles per prism).
  Prisms whose rebuilt footprint area deviates by more than 5 % keep the source mesh.
* **LoD2** (312.5 m cells) replaces LoD1 within roughly 2.5 km of the camera. Meshes are decimated with
  an absolute error bound of 0.10 m (`--lod2-simplify`), which removes about half the triangles.
* GLBs are glTF 2.0 with int16-quantised positions (`KHR_mesh_quantization`), meshopt-compressed vertex
  and index streams (`EXT_meshopt_compression`), a per-vertex feature id (`EXT_mesh_features`) and a
  building property table (`EXT_structural_metadata`: uid, mean/max/min/P90 height, ground level,
  footprint area and perimeter, LoD, source triangle count). No normals are stored; the viewer shades
  from screen-space derivatives in a custom shader.
* Heights: buildings are dropped to ellipsoid height 0 (`--flatten` default) because the viewer has no
  terrain yet. Rebuild with `--no-flatten` and add the ODN-to-ellipsoid geoid offset when adding terrain.
* Coordinates: EPSG:27700 to WGS84 with the transformation available to pyproj on the build machine
  (recorded in each `tileset.json` `asset.extras.crs_note`; the 2 m Helmert unless the OSTN15 grid is installed).

## Rebuilding

One-time setup (macOS, Python 3.8+ with numpy, pyproj, shapely):

```bash
sh tools/build_meshopt.sh                       # clones zeux/meshoptimizer and builds the dylib
pip3 install --target tools/vendor/py mapbox_earcut   # polygon triangulation (earcut)
rm -rf tools/vendor/py/numpy*                   # keep the environment's numpy
```

Then, with the Zenodo data unpacked one level above this folder (`../lod1_3d_building_model`, `../lod2_3d_building_model`):

```bash
sh tools/build_all.sh --workers 8               # ~10 minutes; writes tiles/
```

or a single tile:

```bash
python3 tools/build_tileset.py --tile NS56SE \
  --lod1 ../lod1_3d_building_model/cityjson/lod1_3d_building_model_json/NS56SE \
  --lod2 ../lod2_3d_building_model/cityjson/lod2_3d_building_model_json/NS56SE \
  --out tiles/NS56SE
python3 tools/assemble.py tiles
```

After regenerating, bump `DATA_VERSION` in `index.html` so browsers refetch the tiles.

## Local preview

```bash
python3 -m http.server 8765 --directory .
```

then open http://localhost:8765/.

## Dashboard notes

* Basemap: OpenStreetMap standard tiles (fine for light use with attribution). For a production
  deployment consider the OS Maps API (free OS Data Hub key) or another keyed provider.
* Colours: building height classes use a single-hue ordinal ramp; LoD mode uses two categorical hues.
* Zenodo download/view counts are fetched live from the Zenodo API, with a snapshot fallback.
