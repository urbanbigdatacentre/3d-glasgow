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
* Heights: with `--no-flatten` (the deployed configuration) buildings keep their ODN heights, converted
  to WGS84 ellipsoidal heights with OSTN15 + OSGM15, so they sit on the terrain tiles. The default
  flattens every building to ellipsoid height 0 for a terrain-less globe.
* Coordinates: EPSG:27700 to WGS84 with the transformation available to pyproj on the build machine
  (recorded in each `tileset.json` `asset.extras.crs_note`; the 2 m Helmert unless the OSTN15 grid is installed).

## Terrain

`tools/build_terrain.py` turns the UBDC 0.5 m LiDAR DTM (Zenodo record 13273124, `DTM_5x5km.zip`,
5.3 GB; `tools/fetch_dtm.py` downloads the 23 tiles individually) into Cesium quantized-mesh tiles
under `tiles/terrain/` (TMS geodetic tiling, levels 0-14, Delatin meshing with a 0.5 m error bound
at level 14). Outside the DTM, Copernicus DEM GLO-30 provides the surroundings, feathered over 400 m
so the survey edge is not a cliff. Heights are converted from ODN to WGS84 ellipsoidal with PROJ's
OSGM15 grid (about +54 m in Glasgow); the building tiles are built with `--no-flatten` so they use
the same conversion and sit on the terrain.

```bash
tools/venv/bin/python tools/fetch_dtm.py --out ../DTM_5x5km
tools/venv/bin/python tools/build_terrain.py --dtm ../DTM_5x5km --out tiles/terrain --max-level 14 --base-error 0.5
```

## Rebuilding

One-time setup (macOS; a Python 3.9+ virtualenv in `tools/venv`, gitignored):

```bash
python3 -m venv tools/venv
tools/venv/bin/pip install numpy pyproj shapely scipy mapbox_earcut rasterio pydelatin quantized-mesh-encoder
tools/venv/bin/pyproj sync --file uk_os_OSTN15_NTv2_OSGBtoETRS   # OS horizontal grid
tools/venv/bin/pyproj sync --file uk_os_OSGM15_GB                 # OS geoid grid (ODN -> ellipsoid)
sh tools/build_meshopt.sh                       # clones zeux/meshoptimizer and builds the dylib
```

Then, with the Zenodo data unpacked one level above this folder (`../lod1_3d_building_model`, `../lod2_3d_building_model`):

```bash
PYTHON=tools/venv/bin/python sh tools/build_all.sh --no-flatten --workers 8   # ~5 minutes; writes tiles/
```

(`--no-flatten` keeps real heights for use with the terrain; without it buildings are dropped to
ellipsoid height 0 for a flat globe.)

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

* Basemaps: OpenStreetMap (colour or greyscale) needs no key and is fine for light use with attribution.
  Two optional keys at the top of the script in `index.html` unlock more: `OS_MAPS_API_KEY`
  (free OS OpenData plan at osdatahub.os.uk: OS Light / Road / Outdoor styles) and `CESIUM_ION_TOKEN`
  (free Cesium ion account: Bing aerial imagery, and Cesium World Terrain if terrain is added later).
  The chosen basemap and the panel state are remembered per browser.
* Terrain: self-hosted quantized-mesh tiles from the UBDC DTM (see Terrain above); `CONFIG.TERRAIN_URL`
  in `index.html` points at them (set it to "" for a flat globe, together with tiles built without
  `--no-flatten`).
* Colours: building height classes use a single-hue ordinal ramp; LoD mode uses two categorical hues.
* Detail menu: *auto* swaps LoD2 roofs in within ~2.5 km, *LoD1 only* never does (fast), *LoD2 everywhere*
  loads all 1,011 LoD2 cells at once (~20 MB). A line under the controls shows how many LoD1/LoD2 tiles
  are on screen.
* Terrain controls: relief shading (the globe is lit by the same fixed light as the buildings) and a
  1×/2×/3× vertical exaggeration relative to sea level, which lifts buildings consistently with the ground.
  Clicking open ground reports its height above sea level (ODN) from the DTM.
* Zenodo download/view counts are fetched live from the Zenodo API, with a snapshot fallback.
