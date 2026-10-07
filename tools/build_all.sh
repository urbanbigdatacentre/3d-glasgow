#!/bin/sh
# Rebuild every OS tile (4 core tiles with LoD1+LoD2, 19 outer tiles LoD1 only), then assemble.
# Run from the dashboard folder: sh tools/build_all.sh   (OUT=tiles_new sh tools/build_all.sh to build elsewhere)
set -e
OUT=${OUT:-tiles}
DATA=..
L1=$DATA/lod1_3d_building_model/cityjson/lod1_3d_building_model_json
L1B=$DATA/lod1_3d_building_model/cityjson/lod1_3d_building_model_part2_json
L2=$DATA/lod2_3d_building_model/cityjson/lod2_3d_building_model_json
for T in NS56NE NS56SE NS66NW NS66SW; do
  python3 tools/build_tileset.py --tile $T --lod1 $L1/$T --lod2 $L2/$T --out $OUT/$T "$@"
done
for T in $(ls $L1B); do
  python3 tools/build_tileset.py --tile $T --lod1 $L1B/$T --out $OUT/$T "$@"
done
python3 tools/assemble.py $OUT
