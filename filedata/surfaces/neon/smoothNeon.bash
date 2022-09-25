#!/bin/bash
# Declare an array of string with type
declare -a Runname="./bin/smoothInterfaces"
declare -a Meshes=("2x2" "3x3" "4x4")
declare -a Geometries=("../filedata/surfaces/neon_side_split_pillar" "../filedata/surfaces/neon_side_split_pillar_extended" "../filedata/surfaces/neon_side_split")

for geometry in ${Geometries[@]}; do
    for mesh in ${Meshes[@]}; do
        echo "$Runname" -i "$geometry"_fit_"$mesh".xml -o "$geometry"_fit_"$mesh"_smooth.xml
        eval "$Runname" -i "$geometry"_fit_"$mesh".xml -o "$geometry"_fit_"$mesh"_smooth.xml
    done
done
