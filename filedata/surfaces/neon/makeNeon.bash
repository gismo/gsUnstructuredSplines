#!/bin/bash
# Declare an array of string with type
declare -a Runname="./bin/refit_patches"
declare -a Meshes=("2x2" "3x3" "4x4")
declare -a BVP="../extensions/gsKLShell/filedata/pde/car_bvp.xml"
declare -a k_vec=(1 2 3)
declare -a N=2000
declare -a Geometries=("../filedata/surfaces/neon_side_split_pillar" "../filedata/surfaces/neon_side_split_pillar_extended" "../filedata/surfaces/neon_side_split")
declare -a d=2
declare -a t=1

for geometry in ${Geometries[@]}; do
    for k in ${k_vec[@]}; do
        declare -a kpp=$((k+1))
        declare -a name="$kpp"x"$kpp"
        echo "$Runname" -k $k -N $N -d $d -t $t "$geometry".3dm
        eval "$Runname" -k $k -N $N -d $d -t $t "$geometry".3dm
        # eval "$Runname" -G ../filedata/surfaces/neon_side_fit_"$mesh".xml -B $BVP -m $method -r $r -p $p -s $s --plot -N 16 -S 1e-5 --write -o ModalResults/Car_"$mesh"_r"$r"_p"$p"_s"$s"_m"$method"
        echo mv final.xml $"$geometry"_fit_"$name".xml
        eval mv final.xml $"$geometry"_fit_"$name".xml
    done
done
