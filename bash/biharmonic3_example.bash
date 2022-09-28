#!/bin/bash
# Declare an array of string with type
declare -a Runname="../../../build/bin/biharmonic3_example"
declare -a Geometries=("1121")
# "1121" "1707" "1708" "1402" "1022")
declare -a Methods=(0 1 2)

mkdir -p Output

for geo in ${Geometries[@]}; do
    for (( p=2; p<6; p++)) do
        for m in ${Methods[@]}; do
            # Set smoothness
            declare s=$(($p-1))
            # Options per method
            if (($m==0))
            then
                if (($p<3))
                then
                    continue
                fi
            elif (($m==2))
            then
                if (($p>2))
                then
                    continue
                fi
            fi

            echo "$Runname" -m $m -p $p -s $s -g g"$geo" -r 4

            eval "$Runname" -m $m -p $p -s $s -g g"$geo" -r 4 -w Output/Biharmonic_m"$m"_p"$p"_s"$s"_g"$geo".csv > Output/Biharmonic_m"$m"_p"$p"_s"$s"_g"$geo"
        done
    done
done
