#!/bin/bash
# bash create_input_file.sh D3N N3D

forward=$1
backward=$2
fepout_file=${forward}-${backward}.csv
if [[ -f ${fepout_file} ]]; then rm ${fepout_file}; fi

for file in $(find `pwd` -type f -name "*.fepout" | sort -k 2); do
    if [[ $file == *"/bound/"* ]]; then
        state="bound"
    elif [[ $file == *"/free/"* ]]; then
        state="free"
    else
        direction="#" # comment line
    fi

    if [[ $file == *"/$forward/"* ]]; then
        direction="forward"
    elif [[ $file == *"/$backward/"* ]]; then
        direction="backward"
    else
        continue
    fi

    # one file per line
    echo "$state,$direction,$file" >>"${fepout_file}"
done

echo Done, and run with:
echo python fepout-analysis.py -f ${fepout_file}
