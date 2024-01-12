#!/bin/bash

# declare the arrays for the files and the sorting
declare -A files
declare -A sorting

# get a list of filenames into it, saving number without 0's as key
for file in *; do
    fnum=$(echo "$file" | tr -d -c 0-9 | sed 's/^0*//')
    files[$fnum]="$file"
    sorting[$fnum]=$fnum
done

# sort the array by its numeric key values
IFS=$'\n' sorted=($(sort -n <<<"${sorting[*]}"))
unset IFS

# check for user input and if its numerical
if [[ $1 =~ ^-?[0-9]+$ ]]; then
    # iterate through the array
    for i in "${sorted[@]}"; do
        # only handle files above user input number
        if [[ $i -gt $1 ]]; then
            # execute your sql here, echo is just for debugging
            echo ${files[$i]}
        fi
    done
else
    echo "Please supply a number as argument"
    exit 1
fi
