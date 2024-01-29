#!/bin/bash

# Usage:
# Example 1:  
# :~$ bash rename-png-files-in-order.sh 0
# This will rename files from starting at 1 
# Example 2:
# :~$ bash rename-png-files-in-order.sh 2
# This will rename files from starting at 3	

# After this call, the 'avi' video can be generated: 
# :~$ ffmpeg -framerate 25 -i %03d.png output.avi

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
            mv ${files[$i]} $(printf %03d.%s $i "png")
            # echo $i; echo ${files[$i]} # A lines just for debuging!!! 
        fi
    done
else
    echo "Please supply a number as argument"
    exit 1
fi


# declare the arrays for the files and the sorting
declare -A files
declare -A sorting

# get a list of filenames into it, saving number without 0's as key
for file in *.png_* 
do
    fnum=$(echo "$file" | tr -d -c 0-9 | sed 's/^0*//')
    echo "$fnum"
    files[$fnum]="$file"
    sorting[$fnum]=$fnum
done

# sort the array by its numeric key values
IFS=$'\n' sorted=($(sort -n <<<"${sorting[*]}"))
unset IFS

# check for user input and if its numerical
if [[ $1 =~ ^-?[0-9]+$ ]]; then
    # iterate through the array
    for i in "${sorted[@]}"
    do
        # only handle files above user input number
        if [[ $i -gt $1 ]]; then
            cp -v ${files[$i]} $(printf %03d.%s $i "png")
            echo ${files[$i]}
        fi
    done
else
    echo "Please supply a number as argument"
    exit 1
fi

