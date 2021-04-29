#!/bin/bash

file1=$1
file2=$2

if [ "$file1" = "$file2" ]
then
    echo "Files have the same content"
else
    echo "Files have NOT the same content"
fi

