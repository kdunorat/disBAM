#!/bin/bash

path=$(pwd)

samtools view -h -F 4 $path/input/$1 | awk '(index($6, "S") != 0) || $1 ~ /@/ {print}' | samtools view -bS - > $path/input/soft_cliped.bam

samtools view -o $path/input/soft_cliped.sam $path/input/soft_cliped.bam

rm $path/input/soft_cliped.bam
