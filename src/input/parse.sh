#!/bin/bash

samtools view -h -F 4 src/input/$1 | awk '(index($6, "S") != 0) || $1 ~ /@/ {print}' | samtools view -bS - > src/input/soft_cliped.bam

samtools view -o src/input/soft_cliped.sam src/input/soft_cliped.bam

rm src/input/soft_cliped.bam
