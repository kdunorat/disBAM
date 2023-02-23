#!/bin/bash

samtools view -h -F 4 $2/input/$1 | awk '(index($6, "S") != 0) || $1 ~ /@/ {print}' | samtools view -bS - > $2/input/soft_cliped.bam

samtools view -o $2/input/soft_cliped.sam $2/input/soft_cliped.bam

rm $2/input/soft_cliped.bam
