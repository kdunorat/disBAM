#!/bin/bash

samtools view -h -F 4 input/$1 | awk '(index($6, "S") != 0) || $1 ~ /@/ {print}' | samtools view -bS - > input/soft_cliped.bam

samtools view -o input/soft_cliped.sam input/soft_cliped.bam

rm input/soft_cliped.bam
