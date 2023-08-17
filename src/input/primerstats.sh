#!/bin/bash

path=$(pwd)
# shellcheck disable=SC2164
cd $path/primers
bowtie2-build $path/primers/covid.fasta covid_index
bowtie2 -x covid_index -f ArticV4.1.fasta -S primersmapped.sam --local
samtools view -bS primersmapped.sam > primersmapped.bam
samtools view -bS primersmapped.sam > primersmapped.bam ; samtools sort primersmapped.bam -o primersmapped_sorted.bam
bedtools bamtobed -i primersmapped_sorted.bam > primers.bed
samtools ampliconstats -l 20000 primers.bed $path/input/aln.virus.sorted_NC_045512.2.sorted.bam > stats.txt
rm primersmapped.bam ; rm primersmapped.sam ; rm primersmapped_sorted.bam