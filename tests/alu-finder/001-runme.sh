#!/bin/bash

wd=tmp/001
mkdir -p ${wd}

S=17

etc/mkAluIns ${S} etc/BRCA2-ex3.bed etc/some-alus.fa > ${wd}/var.txt 2> ${wd}/exp.txt

zot spikein -v -g ~/data/hg19 -b etc/BRCA2-exons.bed -N 50000 -S $((${S} + 1)) -V 1.0 -f ${wd}/var.txt -z ${wd}/reads

zot alu-finder -g ~/data/hg19 etc/BRCA2-exons.bed ${wd}/reads*.fastq.gz > ${wd}/out.txt
