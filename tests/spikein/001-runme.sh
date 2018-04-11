#!/bin/bash

wd=tmp/001
mkdir -p ${wd}

S=17

./bin/mkRandomSequence ${S} 4096 > ${wd}/chr13.txt
echo '>chr13' > ${wd}/chr13.fa
cat ${wd}/chr13.txt >> ${wd}/chr13.fa

echo "chr13 1024 4098 NONE" > ${wd}/none.bed

zot spikein -e 0 -g ${wd} -b ${wd}/none.bed -N 100 -S $(($S + 1)) ${wd}/reads

./bin/checkReads ${wd} ${wd}/none.bed ${wd}/reads_*
