#!/bin/bash

wd=tmp/002
mkdir -p ${wd}

S=17

./bin/mkRandomSequence ${S} 4096 > ${wd}/chr13.txt
echo '>chr13' > ${wd}/chr13.fa
cat ${wd}/chr13.txt >> ${wd}/chr13.fa

v='chr13:g.1536_2047del'
./bin/applyVariant ${wd}/chr13.fa ${v} > ${wd}/mut.fa

echo "chr13 1024 4098 NONE" > ${wd}/none.bed

zot spikein -e 0 -g ${wd} -b ${wd}/none.bed -N 100 -S $(($S + 1)) -V 1.0 ${wd}/reads ${v}

if ./bin/checkReads ${wd}/mut.fa ${wd}/reads_* > ${wd}/errors.txt
then
    echo "ACK" > ${wd}/result.txt
else
    echo "NAK" > ${wd}/result.txt
fi
