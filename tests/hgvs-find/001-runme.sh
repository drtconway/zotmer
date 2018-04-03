#!/bin/bash

wd=tmp/001
mkdir -p ${wd}

S=17

v=$(zot hgvs-gen -g ~/data/hg19 -t del -S ${S} etc/TP53-exons.bed)
msg="wild-type del ${v}"

zot spikein -g ~/data/hg19 -b etc/TP53.bed -N 10000 -V 0.5 -S 17 -z ${wd}/reads

zot hgvs-find -X -g ~/data/hg19 ${wd}/variant.idx ${v}

zot hgvs-find -o ${wd}/output.txt ${wd}/variant.idx ${wd}/reads_*.fastq.gz

cat > ${wd}/rules.py << EOF
- if: n == '0'
  then:
    - res == 'wt'
    - 0.7 <= float(wtVaf) and float(wtVaf) <= 1.0
EOF
./bin/check ${wd}/rules.py ${wd}/output.txt > ${wd}/result.txt

echo "$0" $(cat ${wd}/result.txt) ${msg}

cleanup=yes
if [ $(cat ${wd}/result.txt) == 'ACK' ] && [ ${cleanup} == 'yes' ]
then
    rm -rf ${wd}
fi
