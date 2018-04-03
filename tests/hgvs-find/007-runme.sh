#!/bin/bash

wd=tmp/007
mkdir -p ${wd}

S=17

v=$(zot hgvs-gen -g ~/data/hg19 -t ins -S ${S} etc/TP53-exons.bed)
msg="heterozygous ins ${v}"

zot spikein -g ~/data/hg19 -b etc/TP53.bed -N 10000 -V 0.5 -S 21 -z ${wd}/reads ${v}

zot hgvs-find -X -g ~/data/hg19 ${wd}/variant.idx ${v}

zot hgvs-find -o ${wd}/output.txt ${wd}/variant.idx ${wd}/reads_*.fastq.gz

cat > ${wd}/rules.py << EOF
- if: n == '0'
  then:
    - res == 'wt/mut'
    - 0.4 <= float(wtVaf) and float(wtVaf) <= 0.6
    - 0.4 <= float(mutVaf) and float(mutVaf) <= 0.6
EOF
./bin/check ${wd}/rules.py ${wd}/output.txt > ${wd}/result.txt

echo "$0" $(cat ${wd}/result.txt) ${msg}

cleanup=yes
if [ $(cat ${wd}/result.txt) == 'ACK' ] && [ ${cleanup} == 'yes' ]
then
    rm -rf ${wd}
fi
