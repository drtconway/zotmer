#!/bin/bash

wd=tmp/018
mkdir -p ${wd}

v="NC_000017.10:g.41243867_41243868insCCCTAGTAGACTGAGAAGGTATATTGTTTACTTTACCAAATAACAAGTGT"
msg="long ins ${v}"

zot spikein -g ~/data/hg19 -b etc/BRCA1.bed -N 50000 -V 0.5 -S 17 -z ${wd}/reads ${v}

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
