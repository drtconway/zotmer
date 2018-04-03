#!/bin/bash

wd=tmp/016
mkdir -p ${wd}

S=123

zot hgvs-gen -g ~/data/hg19 -N 100 -S ${S} etc/panel.bed > ${wd}/all.txt
zot hgvs-find -X -g ~/data/hg19 -f ${wd}/all.txt ${wd}/all.idx

v=$(head -n 1 < ${wd}/all.txt)
msg="one-among-many ${v}"

zot spikein -g ~/data/hg19 -b etc/panel.bed -N 400000 -V 0.5 -S ${S} -z ${wd}/reads ${v}
zot hgvs-find -o ${wd}/output.txt ${wd}/all.idx ${wd}/reads_*.fastq.gz
cat > ${wd}/rules.py << EOF
- if: n == '$m'
  then:
    - res == 'wt/mut'
    - 0.3 <= float(wtVaf) and float(wtVaf) <= 0.7
    - 0.3 <= float(mutVaf) and float(mutVaf) <= 0.7
- if: n != '0'
  then:
    - res == 'wt'
EOF
./bin/check ${wd}/rules.py ${wd}/output.txt > ${wd}/result.txt

echo "$0" $(cat ${wd}/result.txt) ${msg}

cleanup=yes
if [ $(cat ${wd}/result.txt) == 'ACK' ] && [ ${cleanup} == 'yes' ]
then
    rm -rf ${wd}
fi
