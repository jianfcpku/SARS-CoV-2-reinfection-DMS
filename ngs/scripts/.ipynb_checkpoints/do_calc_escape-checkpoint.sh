#!/bin/bash

set -e

IN=../outputs
INFO=../sampleInfo.csv
scripts=.
tmp=$IN/tmp

mkdir -p $IN/calc

cat $INFO | tail -n +2 | awk 'BEGIN{FS=","}{if($2 != "ref") print $0}' | while read x
do
    file=$(echo $x | cut -d',' -f 3)
    sample=$(echo $x | cut -d',' -f 1)
    lib=2mutBA5Tmerged
    antibody=$(echo $x | cut -d',' -f 2)
    ref=$(echo $x | cut -d',' -f 4)
    
    echo "$sample . Library: $lib . Antibody: $antibody"
    wt=../SARS-CoV-2-BA5.fasta
    table=$(cat $libinfo | grep ^${lib}, | cut -d',' -f 5)
    
    singlemut=$(cat $libinfo | grep ^${lib}, | cut -d',' -f 6)
    bind=$(cat $libinfo | grep ^${lib}, | cut -d',' -f 7)
    expr=$(cat $libinfo | grep ^${lib}, | cut -d',' -f 8)
    
    echo "    ${antibody}-${sample}_${lib}. Using reference $ref"

    cat << EOF > $tmp/sh/_05_${antibody}-${sample}_${lib}.tmp.sh
#!/bin/bash
python $scripts/02-call_escape_score.py \
           -r $ref \
           -e $IN/align/${antibody}-${sample}_${lib}/${antibody}-${sample}_${lib}_variant_counts.csv \
           -exp MACS \
           -w $wt -t $table \
           -o $IN/calc \
           -p SARS-CoV-2-BA5 -S 330 \
           --mutbindexpr=$singlemut --variant_bind=$bind --variant_expr=$expr\
           --exprmin=-1 --bindmin=-999
EOF

    bash $tmp/sh/_05_${antibody}-${sample}_${lib}.tmp.sh > $tmp/outputs/_05_${antibody}-${sample}_${lib}.o.txt
done
