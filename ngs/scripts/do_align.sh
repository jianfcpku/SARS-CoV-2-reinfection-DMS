#!/bin/bash
set -e

IN=../fastq
OUT=../outputs

INFO=../sampleInfo.csv
scripts=.

mkdir -p $OUT/raw
mkdir -p $OUT/align

tmp=$OUT/tmp

mkdir -p $tmp/sh
mkdir -p $tmp/outputs

cat $INFO | tail -n +2 | while read x
do
    sample=$(echo $x | cut -d',' -f 1)
    lib=2mutBA5Tmerged
    antibody=$(echo $x | cut -d',' -f 2)
    file=$(echo $x | cut -d',' -f 3)
    
    if [ ! -f $file ];then
        echo "Empty $sample"
        continue
    fi

    wt=../SARS-CoV-2-BA5.fasta
    bclen=26
    table=../../pacbio/outputs/codon_variant_table_2mutBA5Tmerged.csv

    alignout=$OUT/align/${antibody}-${sample}_${lib}
    mkdir -p $alignout

    echo "#!/bin/bash" > $tmp/sh/_02-align-${sample}.tmp.sh
    echo "python $scripts/01-align.py -i $file -o $alignout -t $table -w $wt -b $bclen" >> $tmp/sh/_02-align-${sample}.tmp.sh
    echo $sample,wt:$wt,bclen:$bclen,table:$table

    bash $tmp/sh/_02-align-${sample}.tmp.sh > $tmp/outputs/_02-align-${sample}.o.txt
done
