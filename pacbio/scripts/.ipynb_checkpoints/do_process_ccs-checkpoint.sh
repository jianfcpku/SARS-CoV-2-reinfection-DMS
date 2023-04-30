#!/bin/bash

ncpu=1
targ=SARS-CoV-2-BA5
runs=../pacbio_runs.csv

IN=../ccsfastq
OUT=../outputs

python process_ccs.py --seqsfile=../pacbio_amplicons.gb \
                      --configfile=../parse_specs.yaml \
                      --threads=$ncpu \
                      --pbruns=$runs \
                      --output=$OUT \
                      --target=$targ