# SARS-CoV-2-reinfection-DMS
This repo contains DMS data analysis code and processed data of mAbs involved in the manuscript **Repeated Omicron infection alleviates SARS-CoV-2 immune imprinting** ([Yisimayi et al. bioRxiv 2023](https://www.biorxiv.org/content/10.1101/2023.05.01.538516))

## Processing raw sequencing data

Raw PacBio sequencing data of variants in the library, and NGS data of barcodes before and after antibody screening have been uploaded to [China National GeneBank](db.cngb.org) with accession number [CNP0004294](https://db.cngb.org/search/project/CNP0004294). The data are also available on [Genome Sequence Archive (GSA) of China National Center for Bioinformation](https://ngdc.cncb.ac.cn/gsa) with Project accession [PRJCA020116](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA020116).

Re-analysis of the raw sequencing data is time-consuming and not necessary for most purposes of reusing our results. You can skip this section if you would like to directly use our processed escape scores for each antibody.

### PacBio sequencing 

We built two independent mutant libraries based on BA.5 RBD, named **2mutBA5lib2T** and **2mutBA5lib4T** respectively. Each of the two libraries were sequenced using PacBio SMRT Platform. The .fastq.gz files (ccs) the should be moved to `pacbio/ccsfastq/2mutBA5lib2T.fastq.gz` and `pacbio/ccsfastq/2mutBA5lib4T.fastq.gz`, respectively.

We used the [pipeline from J. Bloom lab](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS/) to generate the barcode-variant table. You can enter the `pacbio/scripts` directory and run `do_process_ccs.sh`.

The two libraries are then merged for the following experiments for antibody DMS. The merged library is referred to as **2mutBA5Tmerged**, and the corresponding barcode-variant table is `pacbio/outputs/codon_variant_table_2mutBA5Tmerged.csv`. The barcodes that conflict in the two libraries are dropped.

### NGS data alignment

To generate the escape scores of a sample (antibody), two NGS profiles corresponding to the barcodes before and after antibody screening are needed, which are called the reference barcodes and sample barcodes, respectively. Multiple samples may share the same references as they were tested in the same batch. Compressed fastq files of both samples and references downloaded from the database should be moved to `ngs/fastq`.

Run `ngs/scripts/do_align.sh` to get a variant table for each sample or reference. The results are in `ngs/outputs/align` and will be used to calculate the escape scores.

### Calculate escape scores and merge

We will get two files `effects_df.csv` and `effects_df_no_filter.csv` for each sample after the calculation, corresponding to the mutation escape scores with or without the filter on RBD expression. The filtered data will be used for clustering to reduce noise, and the unfiltered data will be used for calculating preferences and visualizing the average DMS profiles for clusters to retain as many mutations as possible.

(will be completed before publication)

## Antibody clustering

All of the processed DMS profiles (filtered by RBD expression) are now merged into one `.csv` file (`antibody_dms_merge.csv.gz`). Please uncompress it before running the following procedures. 

Scripts for merging and clustering the antibody DMS profiles are incorporated in the Jupyter Notebook `clustering/clustering.ipynb`. 

## Calculating the preference of mutations

Neutralizing activities (pseudovirus IC50) of the antibodies against BA.5 and XBB.1.5 (`antibody_info.csv`), and [DMS on ACE2 binding and RBD expression of BA.2 from J. Bloom lab](https://doi.org/10.1371/journal.ppat.1010951) (`calculation/BA2_bind_expr.csv`, the same as [that in previous study](https://github.com/jianfcpku/convergent_RBD_evolution/blob/main/bind_expr/bind_expr_BA2.csv) are needed for the calculation. You can also use other DMS datasets on ACE2 binding and RBD expression with the same format based on BA.5 or XBB.1.5 which might be published in the future to get potentially better results.

We use the unfiltered merged DMS profiles after denoising according to medians (`antibody_dms_merge_no_filter_clean.csv`)

See `calculation/calculate_preference.ipynb` for detail. 

The pipeline is similar to [that in the previous study](https://github.com/jianfcpku/convergent_RBD_evolution), except for a modified weighting strategy for codon constraints, and an additional weight to correct the proportion of Omicron-specific (or cross-reactive) antibodies in the dataset, making the ratio the same as the unbiased ratio determined by ELISA (89% cross-reactive in BA.5/BF.7 BTI, and 51% in reinfection).
