# ONTMethCalling

This repository contains scripts to call DNA methylation from ONT data and  perform statistical comparisions between a smoker and non-smoker. These scripts contain the analyses presented in 

***Evaluation of nanopore sequencing for epigenetic epidemiology: a comparison with DNA methylation microarrays. Flynn et al. 2022*****

Notes:


- The unix scripts are run in concert with a config file that contains the relevant filepaths as bash  variables. An example config file can be found in 
```
Config/example.txt
```
- Scripts are separated into job submission scripts which are submitted to the HPC scheduler and the processing scripts which are called/executed by the job submission script. These can be found in 
```
Scripts/
```
and
```
JobSubmission/
```
respectively.

## Description of scripts to process ONT data into DNA methylation levels

1_Alignment.sh

Takes fastq files from  ONT and merges, aligns with minimap2. Performs downstream prcosessing with samtools including, sorting, fitering to high quality primary alignments and indexing. 

2_CalcCoverage.sh

Takes aligned, filtered reads in bam file and use BEDtools to summarise read coverage across genome.

3_DNAMcallingNanopolish.sh

Takes fast5 files from ONT and aligned, indexed bam file from 1_Alignment.sh and uses nanopolish to call DNA methylation  and summarise as frequencys per position. 

## Description of R scripts to perform statistical comparisions

These scripts contain the code to create the figures included in the manuscript.

4_PlotCoverage.r

Takes bedgraph files from BEDtools (2_CalcCoverage.sh) and table of guide positions to generate figures and summary statistics of sequencing coverage across genome and within target regions.  

5_CompareArrayONT.r

Takes nanopolish DNA methylation frequencys and compares to matched data generated with the EPIC array.

6_CompareSamples.r

Takes nanopolish DNA methylation frequencys and performs statistical comparision between samples. Compares these results to smoking EWAS from large population cohort. 




