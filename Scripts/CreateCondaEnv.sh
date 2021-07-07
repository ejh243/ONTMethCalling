## set up conda environment for isoseq3 pipeline

module load Anaconda3/5.3.0

## had error "InvalidVersionSpecError: Invalid version spec: =2.7" which below fixed
conda config --remove channels conda-forge

conda create --name isoseq python=3.7 anaconda
source activate isoseq

conda install -n isoseq biopython
conda install -n isoseq -c http://conda.anaconda.org/cgat bx-python

conda install -n isoseq -c bioconda isoseq3
conda install -n isoseq -c bioconda lima
conda install -n isoseq -c bioconda pbccs

## below are optional
conda install -n isoseq -c bioconda pbcoretools # for manipulating PacBio datasets
conda install -n isoseq -c bioconda bamtools    # for converting BAM to fasta
conda install -n isoseq -c bioconda pysam       # for making CSV reports


