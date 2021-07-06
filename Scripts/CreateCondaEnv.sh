## set up ONT conda environment

module load Anaconda3/5.2.0

conda create --name ont_meth

conda install -n ont_meth -c bioconda minimap2
conda install -n ont_meth -c bioconda samtools
conda install -n ont_meth -c bioconda nanopolish