## run nanopolish to call dnam
## nanopolish calls dnam at a group/regional level

### index fast5 files
FQFILES=($(ls ${DATADIR}/*/fastq_pass/*merged.fastq.gz))
for f in ${FQFILES[@]};	
  do
    basename=$(basename $f)
	samplename=${basename/_merged.fastq.gz}
	${NANOPOLISH} index -d ${DATADIR}/*${samplename}*/fast5_pass/ ${f}
done
	

## needs aligned, indexed bam file
cd ${ALIGNEDDIR}/
BAMFILES=($(ls */*.sorted.filtered.primary.bam))

echo "Number of bam files found for dnam calling:"" ""${#BAMFILES[@]}""" 
 
mkdir -p ${OUTDIR}/Nanopolish
for bf in ${BAMFILES[@]}
do
    samplename=$(dirname $bf)
	
	mkdir -p ${OUTDIR}/Nanopolish/${samplename}
	
	## find fq files for this samplename
	fqf=${DATADIR}/*${samplename}*/fastq_pass/${samplename}_merged.fastq.gz
	
	${NANOPOLISH} call-methylation -t 16 -r ${fqf} -b ${bf} -g ${REFGENOME} > ${OUTDIR}/Nanopolish/${samplename}/${samplename}_methylation_calls.tsv
	
	${NANOPOLISHSCRIPTS}/calculate_methylation_frequency.py ${OUTDIR}/Nanopolish/${samplename}/${samplename}_methylation_calls.tsv > ${OUTDIR}/Nanopolish/${samplename}/${samplename}_methylation_frequency.tsv
done

