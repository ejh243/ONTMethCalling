## align ONT data with minimap
  

FOLDERS=($(ls ${DATADIR}))

for sample in ${FOLDERS[@]}
do 
	FQFILES=($(ls ${DATADIR}/${sample}/fastq_pass/*_pass*.fastq.gz))
    
	id=$(basename ${FQFILES})
	id=${id/_pass*.fastq.gz}
	
	## Merge into single fastq file for each sample/flowcell
	cat ${FQFILES[@]} > ${DATADIR}/${sample}/fastq_pass/${id}_merged.fastq.gz
	  
	mkdir -p ${ALIGNEDDIR}/${id}	
	echo "Aligning"" ${f}"
	  
    minimap2 -a -t 16 -x map-ont ${REFGENOME} ${DATADIR}/${sample}/fastq_pass/${id}_merged.fastq.gz \
       > ${ALIGNEDDIR}/${id}/${id}.sam \
       2> ${ALIGNEDDIR}/${id}/${id}.sam.log
  
    ## sort aligned file
    samtools sort -T tmp -o ${ALIGNEDDIR}/${id}/${id}.sorted.bam ${ALIGNEDDIR}/${id}/${id}.sam
	
	## index aligned file
    samtools index ${ALIGNEDDIR}/${id}/${id}.sorted.bam

    ## delete sam file
    rm ${ALIGNEDDIR}/${id}/${id}.sam
  
done