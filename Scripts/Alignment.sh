## align ONT data with minimap
  
## find fastq files within multiple folders where each folder is a separate flowcell

FQFILES=($(ls ${DATADIR}/*/fastq_pass/*q.gz))
echo "Number of fastq files found for alignment:"" ""${#FQFILES[@]}"""  
 
		
for f in ${FQFILES[@]};	
  do
    basename=$(basename $f)
	folder=${basename/_pass*q.gz}
	samplename=${basename/.fastq.gz}
	  
	mkdir -p ${ALIGNEDDIR}/${folder}	
	echo "Aligning"" ${basename}"
	  
    minimap2 -a -t 16 -x map-ont ${REFGENOME} ${f} \
       > ${ALIGNEDDIR}/${folder}/${samplename}.sam \
       2> ${ALIGNEDDIR}/${folder}/${samplename}.sam.log
  
    ## sort aligned file
    samtools sort -T tmp -o ${ALIGNEDDIR}/${folder}/${samplename}.sorted.bam ${ALIGNEDDIR}/${folder}/${samplename}.sam
  
    ## delete sam file
    rm ${ALIGNEDDIR}/${folder}/${samplename}.sam
  
done

## merge aligned files within each folder where a folder represents a flowcell/sample 
cd ${ALIGNEDDIR}/
FOLDERS=($(ls ))

for sample in ${FOLDERS[@]}
do 
  BAMFILES=($(ls ${sample}/*.sorted.bam))
  
  samtools merge ${sample}/merged.sorted.bam ${BAMFILES[@]}
  
  ## index aligned file
  samtools index ${sample}/merged.sorted.bam

  ## delete intermediate bam files
  rm ${BAMFILES[@]}
done

  
