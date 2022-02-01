## Summarise coverage of aligned reads, and calculate length of aligned reads for further characterisation


cd ${ALIGNEDDIR}/
## focus on primary, high quality aligned reads
BAMFILES=($(ls */*.sorted.filtered.primary.bam))


## each file separately
for f in ${BAMFILES[@]}
do
  id=$(dirname $f)
  genomeCoverageBed -ibam $f -bg > SumStats/${id}.bedg
  
  ## calculate length of each read in bam file
  samtools view $f | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' > SumStats/${id}_RL.tmp
  samtools view $f | cut -f 1-5 > SumStats/readinfo.tmp
  paste SumStats/readinfo.tmp SumStats/${id}_RL.tmp > SumStats/${id}_RL.txt
  rm SumStats/readinfo.tmp
  rm SumStats/${id}_RL.tmp

done

