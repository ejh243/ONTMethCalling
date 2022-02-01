## Script to summarise and plot coverage across targetted regions
## nb guides were designed in hg19 
## ONT aligned to hg38

convertBamToGranges<-function(list){
	GRanges(seqnames = list$rname, strand = list$strand, ranges = IRanges(start = list$pos, width = list$qwidth))
}


library(Gviz)
library(vioplot)
library(Rsamtools)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# need a biostrings object for reference
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

## define target regions
allGuides<-read.csv("../Resources/SmokingEWAS/GuideRNAsFINAL.csv")

targetRegions<-cbind(aggregate(allGuides$hg38, by = list(allGuides$Chr), min), aggregate(allGuides$hg38, by = list(allGuides$Chr), max)$x)
	
targetGRanges<-GRanges(seqnames = paste0("chr", targetRegions[,1]), strand = "*", ranges = IRanges(start = targetRegions[,2], end = targetRegions[,3]))

## load chr stats
samStats1<-read.table("../Aligned/ONT/FAN62113/FAN62113.sorted.filtered.primary.stats")
samStats2<-read.table("../Aligned/ONT/FAO06764/FAO06764.sorted.filtered.primary.stats")

pdf("Plots/BarplotChromosomeDistribution.pdf", width = 10, height = 6)

barplot(samStats1[1:22,"V3"]+samStats2[1:22,"V3"], xlab = "chr", ylab = "Number of Reads", names = c(1:22))

## calculate expected number of reads based on chr sizes
totReads<-sum(samStats1[1:22,"V3"]+samStats2[1:22,"V3"])
chrProps<-samStats1[1:22,"V2"]/sum(samStats1[1:22,"V2"])
barplot(rbind(chrProps*totReads,samStats1[1:22,"V3"]+samStats2[1:22,"V3"]), beside = TRUE, xlab = "chr", ylab = "Number of Reads", names = c(1:22), col = c("grey", "darkred"))
legend("topright", c("Expected", "Actual"), col = c("grey", "darkred"), pch = 15)
dev.off()

setwd("../Aligned/ONT")
## load coverage stats
bedg1<-read.table("SumStats/FAN62113.bedg")
bedg1<-GRanges(seqnames = bedg1$V1, strand = "*", ranges = IRanges(start = bedg1$V2, end = bedg1$V3), Depth = bedg1$V4)

bedg2<-read.table("SumStats/FAO06764.bedg")
bedg2<-GRanges(seqnames = bedg2$V1, strand = "*", ranges = IRanges(start = bedg2$V2, end = bedg2$V3), Depth = bedg2$V4)

## load bam files
what <- c("rname", "strand", "pos", "qwidth")

param <- ScanBamParam(which = targetGRanges, what = what)
bam1 <- scanBam("FAN62113/FAN62113.sorted.filtered.primary.bam", param=param)
bam2 <- scanBam("FAO06764/FAO06764.sorted.filtered.primary.bam", param=param)

bamGranges1 <- lapply(bam1, convertBamToGranges)
names(bamGranges1) <- names(bam1)
bamGranges2 <- lapply(bam2, convertBamToGranges)
names(bamGranges2) <- names(bam2)

## count number of reads in each target regions
unlist(lapply(bamGranges1, length))+ unlist(lapply(bamGranges2, length))

## create alignment tracks for plotting
alTrack1 <- AlignmentsTrack("FAN62113/FAN62113.sorted.filtered.primary.bam",
  isPaired = FALSE, genome="hg38", fill = "#009E73", col = "#009E73")
alTrack2 <- AlignmentsTrack("FAO06764/FAO06764.sorted.filtered.primary.bam",
  isPaired = FALSE, genome="hg38", fill = "#009E73", col = "#009E73")
  



## plot coverage within target regions
axTrack <- GenomeAxisTrack()

for(i in 1:length(targetGRanges)){
	chr<-as.character(seqnames(targetGRanges)[i])
	start<-start(targetGRanges)[i]
	stop<-end(targetGRanges)[i]
	## expand window slightly for better resolution
	windowSize<-stop-start
	start<-start-windowSize*0.1
	stop<-stop+windowSize*0.1

	idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)

	refGenes <- UcscTrack(genome="hg38", chromosome=chr, table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=stop, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
							  symbol="name2", transcript="name", strand="strand", fill="#56B4E9",
							  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")					  
	cpgIslands <- UcscTrack(genome="hg38", chromosome=chr, track="cpgIslandExt", from=start, to=stop,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#009E73", col = "#009E73", name="CpG Islands")

	subGuides<-allGuides[which(allGuides$Chr == gsub("chr", "",chr)),]
	arrowSize<-windowSize*0.02
	subGuides$Start<-subGuides$hg38
	subGuides$End<-subGuides$hg38
	subGuides$Start[subGuides$Strand == "-1"]<-subGuides$Start[subGuides$Strand == "-1"]-arrowSize
	subGuides$End[subGuides$Strand == "1"]<-subGuides$Start[subGuides$Strand == "1"]+arrowSize

	subGuides<-GRanges(seqnames = chr, strand = subGuides$Strand, ranges = IRanges(start = subGuides$Start, end = subGuides$End))

	aTrack<-AnnotationTrack(subGuides, genome = "hg38", name = "Guides", fill = "#999999", col = "#999999")

	alTrack <- AlignmentsTrack(c(bamGranges1[[i]], bamGranges2[[i]]),
		isPaired = FALSE, genome="hg38", fill = "#E69F00", col = "#E69F00", coverageHeight = 0.2, max.height = 0.01, name = "ONT reads")
		
	pdf(paste0("Plots/CoverageHistogram_", chr, "_", start, "-", stop, ".pdf"), width = 10, height = 6)
	plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, aTrack, alTrack), from=start, to=stop, showTitle=TRUE, sizes = c(1,1.2,2,0.5,0.5,8))
	dev.off()
}

## count reads in regions
inTargets1<-subsetByOverlaps(bedg1, targetGRanges)
inTargets2<-subsetByOverlaps(bedg2, targetGRanges)

## calculate mean coverage within target regions
sum(width(inTargets1)*inTargets1$Depth)/sum(width(inTargets1))
nReads1<-aggregate(width(inTargets1)*inTargets1$Depth, list(as.character(seqnames(inTargets1))), sum)$x
nReads2<-aggregate(width(inTargets2)*inTargets2$Depth, list(as.character(seqnames(inTargets2))), sum)$x
nBases1<-aggregate(width(inTargets1), list(as.character(seqnames(inTargets1))), sum)$x
nBases2<-aggregate(width(inTargets2), list(as.character(seqnames(inTargets2))), sum)$x

((nReads1/nBases1)+(nReads2/nBases2))/2

for(i in 1:length(targetGRanges)){
	sum(width(inTargets1)*inTargets1$Depth)
}


## calculate mean coverage genome-wide
mu<-(sum(width(bedg1)*bedg1$Depth)/sum(width(bedg1)) +sum(width(bedg2)*bedg2$Depth)/sum(width(bedg2)))/2
## calculate the sd genome-wide

sqrt(sum(c(((bedg1$Depth-mu)^2)*width(bedg1), ((bedg2$Depth-mu)^2)*width(bedg2)))/sum(c(width(bedg1), width(bedg2))))

## plot distribution of read lengths within each region
pdf("Plots/HistogramReadLengthsPerTarget.pdf", height = 4, width = 12)
par(mfrow = c(1,3))
for(i in 1:length(targetGRanges)){
	hist(c(width(bamGranges1[[i]]),width(bamGranges2[[i]])), main = paste0(as.character(seqnames(targetGRanges)[i]), ":", start(targetGRanges)[i], "-", end(targetGRanges)[i]), xlab = "Read length(bp)", ylab = "Number of reads", breaks = 20) 
}
dev.off()

for(i in 1:length(targetGRanges)){
	print(mean(c(width(bamGranges1[[i]]),width(bamGranges2[[i]]))))
 }
 
 for(i in 1:length(targetGRanges)){
	print(max(c(width(bamGranges1[[i]]),width(bamGranges2[[i]]))))
 }
 
 for(i in 1:length(targetGRanges)){
	print(sum(c(width(bamGranges1[[i]]),width(bamGranges2[[i]])) > 10000))
 }
 
 lengths(targetGRanges)
# as a proportion of target region

pdf("Plots/HistogramReadLengthsPerTargetProportionOfTarget.pdf", height = 4, width = 12)
par(mfrow = c(1,3))
for(i in 1:length(targetGRanges)){
	hist(c(width(bamGranges1[[i]]),width(bamGranges2[[i]]))/lengths(targetGRanges)[i], main = paste0(as.character(seqnames(targetGRanges)[i]), ":", start(targetGRanges)[i], "-", end(targetGRanges)[i]), xlab = "Proportion of Target Region", ylab = "Number of reads", breaks = 20) 
}
dev.off()

 for(i in 1:length(targetGRanges)){
	print(mean(c(width(bamGranges1[[i]]),width(bamGranges2[[i]]))/lengths(targetGRanges)[i]))
 }
 
 for(i in 1:length(targetGRanges)){
	print(max(c(width(bamGranges1[[i]]),width(bamGranges2[[i]]))))
 }
 
 for(i in 1:length(targetGRanges)){
	print(sum(c(width(bamGranges1[[i]]),width(bamGranges2[[i]])) > 10000))
 }
 
 