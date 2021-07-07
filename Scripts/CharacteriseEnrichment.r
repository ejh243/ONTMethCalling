## characterise efficiency of crispr enrichment
## nb co-ords are in hg38

## define target regions from most extreme guides 
targetGRanges<-GRanges(seqnames = paste0("chr", targetRegions[,1]), strand = "*", ranges = IRanges(start = targetRegions[,2], end = targetRegions[,3]))

## load aligned read details
rl.files<-list.files("SumStats/", pattern = "_RL.txt")
rl.stats<-lapply(paste0("SumStats/", rl.stats), read.table)

##convert into GRanges

convertGRange<-function(data){
	return(GRanges(seqnames = data$V3, strand = "*", ranges = IRanges(start = data$V4, width = data$V6), name = data$V1, length = data$V6))
}
alignedReads<-lapply(rl.stats, convertGRange)

## for increasing windows away from target region count number of reads and mean lengths
windows<-seq(0,100000,1000)
nIntersect<-matrix(data = NA, ncol = 9, nrow = length(windows))
nIntersect[,1]<-windows
colnames(nIntersect)<-c("Distance", "S1:GFI1", "S1:CHR2I", "S1:AHRR",  "S2:GFI1", "S2:CHR2I", "S2:AHRR",  "S1:Outside", "S2:Outside")

medianRL<-matrix(data = NA, ncol = 9, nrow = length(windows))
medianRL[,1]<-windows
colnames(medianRL)<-c("Distance", "S1:GFI1", "S1:CHR2I", "S1:AHRR",  "S2:GFI1", "S2:CHR2I", "S2:AHRR",  "S1:Outside", "S2:Outside")

maxRL<-matrix(data = NA, ncol = 9, nrow = length(windows))
maxRL[,1]<-windows
colnames(maxRL)<-c("Distance", "S1:GFI1", "S1:CHR2I", "S1:AHRR", "S2:GFI1", "S2:CHR2I", "S2:AHRR",  "S1:Outside", "S2:Outside")

nLarge<-matrix(data = NA, ncol = 9, nrow = length(windows))
nLarge[,1]<-windows
colnames(nLarge)<-c("Distance", "S1:GFI1", "S1:CHR2I", "S1:AHRR",  "S2:GFI1", "S2:CHR2I", "S2:AHRR",  "S1:Outside", "S2:Outside")

## define size of large reads
large<-10000

totReads<-unlist(lapply(alignedReads, length))
rowNum<-1
for(step in seq(0,100000,1000)){
	## adjust target windows
	targetGRanges.adj<-GRanges(seqnames = paste0("chr", targetRegions[,1]), strand = "*", ranges = IRanges(start = targetRegions[,2]-step, end = targetRegions[,3]+step))

	subsetReads<-lapply(alignedReads, findOverlaps, targetGRanges.adj)
	## count
	nIntersect[rowNum,2:7]<-unlist(lapply(lapply(subsetReads, subjectHits), table))
	nIntersect[rowNum,8:9]<-totReads-unlist(lapply(lapply(subsetReads, subjectHits), length))

	## pull out readlength stats
	collateMed<-rep(NA, 4*length(alignedReads))
	collateMax<-rep(NA, 4*length(alignedReads))
	collateLarge<-rep(NA, 4*length(alignedReads))
	for(i in 1:length(alignedReads)){
		RLall<-mcols(alignedReads[[i]])[lapply(subsetReads, queryHits)[[i]],]$length
		whichRegion<-lapply(subsetReads, subjectHits)[[i]]
		collateMed[(1:3)+((i-1)*4)]<-aggregate(RLall, by = list(whichRegion), FUN = median)$x
		collateMax[(1:3)+((i-1)*4)]<-aggregate(RLall, by = list(whichRegion), FUN = max)$x
		collateLarge[(1:3)+((i-1)*4)]<-aggregate(RLall > large, by = list(whichRegion), FUN = sum)$x
		RLout<-mcols(alignedReads[[i]])[-lapply(subsetReads, queryHits)[[i]],]$length
		collateMed[4+((i-1)*4)]<-median(RLout)
		collateMax[4+((i-1)*4)]<-max(RLout)
		collateLarge[4+((i-1)*4)]<-sum(RLout > large)
	}
	medianRL[rowNum, -1]<-collateMed
	maxRL[rowNum, -1]<-collateMax
	
	rowNum<-rowNum+1
}



## create plots
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## convert to % for barplot
nIntersect<-nIntersect[,c(1:4,8,5:7,9)]
nIntersect.per<-nIntersect
nIntersect.per[,2:5]<-nIntersect.per[,2:5]/rowSums(nIntersect[,c(2:5)])*100
nIntersect.per[,6:9]<-nIntersect.per[,6:9]/rowSums(nIntersect[,c(6:9)])*100

pdf("Plots/BarplotPercentageInRegion.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
barplot(t(nIntersect.per[,2:5]), col = colorBlindBlack8[1:4], ylab = "% of aligned reads", main = "Sample 1")
barplot(t(nIntersect.per[,6:9]), col = colorBlindBlack8[1:4], ylab = "% of aligned reads", main = "Sample 2", legend.text = c("GFI1", "Chr2_I", "AHRR", "Outside"))
dev.off()

pdf("Plots/LineGraphPercentageInRegionAgainstDistance.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(nIntersect.per[,1], nIntersect.per[,2], type = "l", xlab = "Distance from target", ylab = "% Aligned Reads", lwd = 2, col = colorBlindBlack8[1], main = "Sample 1", ylim = c(0,2.5))
lines(nIntersect.per[,1], nIntersect.per[,3], type = "l", lwd = 2, col = colorBlindBlack8[2])
lines(nIntersect.per[,1], nIntersect.per[,4], type = "l", lwd = 2, col = colorBlindBlack8[3])
lines(nIntersect.per[,1], rowSums(nIntersect.per[,2:4]), type = "l", lwd = 2, col = colorBlindBlack8[7])
legend("topleft", c("GFI1", "Chr2I", "AHRR", "All"), col = colorBlindBlack8[c(1:3,7)], lwd = 2)

plot(nIntersect.per[,1], nIntersect.per[,6], type = "l", xlab = "Distance from target", ylab = "% Aligned Reads", lwd = 2, col = colorBlindBlack8[1], main = "Sample 2", ylim = c(0,2.5))
lines(nIntersect.per[,1], nIntersect.per[,7], type = "l", lwd = 2, col = colorBlindBlack8[2])
lines(nIntersect.per[,1], nIntersect.per[,8], type = "l", lwd = 2, col = colorBlindBlack8[3])
lines(nIntersect.per[,1], rowSums(nIntersect.per[,6:8]), type = "l", lwd = 2, col = colorBlindBlack8[7])
dev.off()

## plot RL stats against distance from target
pdf("Plots/LineGraphMedianRLAgainstDistance.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(medianRL[,1], medianRL[,2], type = "l", xlab = "Distance from target", ylab = "Median read length", lwd = 2, col = colorBlindBlack8[1], main = "Sample 1", ylim = range(medianRL[,-1]))
lines(medianRL[,1], medianRL[,3], type = "l", lwd = 2, col = colorBlindBlack8[2])
lines(medianRL[,1], medianRL[,4], type = "l", lwd = 2, col = colorBlindBlack8[3])
lines(medianRL[,1], medianRL[,8], type = "l", lwd = 2, col = colorBlindBlack8[4])


plot(medianRL[,1], medianRL[,5], type = "l", xlab = "Distance from target", ylab = "Median read length", lwd = 2, col = colorBlindBlack8[1], main = "Sample 2", ylim = range(medianRL[,-1]))
lines(medianRL[,1], medianRL[,6], type = "l", lwd = 2, col = colorBlindBlack8[2])
lines(medianRL[,1], medianRL[,7], type = "l", lwd = 2, col = colorBlindBlack8[3])
lines(medianRL[,1], medianRL[,9], type = "l", lwd = 2, col = colorBlindBlack8[4])
legend("topleft", c("GFI1", "Chr2I", "AHRR", "Outside"), col = colorBlindBlack8[c(1:4)], lwd = 2)
dev.off()

pdf("Plots/LineGraphMaxRLAgainstDistance.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(maxRL[,1], maxRL[,2], type = "l", xlab = "Distance from target", ylab = "Median read length", lwd = 2, col = colorBlindBlack8[1], main = "Sample 1", ylim = range(maxRL[,-1]))
lines(maxRL[,1], maxRL[,3], type = "l", lwd = 2, col = colorBlindBlack8[2])
lines(maxRL[,1], maxRL[,4], type = "l", lwd = 2, col = colorBlindBlack8[3])
lines(maxRL[,1], maxRL[,8], type = "l", lwd = 2, col = colorBlindBlack8[4])


plot(maxRL[,1], maxRL[,5], type = "l", xlab = "Distance from target", ylab = "Maximum read length", lwd = 2, col = colorBlindBlack8[1], main = "Sample 2", ylim = range(maxRL[,-1]))
lines(maxRL[,1], maxRL[,6], type = "l", lwd = 2, col = colorBlindBlack8[2])
lines(maxRL[,1], maxRL[,7], type = "l", lwd = 2, col = colorBlindBlack8[3])
lines(maxRL[,1], maxRL[,9], type = "l", lwd = 2, col = colorBlindBlack8[4])
legend("topleft", c("GFI1", "Chr2I", "AHRR", "Outside"), col = colorBlindBlack8[c(1:4)], lwd = 2)
dev.off()

### boxplots of read lengths within target regions + 10000
step<-10000
targetGRanges.adj<-GRanges(seqnames = paste0("chr", targetRegions[,1]), strand = "*", ranges = IRanges(start = targetRegions[,2]-step, end = targetRegions[,3]+step))

subsetReads<-lapply(alignedReads, findOverlaps, targetGRanges.adj)

pdf("Plots/BoxplotRLReadsinTargets.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
## pull out readlength stats
for(i in 1:length(alignedReads)){
	RLall<-mcols(alignedReads[[i]])[lapply(subsetReads, queryHits)[[i]],]$length
	whichRegion<-lapply(subsetReads, subjectHits)[[i]]
	RLout<-mcols(alignedReads[[i]])[-lapply(subsetReads, queryHits)[[i]],]$length
	boxplot(RLall ~ whichRegion, col = colorBlindBlack8[c(1:3)], xlim = c(0.5,4.5), ylab = "Read Length", xlab = "", names = c("GFI1", "Chr2_I", "AHRR"))
	boxplot(RLout, col = colorBlindBlack8[4], add = TRUE, at = 4, names = "Outside")
}
dev.off()
	

## for each guide count number of reads within X bp.
allGuides<-read.csv("../../Resources/SmokingEWAS/GuideRNAsFINAL.csv")
windows<-c(1,seq(50,10000,50))
nIntersect<-matrix(data = NA, ncol = 1+2*nrow(allGuides), nrow = length(windows))
nIntersect[,1]<-windows

rowNum<-1
for(step in windows){
	## adjust target windows
	guideGRanges.adj<-GRanges(seqnames = paste0("chr", allGuides$Chr), strand = "*", ranges = IRanges(start = allGuides$hg38-step, end = allGuides$hg38+step))

	subsetReads<-lapply(alignedReads, findOverlaps, guideGRanges.adj)
	## count
	counts<-unlist(lapply(lapply(subsetReads, subjectHits), table))
	while(length(counts) < nrow(allGuides)*2){
		insert<-min(which(names(counts) != rep(c(1:nrow(allGuides)),2)))
		counts<-c(counts[1:(insert-1)], NA, counts[insert:length(counts)])
	}
	names(counts) == rep(c(1:nrow(allGuides)),2)
	nIntersect[rowNum,-1]<-counts
	rowNum<-rowNum+1
}

pdf("Plots/LineGraphReadsPerGuideAgainstDistance.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(nIntersect[,1], nIntersect[,2], xlab = "Distance from cut", ylab = "Number of reads", type = "l", lwd = 2, ylim = range(nIntersect[,-1], na.rm = TRUE), main = "Sample 1", col = colorBlindBlack8[as.numeric(as.factor(allGuides$Chr))][1])
for(i in 3:(nrow(allGuides)+1)){
	lines(nIntersect[,1], nIntersect[,i], lwd = 2, col = colorBlindBlack8[as.numeric(as.factor(allGuides$Chr))][i])
}

plot(nIntersect[,1], nIntersect[,20], xlab = "Distance from cut", ylab = "Number of reads", type = "l", lwd = 2, ylim = range(nIntersect[,-1], na.rm = TRUE), main = "Sample 2", col = colorBlindBlack8[as.numeric(as.factor(allGuides$Chr))][1])
for(i in 3:(nrow(allGuides)+1)){
	lines(nIntersect[,1], nIntersect[,i+18], lwd = 2, col = colorBlindBlack8[as.numeric(as.factor(allGuides$Chr))][i])
}
dev.off()

pdf("Plots/ScatterplotCompareNumberGuidesPerSample.pdf", width = 5, height = 5)
plot(nIntersect[21,2:19], nIntersect[21,20:37], pch = 16, col = colorBlindBlack8[as.numeric(as.factor(allGuides$Chr))], xlab = "Sample 1", ylab = "Sample 2")
abline(a=0, b = 1, lty = 2)
dev.off()

step<-1000
	guideGRanges.adj<-GRanges(seqnames = paste0("chr", allGuides$Chr), strand = "*", ranges = IRanges(start = allGuides$hg38-step, end = allGuides$hg38+step))

	subsetReads<-lapply(alignedReads, findOverlaps, guideGRanges.adj)
pdf("Plots/BoxplotReadLengthPerGuide.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
for(i in 1:length(alignedReads)){
   RLall<-mcols(alignedReads[[i]])[lapply(subsetReads, queryHits)[[i]],]$length
   whichRegion<-lapply(subsetReads, subjectHits)[[i]]
	boxplot(RLall ~ whichRegion, col = colorBlindBlack8[as.numeric(as.factor(allGuides$Chr))])
}
dev.off()

### plot against guide location

library(Gviz)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# need a biostrings object for reference
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

axTrack <- GenomeAxisTrack()

dTrack1 <- DataTrack(start = allGuides$hg38, width = 1000, genome = "hg38", type = "histogram", 
                     chromosome = allGuides$Chr, name = "Sample 1", data = nIntersect[21,2:19])
dTrack2 <- DataTrack(start = allGuides$hg38, width = 1000, genome = "hg38", type = "histogram", 
                     chromosome = allGuides$Chr, name = "Sample 1", data = nIntersect[21,20:37])

## these coords have been lifted over
chr<-5
start<-304177-10000
stop<-438290+10000

idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)

refGenes <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=stop, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
						  symbol="name2", transcript="name", strand="strand", fill="#8282d2",
						  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")					  
cpgIslands <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), track="cpgIslandExt", from=start, to=stop,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#006400", name="CpG Islands")

					 					 

pdf("Plots/NumberReadsByGuide_AHRR.pdf", width = 10, height = 4)
plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, dTrack1, dTrack2), from=start, to=stop, showTitle=TRUE)
dev.off()



chr<-2
start<-232418619-10000
stop<-232420224+10000


idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)

refGenes <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=stop, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
						  symbol="name2", transcript="name", strand="strand", fill="#8282d2",
						  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")					  
cpgIslands <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), track="cpgIslandExt", from=start, to=stop,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#006400", name="CpG Islands")

				 					 

pdf("Plots/NumberReadsByGuide_GFI1.pdf", width = 10, height = 4)
plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, dTrack1, dTrack2), from=start, to=stop, showTitle=TRUE)
dev.off()


chr<-1
start<-92474761-10000
stop<-92486876+10000

idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)

refGenes <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=stop, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
						  symbol="name2", transcript="name", strand="strand", fill="#8282d2",
						  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")					  
cpgIslands <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), track="cpgIslandExt", from=start, to=stop,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#006400", name="CpG Islands")

				 					 

pdf("Plots/NumberReadsByGuide_Chr2_I.pdf", width = 10, height = 4)
plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, dTrack1, dTrack2), from=start, to=stop, showTitle=TRUE)
dev.off()


## compare number of reads to IDT scores
pdf("Plots/ScatterplotNumberReadsIDTScores.pdf", width = 10, height = 4)
par(mfrow = c(1,2))
plot(nIntersect[21,2:19], allGuides$Specificity.Score, pch = 16, xlab = "Number of Reads", ylab = "Specificity score", xlim = range(nIntersect[21,-1], na.rm = TRUE))
points(nIntersect[21,20:37], allGuides$Specificity.Score, pch = 16, col = "red")

plot(nIntersect[21,2:19], allGuides$Efficiency.Score, pch = 16, xlab = "Number of Reads", ylab = "Efficiency score", xlim = range(nIntersect[21,-1], na.rm = TRUE))
points(nIntersect[21,20:37], allGuides$Efficiency.Score, pch = 16, col = "red")
dev.off()
