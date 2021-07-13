## compare DNAm values from ONT to Array values
setwd("DNAmethylation/")

files<-list.files("Nanopolish/", pattern = "_methylation_frequency.tsv", recursive=TRUE)

ont.nano<-lapply(paste0("Nanopolish/",files), read.table, header = TRUE)

## plot DNAm distribution and coverage

pdf("Plots/HistogramDNAmAll.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
hist(ont.nano[[1]]$methylated_frequency, xlab = "% DNA methylation", ylab = "Number of Regions", main = "Sample 1")
hist(ont.nano[[2]]$methylated_frequency, xlab = "% DNA methylation", ylab = "Number of Regions", main = "Sample 2")
dev.off()

pdf("Plots/HistogramReadDepthAll.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
hist(ont.nano[[1]]$called_sites, xlab = "Number of Reads*CpGs", ylab = "Number of Regions", main = "Sample 1", breaks = seq(0,660, 10))
hist(ont.nano[[2]]$called_sites, xlab = "Number of Reads*CpGs", ylab = "Number of Regions", main = "Sample 2", breaks = seq(0,660, 10))

par(mfrow = c(1,2))
hist(ont.nano[[1]]$called_sites, xlab = "Number of Reads*CpGs", ylab = "Number of Regions", main = "Sample 1", breaks = seq(0,660, 1), xlim = c(0,10))
hist(ont.nano[[2]]$called_sites, xlab = "Number of Reads*CpGs", ylab = "Number of Regions", main = "Sample 2", breaks = seq(0,660, 1), xlim = c(0,10))


par(mfrow = c(1,2))
hist(ont.nano[[1]]$called_sites/ont.nano[[1]]$num_motifs_in_group, xlab = "Number of Reads", ylab = "Number of Regions", main = "Sample 1", breaks = seq(0,660, 1), xlim = c(0,10))
hist(ont.nano[[2]]$called_sites/ont.nano[[2]]$num_motifs_in_group, xlab = "Number of Reads", ylab = "Number of Regions", main = "Sample 2", breaks = seq(0,660, 1), xlim = c(0,10))
dev.off()


## to speed up comparison, start by removing anything with less than 5 reads.
filterRD<-function(data, threshold){
	return(data[which(data$called_sites > threshold),])
}

ont.nano.filt<-lapply(ont.nano, filterRD, 5)


pdf("Plots/HistogramDNAmAllRD5.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
hist(ont.nano.filt[[1]]$methylated_frequency, xlab = "% DNA methylation", ylab = "Number of Regions", main = "Sample 1")
hist(ont.nano.filt[[2]]$methylated_frequency, xlab = "% DNA methylation", ylab = "Number of Regions", main = "Sample 2")
dev.off()

## convert to GRanges
sam1<-GRanges(seqnames = ont.nano.filt[[1]]$chromosome, strand = "*", ranges = IRanges(start = ont.nano.filt[[1]]$start, end = ont.nano.filt[[1]]$end), DNAm = ont.nano.filt[[1]]$methylated_frequency, Reads = ont.nano.filt[[1]]$called_sites, Coverage = ont.nano.filt[[1]]$called_sites/ont.nano.filt[[1]]$num_motifs_in_group)

sam2<-GRanges(seqnames = ont.nano.filt[[2]]$chromosome, strand = "*", ranges = IRanges(start = ont.nano.filt[[2]]$start, end = ont.nano.filt[[2]]$end), DNAm = ont.nano.filt[[2]]$methylated_frequency, Reads = ont.nano.filt[[2]]$called_sites, Coverage = ont.nano.filt[[2]]$called_sites/ont.nano.filt[[2]]$num_motifs_in_group)

### how many overlapping sites GW?
sharedSites<-findOverlaps(sam1, sam2)

pdf("Plots/ScatterplotDNAmGenomewide.pdf")
plot(sam1$DNAm[queryHits(sharedSites)], sam2$DNAm[subjectHits(sharedSites)], pch = 16, xlab = "Sample 1", ylab = "Sample 2", main = "CpGs with > 5 reads")
dev.off()



## focus in on target regions

allGuides<-read.csv("../Resources/SmokingEWAS/GuideRNAsFINAL.csv")

targetRegions<-cbind(aggregate(allGuides$hg38, by = list(allGuides$Chr), min), aggregate(allGuides$hg38, by = list(allGuides$Chr), max)$x)
	
targetGRanges<-GRanges(seqnames = paste0("chr", targetRegions[,1]), strand = "*", ranges = IRanges(start = targetRegions[,2], end = targetRegions[,3]))

inTargets1<-subsetByOverlaps(sam1, targetGRanges)
inTargets2<-subsetByOverlaps(sam2, targetGRanges)

### how many overlapping sites ?
sharedSites<-findOverlaps(inTargets1, inTargets2)

pdf("Plots/ScatterplotDNAmInCRISPRRegions.pdf")
plot(inTargets1$DNAm[queryHits(sharedSites)], inTargets2$DNAm[subjectHits(sharedSites)], pch = 16, xlab = "Sample 1", ylab = "Sample 2", main = "CpGs with > 5 reads")
dev.off()

## first genome-wide correlations and error metrics
load("../Resources/SmokingEWAS/ArrayData.rda")

probeAnno<-read.table("/gpfs/mrc0/projects/Research_Project-MRC190311/References/EPICArray/EPIC.anno.GRCh38.tsv", header = TRUE, fill = TRUE)
probeAnno<-probeAnno[match(rownames(smokebetas), probeAnno$probeID),]
arrayData<-GRanges(seqnames = probeAnno$chrm, strand = "*", ranges = IRanges(start = probeAnno$start, end = probeAnno$end), smokebetas)

overlapArray1<-findOverlaps(sam1, arrayData)
overlapArray2<-findOverlaps(sam2, arrayData)

## need to double check which sample is which.
pdf("Plots/ScatterplotDNAmArrayvsNanopolishGenomewide.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(sam1$DNAm[queryHits(overlapArray1)], arrayData$"Non-Smoker"[subjectHits(overlapArray1)], pch = 16, xlab = "ONT", ylab = "EPIC array")
mtext(side = 3, line = 0, adj = 1, paste0("cor = ", signif(cor(sam1$DNAm[queryHits(overlapArray1)], arrayData$"Non-Smoker"[subjectHits(overlapArray1)]),3), "; n = ", length(overlapArray1), " sites"))

plot(sam2$DNAm[queryHits(overlapArray2)], arrayData$Smoker[subjectHits(overlapArray2)], pch = 16, xlab = "ONT", ylab = "EPIC array")
mtext(side = 3, line = 0, adj = 1, paste0("cor = ", signif(cor(sam2$DNAm[queryHits(overlapArray2)], arrayData$Smoker[subjectHits(overlapArray2)]),3), "; n = ", length(overlapArray2), " sites"))
dev.off()

## cor as a function of read depth
thresholds<-seq(5,100,1)
errorMetrics<-matrix(data = NA, nrow = length(thresholds), ncol = 7)
errorMetrics[,1]<-thresholds
rowNum<-1
for(rdThres in thresholds){
	keep<-which(sam2$Coverage[queryHits(overlapArray2)] > rdThres)
	r1<-cor(sam2$DNAm[queryHits(overlapArray2)[keep]], arrayData$Smoker[subjectHits(overlapArray2)[keep]])
	rmse1<-sqrt(median((sam2$DNAm[queryHits(overlapArray2)[keep]] - arrayData$Smoker[subjectHits(overlapArray2)[keep]])^2))
	ncpg1<-length(keep)
	
	keep<-which(sam1$Coverage[queryHits(overlapArray1)] > rdThres)
	r2<-cor(sam1$DNAm[queryHits(overlapArray1)[keep]], arrayData$"Non-Smoker"[subjectHits(overlapArray1)[keep]])
	rmse2<-sqrt(median((sam1$DNAm[queryHits(overlapArray1)[keep]] - arrayData$"Non-Smoker"[subjectHits(overlapArray1)[keep]])^2))
	ncpg2<-length(keep)
	
	errorMetrics[rowNum,2:7]<-c(r1, rmse1, ncpg1, r2, rmse2, ncpg2)
	rowNum<-rowNum+1
}

pdf("Plots/LineGraphArrayONTComparisionAgainstReadDepth.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(errorMetrics[,1], errorMetrics[,2], xlab = "RD Threshold", ylab = "r", type = "l", lwd = 2, ylim = c(0.8,1), xlim = c(0,45))
lines(errorMetrics[,1], errorMetrics[,5], lwd = 2)

plot(errorMetrics[,1], errorMetrics[,3], xlab = "RD Threshold", ylab = "RMSE", type = "l", lwd = 2, ylim = c(0,0.18), xlim = c(0,45))
lines(errorMetrics[,1], errorMetrics[,6], lwd = 2)
dev.off()

## as a function of DNAm level

errorMetrics.dnam<-matrix(data = NA, ncol = 7, nrow = 10)
errorMetrics.dnam[,1]<-seq(0.05, 0.95, 0.1)

dnamFilter<-cut(arrayData$Smoker[subjectHits(overlapArray1)], breaks = seq(0,1,0.1))
rowNum<-1
for(each in levels(dnamFilter)){
	keep<-which(dnamFilter == each)
	r1<-cor(sam2$DNAm[queryHits(overlapArray2)[keep]], arrayData$Smoker[subjectHits(overlapArray2)[keep]])
	rmse1<-sqrt(median((sam2$DNAm[queryHits(overlapArray2)[keep]] - arrayData$Smoker[subjectHits(overlapArray2)[keep]])^2))
	ncpg1<-length(keep)
	errorMetrics.dnam[rowNum,2:4]<-c(r1, rmse1,ncpg1)
	
	rowNum<-rowNum+1
}

dnamFilter<-cut(arrayData$"Non-Smoker"[subjectHits(overlapArray2)], breaks = seq(0,1,0.1))
rowNum<-1
for(each in levels(dnamFilter)){	
	keep<-which(dnamFilter == each)
	r2<-cor(sam1$DNAm[queryHits(overlapArray1)[keep]], arrayData$"Non-Smoker"[subjectHits(overlapArray1)[keep]])
	rmse2<-sqrt(median((sam1$DNAm[queryHits(overlapArray1)[keep]] - arrayData$"Non-Smoker"[subjectHits(overlapArray1)[keep]])^2))
	ncpg2<-length(keep)
	
	errorMetrics.dnam[rowNum,5:7]<-c(r2, rmse2,ncpg2)
		rowNum<-rowNum+1
}


pdf("Plots/LineGraphArrayONTComparisionAgainstDNAmLevel.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(errorMetrics.dnam[,1], errorMetrics.dnam[,2], xlab = "DNAm mean", ylab = "r", type = "l", lwd = 2, ylim = c(0.8,1))
lines(errorMetrics.dnam[,1], errorMetrics.dnam[,5], lwd = 2)

plot(errorMetrics.dnam[,1], errorMetrics.dnam[,3], xlab = "DNAm mean", ylab = "RMSE", type = "l", lwd = 2, ylim = c(0,0.18))
lines(errorMetrics.dnam[,1], errorMetrics.dnam[,6], lwd = 2)
dev.off()




