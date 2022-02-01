## compare DNAm values from ONT to Array values

filterRD<-function(data, threshold){
	return(data[which(data$called_sites/data$num_motifs_in_group > threshold),])
}

library(GenomicRanges)

setwd("DNAmethylation/")

files<-list.files("Nanopolish/", pattern = "_methylation_frequency.tsv", recursive=TRUE)

ont.nano<-lapply(paste0("Nanopolish/",files), read.table, header = TRUE)

## only consider CpGs with at least 10 reads.
ont.nano.filt<-lapply(ont.nano, filterRD, 10)

## convert to GRanges
sam1<-GRanges(seqnames = ont.nano.filt[[1]]$chromosome, strand = "*", ranges = IRanges(start = ont.nano.filt[[1]]$start, end = ont.nano.filt[[1]]$end), DNAm = ont.nano.filt[[1]]$methylated_frequency, Reads = ont.nano.filt[[1]]$called_sites, Coverage = ont.nano.filt[[1]]$called_sites/ont.nano.filt[[1]]$num_motifs_in_group, nCpG = ont.nano.filt[[1]]$num_motifs_in_group)

sam2<-GRanges(seqnames = ont.nano.filt[[2]]$chromosome, strand = "*", ranges = IRanges(start = ont.nano.filt[[2]]$start, end = ont.nano.filt[[2]]$end), DNAm = ont.nano.filt[[2]]$methylated_frequency, Reads = ont.nano.filt[[2]]$called_sites, Coverage = ont.nano.filt[[2]]$called_sites/ont.nano.filt[[2]]$num_motifs_in_group, nCpG = ont.nano.filt[[2]]$num_motifs_in_group)

## focus in on target regions

allGuides<-read.csv("../Resources/SmokingEWAS/GuideRNAsFINAL.csv")

targetRegions<-cbind(aggregate(allGuides$hg38, by = list(allGuides$Chr), min), aggregate(allGuides$hg38, by = list(allGuides$Chr), max)$x)
	
targetGRanges<-GRanges(seqnames = paste0("chr", targetRegions[,1]), strand = "*", ranges = IRanges(start = targetRegions[,2], end = targetRegions[,3]))

inTargets1<-subsetByOverlaps(sam1, targetGRanges)
inTargets2<-subsetByOverlaps(sam2, targetGRanges)

### how many overlapping sites ?
sharedSites<-findOverlaps(inTargets1, inTargets2)
allCpGs<-union(inTargets1, inTargets2)


aggregate(inTargets1$nCpG, by = list(as.character(seqnames(inTargets1))), sum)
aggregate(inTargets2$nCpG, by = list(as.character(seqnames(inTargets2))), sum)
aggregate(inTargets1$nCpG[queryHits(sharedSites)], by = list(as.character(seqnames(inTargets1[queryHits(sharedSites)]))), sum)


## first genome-wide correlations and error metrics
load("../Resources/SmokingEWAS/ArrayData.rda")

probeAnno<-read.table("/gpfs/mrc0/projects/Research_Project-MRC190311/References/EPICArray/EPIC.anno.GRCh38.tsv", header = TRUE, fill = TRUE)
probeAnno<-probeAnno[match(rownames(smokebetas), probeAnno$probeID),]
arrayData<-GRanges(seqnames = probeAnno$chrm, strand = "*", ranges = IRanges(start = probeAnno$start, end = probeAnno$end), smokebetas)

## count number of sites on EPIC array within these regions
epicOverlap<-subsetByOverlaps(arrayData, targetGRanges)
table(seqnames(epicOverlap))

## compare the spacing profile
intraDist.ont<-NULL
intraDist.array<-NULL
for(i in 1:length(targetGRanges)){
	## first for ONT
	subCpGs<-allCpGs[which(as.character(seqnames(allCpGs)) == as.character(seqnames(targetGRanges)[i])),]
	intraDist.ont<-rbind(intraDist.ont, cbind(as.character(seqnames(targetGRanges)[i]), width(gaps(subCpGs)[-1])))
	## Need to add in distances between CpGs called in a single region
	## calculate average distance	
	regionsIndex<-which(subCpGs$nCpG > 1)
	if(length(regionsIndex)> 0){
		## add multiple times to vector of distances
		intraDist.ont<-rbind(intraDist.ont, cbind(as.character(seqnames(targetGRanges)[i]),rep(width(subCpGs[regionsIndex])/(subCpGs$nCpG -1), (subCpGs$nCpG -1))))
	}
	
	## second for EPIC array
	subCpGs<-epicOverlap[which(as.character(seqnames(epicOverlap)) == as.character(seqnames(targetGRanges)[i])),]
	intraDist.array<-rbind(intraDist.array, cbind(as.character(seqnames(targetGRanges)[i]), width(gaps(subCpGs)[-1])))

}

aggregate(as.numeric(intraDist.array[,2]), by = list(intraDist.array[,1]), median)
aggregate(as.numeric(intraDist.ont[,2]), by = list(intraDist.ont[,1]), median)

## identify sites measured in both
overlapArray1<-findOverlaps(sam1, arrayData)
overlapArray2<-findOverlaps(sam2, arrayData)

dnam.ont<-c(sam1$DNAm[queryHits(overlapArray1)],sam2$DNAm[queryHits(overlapArray2)])
dnam.epic<-c(arrayData$"Non-Smoker"[subjectHits(overlapArray1)], arrayData$Smoker[subjectHits(overlapArray2)])
ind<-c(rep(1, length(overlapArray1)), rep(2, length(overlapArray2)))

pdf("Plots/ScatterplotDNAmArrayvsNanopolishWithinTargettedRegions.pdf")
plot(dnam.ont, dnam.epic, pch = 16, xlab = "ONT", ylab = "EPIC array", col = c("#009E73", "#56B4E9")[ind])
mtext(side = 3, line = 0, adj = 1, paste0("cor = ", signif(cor(dnam.ont, dnam.epic),3)))
dev.off()

## calculate RMSE
sqrt(median((dnam.ont - dnam.epic)^2))


