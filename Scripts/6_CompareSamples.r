## compare DNAm values from ONT to Array values


filterRD<-function(data, threshold){
	return(data[which(data$called_sites/data$num_motifs_in_group > threshold),])
}

library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPseeker)

setwd("DNAmethylation/")

files<-list.files("Nanopolish/", pattern = "_methylation_frequency.tsv", recursive=TRUE)

ont.nano<-lapply(paste0("Nanopolish/",files), read.table, header = TRUE)

## only consider CpGs with at least 10 reads.
ont.nano.filt<-lapply(ont.nano, filterRD, 10)

## convert to GRanges
sam1<-GRanges(seqnames = ont.nano.filt[[1]]$chromosome, strand = "*", ranges = IRanges(start = ont.nano.filt[[1]]$start, end = ont.nano.filt[[1]]$end), DNAm = ont.nano.filt[[1]]$methylated_frequency, Reads = ont.nano.filt[[1]]$called_sites, Coverage = ont.nano.filt[[1]]$called_sites/ont.nano.filt[[1]]$num_motifs_in_group, nCpG = ont.nano.filt[[1]]$num_motifs_in_group, MReads = ont.nano.filt[[1]]$called_sites_methylated)

sam2<-GRanges(seqnames = ont.nano.filt[[2]]$chromosome, strand = "*", ranges = IRanges(start = ont.nano.filt[[2]]$start, end = ont.nano.filt[[2]]$end), DNAm = ont.nano.filt[[2]]$methylated_frequency, Reads = ont.nano.filt[[2]]$called_sites, Coverage = ont.nano.filt[[2]]$called_sites/ont.nano.filt[[2]]$num_motifs_in_group, nCpG = ont.nano.filt[[2]]$num_motifs_in_group, MReads = ont.nano.filt[[2]]$called_sites_methylated)

## Filter to sites within targetted regions


allGuides<-read.csv("../Resources/SmokingEWAS/GuideRNAsFINAL.csv")

targetRegions<-cbind(aggregate(allGuides$hg38, by = list(allGuides$Chr), min), aggregate(allGuides$hg38, by = list(allGuides$Chr), max)$x)
	
targetGRanges<-GRanges(seqnames = paste0("chr", targetRegions[,1]), strand = "*", ranges = IRanges(start = targetRegions[,2], end = targetRegions[,3]))



inTargets1<-subsetByOverlaps(sam1, targetGRanges)
inTargets2<-subsetByOverlaps(sam2, targetGRanges)

## can we detect differences between samples across smoking regions

## merge into granges object with shared sites
sharedSites<-findOverlaps(inTargets1, inTargets2)

combined<-inTargets1[queryHits(sharedSites),]
mcols(combined)<-data.frame("Non-Smoker" = combined$DNAm, "Smoker" = inTargets2$DNAm[subjectHits(sharedSites)], "Diff" =(inTargets2$DNAm[subjectHits(sharedSites)]-combined$DNAm), "Reads_NS" = combined$Reads, "Reads_Sm" = inTargets2$Reads[subjectHits(sharedSites)], "ncpg" = inTargets1$nCpG[queryHits(sharedSites)])

readCounts<-cbind(inTargets2$MReads[subjectHits(sharedSites)], 
		inTargets2$Reads[subjectHits(sharedSites)] - inTargets2$MReads[subjectHits(sharedSites)],
		inTargets1$MReads[queryHits(sharedSites)], 
		inTargets1$Reads[queryHits(sharedSites)] - inTargets1$MReads[queryHits(sharedSites)])
		
fishP<-rep(NA, nrow(readCounts))
for(i in 1:nrow(readCounts)){
	mat<-matrix(data = readCounts[i,], nrow = 2)
	fishP[i]<-fisher.test(mat)$p.value
	
}


fishP<-cbind(fishP, "Diff" = readCounts[,1]/(readCounts[,1]+readCounts[,2])-readCounts[,3]/(readCounts[,3]+readCounts[,4]), readCounts, mcols(combined), seqnames(combined), ranges(combined))

write.csv(fishP, "SmokingAnalysisNanopolishDNAm.csv")


pdf("Plots/FishersPAgainstReadDepth.pdf",height = 4, width = 8)
par(mfrow = c(1,2))
par(mar = c(4,4,1,1))
plot(fishP$Reads_NS+fishP$Reads_Sm, -log10(fishP$fishP), xlab = "Reads", ylab = "-log10P", pch = 16)
abline(h = -log10(0.05/514), lty = 2)
plot(fishP$Reads_NS+fishP$Reads_Sm, -log10(fishP$fishP), xlab = "Reads", ylab = "-log10P", pch = 16, xlim = c(0,100))
abline(h = -log10(0.05/514), lty = 2)
dev.off()

## cluster significant sites for plotting purposes
fishP.toplot<-fishP[which(fishP$fishP < 0.05/nrow(fishP)),]
toplot.granges<-GRanges(seqnames = fishP.toplot$"seqnames(combined)", ranges = fishP.toplot$"ranges(combined)")
toplot.granges<-reduce(toplot.granges, min.gapwidth=5000)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg38

axTrack <- GenomeAxisTrack()

## load EPIC EWAS smoking results
res.smoke<-read.csv("../Resources/SmokingEWAS/sst_bonf_all_dmps_USM1.csv", row.names = 1, stringsAsFactors =FALSE)
probeAnno<-read.table("/gpfs/mrc0/projects/Research_Project-MRC190311/references/EPICArray/EPIC.anno.GRCh38.tsv", header = TRUE, fill = TRUE)
probeAnno<-probeAnno[match(rownames(res.smoke), probeAnno$probeID),]
arrayData<-GRanges(seqnames = probeAnno$chrm, strand = "*", ranges = IRanges(start = probeAnno$start, end = probeAnno$end))
res.smoke<-cbind(res.smoke, probeAnno[match(rownames(res.smoke), probeAnno$probeID),c("chrm", "start")])
## load EPIC EWAS results
res.smoke<-res.smoke[which(res.smoke$chrm != "*"),]
res.smoke$chrm<-factor(as.character(res.smoke$chrm))
arrayData<-arrayData[which(seqnames(arrayData) != "*"),]

## load chromHMM annotations
stateAnno<-readPeakFile("/gpfs/mrc0/projects/Research_Project-MRC190311/references/ChromHMM/E062_15_coreMarks_hg38lift_mnemonics.bed.gz")
stateAnno<-subsetByOverlaps(stateAnno, targetGRanges)
chromHMM<-AnnotationTrack(stateAnno)
feature(chromHMM)<-unlist(lapply(strsplit(as.character(stateAnno$V4), "_"), tail, n = 1))

## plot full targetted regions

for(i in 1:length(targetGRanges)){

	chr<-as.character(seqnames(targetGRanges)[i])
	start<-start(targetGRanges)[i]
	stop<-end(targetGRanges)[i]

	windowSize<-stop-start
	if(windowSize > 1){
		start<-start-windowSize*0.05
		stop<-stop+windowSize*0.05

		idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)
		

		## create Grange sobject of EWAS results
		res.sub<-res.smoke[which(res.smoke$chrm == chr & res.smoke$start <= stop & res.smoke$start >= start),]

		## sensor pval that are 0
		res.sub$pval[which(res.sub$pval == 0)]<-min(res.sub$pval[which(res.sub$pval != 0)])/(10^10)

		res.sub<-res.sub[order(res.sub$start),]
		smokeP<-GRanges(seqnames = res.sub$chrm, strand = "*", ranges = IRanges(start = res.sub$start, end = res.sub$start), PackYearsP = -log10(res.sub$pval))
		smokeEffect<-GRanges(seqnames = res.sub$chrm, strand = "*", ranges = IRanges(start = res.sub$start, end = res.sub$start), Effect = res.sub$effect)

		refGenes <- UcscTrack(genome="hg38", chromosome=chr, table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=stop, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
								  symbol="name2", transcript="name", strand="strand", fill="#56B4E9",
								  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")					  
		cpgIslands <- UcscTrack(genome="hg38", chromosome=chr, track="cpgIslandExt", from=start, to=stop,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#009E73", name="CpG Islands")

		dnamDiff <- DataTrack(range = combined, genome = "hg38", type = "l", 
							 chromosome = chr, name = "Difference", data = fishP$Diff, baseline = 0, lwd = 2, col = "#D55E00", ylim = c(min(fishP$Diff),max(fishP$Diff)))
		dnamP <- DataTrack(range = combined, genome = "hg38", type = "p", 
							 chromosome = chr, name = "-log10P", data = -log10(fishP$fishP), col = "#D55E00", ylim = c(0, max(-log10(fishP$fishP))))				 

		if(nrow(res.sub) > 0){
		
			## create tracks of EPIC EWAS
			dTrack<-DataTrack(range = smokeP, genome = "hg38", type = "p", chromosome=chr, name = "-log10P", col = "#0072B2", ylim = c(0,max(-log10(res.sub$pval))))
			dTrack2<-DataTrack(range = smokeEffect, genome = "hg38", type = "l", chromosome=chr, name = "Difference", baseline = 0, col = "#0072B2", lwd = 2, ylim = c(min(c(0,res.sub$effect)), max(c(0,res.sub$effect))))
			
			pdf(paste0("Plots/SmokingAnalysisNanopolish", chr, "_FullRegion.pdf"), width = 10, height = 5)
			plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, chromHMM, dnamP,dnamDiff, dTrack, dTrack2), from=start, to=stop, showTitle=TRUE,  
			"TxWk" = rgb(0,0.4,0),     
			"Tx" = rgb(0,0.5,0),
			"EnhG" = rgb(0.76,0.88,0.5),
			"TssAFlnk" = rgb(1,0.69,0),
			"TssA" = rgb(1,0,0),
			"Quies" = rgb(1,1,1),
			"ReprPCWk" = rgb(0.75,0.75,0.75),
			"ReprPC" = rgb(0.5,0.5,0.5),
			"BivFlnk" = rgb(0.91,0.59,0.48),
			"EnhBiv"  = rgb(0.74,0.72,0.42),
			"TssBiv"  = rgb(0.8,0.36,0.36),
			"Enh" = rgb(1,1,0),
			"Het" = rgb(0.54,0.57,0.8), sizes = c(1,1.2,1.5,0.5,0.75,1.5,1.5,1.5,1.5))
			plot(0,1, type = "n", axes = FALSE, xlab = "", ylab = "")
			legend("center", box.col = "white", pch = 15, col = c(rgb(0,0.4,0),rgb(0,0.5,0),
			"EnhG" = rgb(0.76,0.88,0.5),rgb(1,0.69,0),
			"TssA" = rgb(1,0,0),rgb(1,1,1),rgb(0.75,0.75,0.75),rgb(0.5,0.5,0.5),rgb(0.91,0.59,0.48),rgb(0.74,0.72,0.42),rgb(0.8,0.36,0.36),rgb(1,1,0),rgb(0.54,0.57,0.8)), c("TxWk","Tx","EnhG","TssAFlnk","TssA","Quies","ReprPCWk","ReprPC","BivFlnk","EnhBiv","TssBiv","Enh","Het"),ncol = 5)
			dev.off()

	}

}

}

## plot zoomed in regions
for(i in 1:length(toplot.granges)){

	chr<-as.character(seqnames(toplot.granges)[i])
	start<-start(toplot.granges)[i]
	stop<-end(toplot.granges)[i]

	windowSize<-stop-start
	if(windowSize > 1){
		start<-start-windowSize*0.1
		stop<-stop+windowSize*0.1

		idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)
		

		## create Grange sobject of EWAS results
		res.sub<-res.smoke[which(res.smoke$chrm == chr & res.smoke$start <= stop & res.smoke$start >= start),]

		## sensor pval that are 0
		res.sub$pval[which(res.sub$pval == 0)]<-min(res.sub$pval[which(res.sub$pval != 0)])/(10^10)

		res.sub<-res.sub[order(res.sub$start),]
		smokeP<-GRanges(seqnames = res.sub$chrm, strand = "*", ranges = IRanges(start = res.sub$start, end = res.sub$start), PackYearsP = -log10(res.sub$pval))
		smokeEffect<-GRanges(seqnames = res.sub$chrm, strand = "*", ranges = IRanges(start = res.sub$start, end = res.sub$start), Effect = res.sub$effect)

		refGenes <- UcscTrack(genome="hg38", chromosome=chr, table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=stop, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
								  symbol="name2", transcript="name", strand="strand", fill="#56B4E9",
								  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")					  
		cpgIslands <- UcscTrack(genome="hg38", chromosome=chr, track="cpgIslandExt", from=start, to=stop,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#009E73", name="CpG Islands")

		dnamDiff <- DataTrack(range = combined, genome = "hg38", type = "l", 
							 chromosome = chr, name = "Difference", data = fishP$Diff, baseline = 0, lwd = 2, col = "#D55E00", ylim = c(min(fishP$Diff),max(fishP$Diff)))
		dnamP <- DataTrack(range = combined, genome = "hg38", type = "p", 
							 chromosome = chr, name = "-log10P", data = -log10(fishP$fishP), col = "#D55E00", ylim = c(0, max(-log10(fishP$fishP))))				 

		if(nrow(res.sub) > 0){
		
			## create tracks of EPIC EWAS
			dTrack<-DataTrack(range = smokeP, genome = "hg38", type = "p", chromosome=chr, name = "-log10P", col = "#0072B2", ylim = c(0,max(-log10(res.sub$pval))))
			dTrack2<-DataTrack(range = smokeEffect, genome = "hg38", type = "l", chromosome=chr, name = "Difference", baseline = 0, col = "#0072B2", lwd = 2, ylim = c(min(c(0,res.sub$effect)), max(c(0,res.sub$effect))))
			
			pdf(paste0("Plots/SmokingAnalysisNanopolish", chr, "_", floor(start), "_", floor(stop), ".pdf"), width = 10, height = 5)
			plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, chromHMM, dnamP,dnamDiff, dTrack, dTrack2), from=start, to=stop, showTitle=TRUE,  
			"TxWk" = rgb(0,0.4,0),     
			"Tx" = rgb(0,0.5,0),
			"EnhG" = rgb(0.76,0.88,0.5),
			"TssAFlnk" = rgb(1,0.69,0),
			"TssA" = rgb(1,0,0),
			"Quies" = rgb(1,1,1),
			"ReprPCWk" = rgb(0.75,0.75,0.75),
			"ReprPC" = rgb(0.5,0.5,0.5),
			"BivFlnk" = rgb(0.91,0.59,0.48),
			"EnhBiv"  = rgb(0.74,0.72,0.42),
			"TssBiv"  = rgb(0.8,0.36,0.36),
			"Enh" = rgb(1,1,0),
			"Het" = rgb(0.54,0.57,0.8), sizes = c(1,1.2,1.5,0.5,0.75,1.5,1.5,1.5,1.5))

			dev.off()
		} else {
			pdf(paste0("Plots/SmokingAnalysisNanopolish", chr, "_", floor(start), "_", floor(stop), ".pdf"), width = 10, height = 4)
			plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, chromHMM, dnamP,dnamDiff), from=start, to=stop, showTitle=TRUE,  
			"TxWk" = rgb(0,0.4,0),     
			"Tx" = rgb(0,0.5,0),
			"EnhG" = rgb(0.76,0.88,0.5),
			"TssAFlnk" = rgb(1,0.69,0),
			"TssA" = rgb(1,0,0),
			"Quies" = rgb(1,1,1),
			"ReprPCWk" = rgb(0.75,0.75,0.75),
			"ReprPC" = rgb(0.5,0.5,0.5),
			"BivFlnk" = rgb(0.91,0.59,0.48),
			"EnhBiv"  = rgb(0.74,0.72,0.42),
			"TssBiv"  = rgb(0.8,0.36,0.36),
			"Enh" = rgb(1,1,0),
			"Het" = rgb(0.54,0.57,0.8), sizes = c(1,1.2,1.5,0.5,0.75,1.5,1.5))
		}
	}

}


## match up EPIC EWAS and ONT analyses
ewasOverlaps<-findOverlaps(combined, arrayData)
matchedSites<-cbind(fishP[queryHits(ewasOverlaps),c("seqnames(combined)","ranges(combined)","Diff","fishP", "Reads_Sm", "Reads_NS")], rownames(res.smoke)[subjectHits(ewasOverlaps)], res.smoke[subjectHits(ewasOverlaps),c("chrm","start","effect","pval")])

write.csv(matchedSites, "OverlappingSitesONTEPICEWAS.csv")


## compare results between technologies
## colour points by which data significant in
colPoint<-rep("black", nrow(matchedSites))
colPoint[which(matchedSites$pval < 9e-8)]<-"#009E73"
colPoint[which(matchedSites$fishP < 0.05/nrow(fishP))]<-"#56B4E9"
colPoint[which(matchedSites$pval < 9e-8 & matchedSites$fishP < 0.05/nrow(fishP))]<-"#D55E00"

pdf("Plots/PvaluesONTAgainstEPIC.pdf",height = 4, width = 8)
par(mfrow = c(1,2))
plot(-log10(matchedSites$fishP), -log10(matchedSites$pval), pch = 16, xlab = "ONT", ylab = "EPIC array", col = colPoint)
abline(v = -log10(0.05/nrow(fishP)), lty = 2)
abline(h = -log10(9e-8), lty = 2)
plot(matchedSites$Diff, matchedSites$effect, pch = 16, xlab = "ONT", ylab = "EPIC array", col = colPoint)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1, lty = 2)
dev.off()

## are sites significant on EPIC array but not ONT have lower read depths?
pdf("Plots/ONTReadDepthAgainstEPICPvals.pdf",height = 4, width = 4)
plot(matchedSites$Reads_NS+matchedSites$Reads_Sm,-log10(matchedSites$pval),  pch = 16, xlab = "ONT Read depth", ylab = "EPIC array", col = colPoint)
abline(h = -log10(9e-8), lty = 2)
dev.off()

pdf("Plots/Figure4.pdf",height = 3, width = 9)
par(mfrow = c(1,3))
par(mar = c(4,4,1.5,1))
plot(-log10(matchedSites$fishP), -log10(matchedSites$pval), pch = 16, xlab = "ONT", ylab = "EPIC array", col = colPoint, cex.axis = 1.2, cex.lab = 1.2, cex = 1.2)
abline(v = -log10(0.05/nrow(fishP)), lty = 2)
abline(h = -log10(9e-8), lty = 2)
mtext("A", line = 0, adj = 0, side = 3, font = 2)
plot(matchedSites$Diff, matchedSites$effect, pch = 16, xlab = "ONT", ylab = "EPIC array", col = colPoint, cex.axis = 1.2, cex.lab = 1.2, cex = 1.2)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1, lty = 2)
mtext("B", line = 0, adj = 0, side = 3, font = 2)
plot(matchedSites$Reads_NS+matchedSites$Reads_Sm,-log10(matchedSites$pval),  pch = 16, xlab = "ONT Read depth", ylab = "EPIC array", col = colPoint, cex.axis = 1.2, cex.lab = 1.2, cex = 1.2)
abline(h = -log10(9e-8), lty = 2)
mtext("C", line = 0, adj = 0, side = 3, font = 2)
legend("topright", c("None", "EPIC", "ONT", "Both"), pch = 16, col = c("black", "#009E73", "#56B4E9", "#D55E00"), cex = 1.2)
dev.off()

