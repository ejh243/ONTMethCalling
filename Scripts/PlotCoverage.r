
## nb guides were designed in hg19 
## ONT aligned to hg38

library(Gviz)
library(vioplot)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# need a biostrings object for reference
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

sam1<-read.table("SumStats/FAN62113.bedg")
sam1<-GRanges(seqnames = sam1$V1, strand = "*", ranges = IRanges(start = sam1$V2, end = sam1$V3), Depth = sam1$V4)

sam2<-read.table("SumStats/FAO06764.bedg")
sam2<-GRanges(seqnames = sam2$V1, strand = "*", ranges = IRanges(start = sam2$V2, end = sam2$V3), Depth = sam2$V4)

## plot coverage within target regions
allGuides<-read.csv("../../Resources/SmokingEWAS/GuideRNAsFINAL.csv")

axTrack <- GenomeAxisTrack()

## these coords have been lifted over
chr<-5
start<-304177-10000
stop<-438290+10000

idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)

refGenes <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=stop, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
						  symbol="name2", transcript="name", strand="strand", fill="#8282d2",
						  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")					  
cpgIslands <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), track="cpgIslandExt", from=start, to=stop,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#006400", name="CpG Islands")

dTrack1 <- DataTrack(range = sam1, genome = "hg38", type = "l", 
                     chromosome = chr, name = "Sample 1")
dTrack2 <- DataTrack(range = sam2, genome = "hg38", type = "l", 
                     chromosome = chr, name = "Sample 2")
					 
ahrr<-allGuides[which(allGuides$Chr == 5),]
arrowSize<-10000
ahrr$Start<-ahrr$hg38
ahrr$End<-ahrr$hg38
ahrr$Start[ahrr$Strand == "-1"]<-ahrr$Start[ahrr$Strand == "-1"]-arrowSize
ahrr$End[ahrr$Strand == "1"]<-ahrr$Start[ahrr$Strand == "1"]+arrowSize

ahrrGuides<-GRanges(seqnames = chr, strand = ahrr$Strand, ranges = IRanges(start = ahrr$Start, end = ahrr$End))

aTrack<-AnnotationTrack(ahrrGuides, genome = "hg38", name = "Crispr guides")

pdf("Plots/CoverageHistogram_AHRR.pdf", width = 10, height = 4)
plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, dTrack1, dTrack2, aTrack), from=start, to=stop, showTitle=TRUE)
dev.off()



chr<-2
start<-232418619-10000
stop<-232420224+10000


idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)

refGenes <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=stop, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
						  symbol="name2", transcript="name", strand="strand", fill="#8282d2",
						  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")					  
cpgIslands <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), track="cpgIslandExt", from=start, to=stop,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#006400", name="CpG Islands")

dTrack1 <- DataTrack(range = sam1, genome = "hg38", type = "l", 
                     chromosome = chr, name = "Sample 1")
dTrack2 <- DataTrack(range = sam2, genome = "hg38", type = "l", 
                     chromosome = chr, name = "Sample 2")
					 
inter<-allGuides[which(allGuides$Chr == 2),]

## in order to draw arrows need to create an arbitraryily sized interval. Change this parameter to influence size of arrow that is drawn

arrowSize<-1000
inter$Start<-inter$hg38
inter$End<-inter$hg38
inter$Start[inter$Strand == "-1"]<-inter$Start[inter$Strand == "-1"]-arrowSize
inter$End[inter$Strand == "1"]<-inter$Start[inter$Strand == "1"]+arrowSize

interGuides<-GRanges(seqnames = chr, strand = inter$Strand, ranges = IRanges(start = inter$Start, end = inter$End))

aTrack<-AnnotationTrack(interGuides, genome = "hg38", name = "Crispr guides")

pdf("Plots/CoverageHistogram_Chr2.pdf", width = 10, height = 4)
plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, dTrack1, dTrack2, aTrack), from=start, to=stop, showTitle=TRUE)
dev.off()

chr<-1
start<-92474761-10000
stop<-92486876+10000


idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)

refGenes <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=stop, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
						  symbol="name2", transcript="name", strand="strand", fill="#8282d2",
						  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")					  
cpgIslands <- UcscTrack(genome="hg38", chromosome=paste("chr", chr, sep = ""), track="cpgIslandExt", from=start, to=stop,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#006400", name="CpG Islands")

dTrack1 <- DataTrack(range = sam1, genome = "hg38", type = "l", 
                     chromosome = chr, name = "Sample 1")
dTrack2 <- DataTrack(range = sam2, genome = "hg38", type = "l", 
                     chromosome = chr, name = "Sample 2")
					 
gfi1<-allGuides[which(allGuides$Chr == 1),]
arrowSize<-1000
gfi1$Start<-gfi1$hg38
gfi1$End<-gfi1$hg38
gfi1$Start[gfi1$Strand == "-1"]<-gfi1$Start[gfi1$Strand == "-1"]-arrowSize
gfi1$End[gfi1$Strand == "1"]<-gfi1$Start[gfi1$Strand == "1"]+arrowSize

gfi1Guides<-GRanges(seqnames = chr, strand = gfi1$Strand, ranges = IRanges(start = gfi1$Start, end = gfi1$End))

aTrack<-AnnotationTrack(gfi1Guides, genome = "hg38", name = "Crispr guides")

pdf("Plots/CoverageHistogram_GFI1.pdf", width = 10, height = 4)
plotTracks(list(idxTrack, axTrack, refGenes,cpgIslands, dTrack1, dTrack2, aTrack), from=start, to=stop, showTitle=TRUE)
dev.off()


## calc coverage by 
targetRegions<-cbind(aggregate(allGuides$hg38, by = list(allGuides$Chr), min), aggregate(allGuides$hg38, by = list(allGuides$Chr), max)$x)
	
targetGRanges<-GRanges(seqnames = paste0("chr", targetRegions[,1]), strand = "*", ranges = IRanges(start = targetRegions[,2], end = targetRegions[,3]))

inTargets1<-subsetByOverlaps(sam1, targetGRanges)
inTargets2<-subsetByOverlaps(sam2, targetGRanges)


pdf("Plots/BoxplotDepthGenome-wide.pdf")

par(mar = c(4, 4, 4, 0.5))
boxplot(sam1$Depth, sam2$Depth, ylab = "Depth", names = c("Sample 1", "Sample 2"), xlab = "", main = "Genome-wide")


dev.off()

pdf("Plots/BoxplotDepthWithTargetedRegions.pdf")
par(mfrow = c(2,1))
par(mar = c(0.5, 4, 4, 0.5))
boxplot(as.numeric(inTargets1$Depth) ~ as.character(seqnames(inTargets1)), ylab = "Depth", names = c("", "", ""), xlab = "", main = "Sample 1" , at = 2:4, xlim = c(0.5, 4.5))
boxplot(sam1$Depth, add = TRUE)
par(mar = c(4, 4, 4, 0.5))
boxplot(as.numeric(inTargets2$Depth) ~ as.character(seqnames(inTargets2)), ylab = "Depth", names = c("Genome-wide", "GFI1", "Chr2_I", "AHRR"), xlab = "", main = "Sample 2", at = 2:4, xlim = c(0.5, 4.5))
boxplot(sam2$Depth, add = TRUE)
dev.off()

pdf("Plots/ViolinplotDepthWithTargetedRegions.pdf")
par(mfrow = c(2,1))
par(mar = c(0.5, 4, 4, 0.5))
vioplot(as.numeric(inTargets1$Depth) ~ as.character(seqnames(inTargets1)), ylab = "Depth", names = c("", "", ""), xlab = "", main = "Sample 1" , at = 2:4, xlim = c(0.5, 4.5))
vioplot(sam1$Depth, add = TRUE)
par(mar = c(4, 4, 4, 0.5))
vioplot(as.numeric(inTargets2$Depth) ~ as.character(seqnames(inTargets2)), ylab = "Depth", names = c("GFI1", "Chr2_I", "AHRR"), xlab = "", main = "Sample 2", at = 2:4, xlim = c(0.5, 4.5))
vioplot(sam2$Depth, add = TRUE)
dev.off()



