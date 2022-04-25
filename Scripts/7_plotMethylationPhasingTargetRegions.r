## calculate phasing statistics

calcLDPair<-function(site1,site2, minObs){
	
	## calculate number of shared reads
	n<-sum((!is.na(site1) & !is.na(site2)))
	
	if(n >= minObs){
		

		countSame<-sum(site1 == site2, na.rm = TRUE)
		countDiff<-sum(site1 != site2, na.rm = TRUE)

		probS1<-sum(site1, na.rm = TRUE)/sum(!is.na(site1))
		probS2<-sum(site2, na.rm = TRUE)/sum(!is.na(site2))

		probSame<-countSame/(countSame+countDiff)
		probDiff<-countSame/(countSame+countDiff)

		probMM<-probS1*probS2
		probUU<-(1-probS1)*(1-probS2)
		expectedSame<-probMM+probUU

		D<-probSame-expectedSame

		# standardize by maximum value
		if( D < 0){
			Dprime<-abs(D)/(probMM+probUU)
		} else {
			Dprime<-D/(1-(probMM+probUU))
		}

		r2<-cor(site1, site2, use="p")
		
		return(c(D,Dprime,r2, n))
	} else {
		return(c(NA,NA,NA, n))
	}
	
}

calcLDMat<-function(mat, minObs=10){
	nsites<-ncol(mat)
	print(paste(nsites, "sites to process"))
	out<-list("D" = matrix(data = NA, nrow = nsites, ncol = nsites),
	"Dprime" = matrix(data = NA, nrow = nsites, ncol = nsites),
	"r2" = matrix(data = NA, nrow = nsites, ncol = nsites),
	"N" = matrix(data = NA, nrow = nsites, ncol = nsites))
	
	for(i in 1:(nsites-1)){
		for(j in (i+1):nsites){
			tmp<-calcLDPair(mat[,i],mat[,j], minObs)
			out[["D"]][i,j]<-tmp[1]
			out[["Dprime"]][i,j]<-tmp[2]
			out[["r2"]][i,j]<-tmp[3]
			out[["N"]][i,j]<-tmp[4]
		}
	}
	return(out)
}

## create a matrix of reads against positions

## positive logLik indicates methylated
methClassifier<-function(value){
	if(length(value) == 1){
		if(value > 0){
			return(1)
		} else {
			return(0)
		}
	} else {
		return(NA)
	}
}

createMethMatrix<-function(dat, chr, start, stop){
	## filter dat to region requested
	dat<-dat[which(dat$chromosome == chr & dat$start <= stop & dat$start >= start),]
	mmat<-cast(dat, read_name ~ start, value="log_lik_ratio", methClassifier)
	rowNames<-mmat$read_name
	colNames<-colnames(mmat)[-1]
	mmat<-as.matrix(mmat[,-1])
	mode(mmat)<-"numeric"
	rownames(mmat)<-rowNames
	colnames(mmat)<-colNames
	return(mmat)
}

library(sqldf)
library(reshape)
library(GenomicRanges)
library(gplots)
library(LDheatmap)

setwd("DNAmethylation/")
allGuides<-read.csv("../Resources/SmokingEWAS/GuideRNAsFINAL.csv")

targetRegions<-cbind(aggregate(allGuides$hg38, by = list(allGuides$Chr), min), aggregate(allGuides$hg38, by = list(allGuides$Chr), max)$x)
	
targetGRanges<-GRanges(seqnames = paste0("chr", targetRegions[,1]), strand = "*", ranges = IRanges(start = targetRegions[,2], end = targetRegions[,3]))

ldStatsAll<-NULL
for(i in 1:length(targetGRanges)){

	chr<-as.character(seqnames(targetGRanges)[i])
	start<-start(targetGRanges)[i]
	stop<-end(targetGRanges)[i]

	windowSize<-stop-start

	sqlQuery<-paste0("select * from file where chromosome = '", chr,"' and start <= ", as.character(stop), " and start >= ", as.character(start))

	sam1<-read.csv.sql("Nanopolish/FAN62113/FAN62113_methylation_calls.tsv", sql = sqlQuery,  header = TRUE, sep = "\t")
	## filter calls with abs(loglik) < 2* numCpGs to be consistent with naonpolish
	sam1<-sam1[abs(sam1$log_lik_ratio)> 2*sam1$num_motifs,]

	sam2<-read.csv.sql("Nanopolish/FAO06764/FAO06764_methylation_calls.tsv", sql = sqlQuery,  header = TRUE, sep = "\t")
	sam2<-sam2[abs(sam2$log_lik_ratio)> 2*sam2$num_motifs,]
	
	sam<-rbind(sam1,sam2)
	
	mat<-createMethMatrix(sam, chr,start,stop)
	
	## work out read start and end CpGs
	readRange<-aggregate(start ~ read_name,sam,  range)

	## order
	readRange<-readRange[order(readRange$start[,1], readRange$start[,2]),]
	mat<-mat[readRange$read_name,]
		
	## filter out positions with < 10 reads
	mat<-mat[,which(colSums(!is.na(mat)) > 10)]
	## then filter reads left with no positions
	mat<-mat[which(rowSums(!is.na(mat)) > 0),]
	ldDat<-calcLDMat(mat, 10)
	colnames(ldDat$Dprime)<-colnames(mat)
	rownames(ldDat$Dprime)<-colnames(mat)
		
	write.csv(ldDat$Dprime, paste0("MethPhasingStatsTargettedRegions", chr, "_", floor(start), "_", floor(stop), ".csv"))
		
	ll<-LDheatmap(ldDat$Dprime, as.numeric(colnames(mat)), LDmeasure="D", flip = TRUE)

	pdf(paste0("Plots/LDHeatmapMethCallsTargettedRegions", chr, "_", floor(start), "_", floor(stop), ".pdf"), width = 10, height = 4)	
	LDheatmap.addGenes(ll, chr=chr, genome="hg38")
	dev.off()
	
	## split region into 15kb windows for plotting purposes
	if(windowSize > 15000){
		plotStart<-seq(min(as.numeric(colnames(ldDat$Dprime))), max(as.numeric(colnames(ldDat$Dprime))), 15000)
		for(i in 1:(length(plotStart)-1)){
			sitesIndex<-which(as.numeric(colnames(ldDat$Dprime)) <= plotStart[i]+15000 & as.numeric(colnames(ldDat$Dprime)) >= plotStart[i])
			ll<-LDheatmap(ldDat$Dprime[sitesIndex,sitesIndex], as.numeric(colnames(mat))[sitesIndex], LDmeasure="D", flip = TRUE)
			pdf(paste0("Plots/LDHeatmapMethCallsTargettedRegions", chr, "_", floor(start), "_", floor(stop), "Section", i, ".pdf"), width = 10, height = 4)
			LDheatmap.addGenes(ll, chr=chr, genome="hg38")
			dev.off()
		}
	}

	## calculate moving average
	distMat<-abs(outer(as.numeric(colnames(mat)), as.numeric(colnames(mat)), FUN="-"))

	ldStats<-unique(cbind("Dist" = c(distMat), "Dprime" = c(ldDat$Dprime), "r2" = c(ldDat$r2)))
	ldStats<-data.frame(ldStats[!is.na(ldStats[,2]),])

	## calculating rolling mean
	meanWindow<-100
	startPoints<-seq(0,max(ldStats$Dist), meanWindow/4)
	rollingMedian<-rep(NA, length(startPoints))
	for(i in 1:length(startPoints)){
		rollingMedian[i]<-median(ldStats$Dprime[which(ldStats$Dist <= startPoints[i]+meanWindow & ldStats$Dist >= startPoints[i])], na.rm = TRUE)
	}
	pdf(paste0("Plots/ScatterplotDprimeAgainstDistTargettedRegions", chr, "_", floor(start), "_", floor(stop), ".pdf"), width = 8, height = 6)	
	plot(ldStats$Dist, ldStats$Dprime, pch = 16, cex = 0.5, col = "grey", xlab = "Distance between sites (bp)", ylab =  "D'", xlim = c(0,10000))
	lines(startPoints+(meanWindow)/2, rollingMedian, lwd = 2)
	dev.off()
	}
	ldStatsAll<-rbind(ldStatsAll, ldStats)
	
}

## calculating rolling mean of D'
meanWindow<-50
startPoints<-seq(0,max(ldStatsAll$Dist), meanWindow/4)
rollingMedian<-rep(NA, length(startPoints))
for(i in 1:length(startPoints)){
	rollingMedian[i]<-median(ldStatsAll$Dprime[which(ldStatsAll$Dist <= startPoints[i]+meanWindow & ldStatsAll$Dist >= startPoints[i])], na.rm = TRUE)
}

pdf("Plots/ScatterplotDprimeAgainstDist.pdf", width = 8, height = 6)	
plot(ldStatsAll$Dist, ldStatsAll$Dprime, pch = 16, cex = 0.5, col = "grey", xlab = "Distance between sites (bp)", ylab =  "D'", xlim = c(0,10000))
lines(startPoints+(meanWindow)/2, rollingMedian, lwd = 2)
dev.off()

## calculating rolling mean of r2
startPoints<-seq(0,max(ldStatsAll$Dist), meanWindow/4)
rollingMedian<-rep(NA, length(startPoints))
for(i in 1:length(startPoints)){
	rollingMedian[i]<-median(abs(ldStatsAll$r2)[which(ldStatsAll$Dist <= startPoints[i]+meanWindow & ldStatsAll$Dist >= startPoints[i])], na.rm = TRUE)
}

pdf("Plots/ScatterplotR2AgainstDist.pdf", width = 8, height = 6)	
plot(ldStatsAll$Dist, abs(ldStatsAll$r2), pch = 16, cex = 0.5, col = "grey", xlab = "Distance between sites (bp)", ylab =  "R2", xlim = c(0,10000))
lines(startPoints+(meanWindow)/2, rollingMedian, lwd = 2)
dev.off()