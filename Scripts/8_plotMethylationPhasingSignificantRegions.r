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

## focus on regions detected with diffs
fishP<-read.csv("SmokingAnalysisNanopolishDNAm.csv", row.names = 1)

fishP.toplot<-fishP[which(fishP$fishP < 0.05/nrow(fishP)),]
fishP.granges<-GRanges(seqnames = fishP.toplot$seqnames.combined., ranges = IRanges(start = fishP.toplot$ranges.combined..start, end = fishP.toplot$ranges.combined..end))
toplot.granges<-reduce(fishP.granges, min.gapwidth=2000)


for(i in 1:length(toplot.granges)){

	chr<-as.character(seqnames(toplot.granges)[i])
	start<-start(toplot.granges)[i]
	stop<-end(toplot.granges)[i]

	windowSize<-stop-start
	if(windowSize > 100){
		start<-start-windowSize*0.1
		stop<-stop+windowSize*0.1
	} else {
		windowSize<-100
		start<-start-windowSize
		stop<-stop+windowSize
	}
	sqlQuery<-paste0("select * from file where chromosome = '", chr,"' and start <= ", as.character(stop), " and start >= ", as.character(start))

	sam1<-read.csv.sql("Nanopolish/FAN62113/FAN62113_methylation_calls.tsv", sql = sqlQuery,  header = TRUE, sep = "\t")
	## filter calls with abs(loglik) < 2* numCpGs to be consistent with naonpolish
	sam1<-sam1[abs(sam1$log_lik_ratio)> 2*sam1$num_motifs,]
	mat1<-createMethMatrix(sam1, chr,start,stop)

	sam2<-read.csv.sql("Nanopolish/FAO06764/FAO06764_methylation_calls.tsv", sql = sqlQuery,  header = TRUE, sep = "\t")
	sam2<-sam2[abs(sam2$log_lik_ratio)> 2*sam2$num_motifs,]
	mat2<-createMethMatrix(sam2, chr,start,stop)
	
	if(length(unique(sam1$start)) > 2 && length(unique(sam2$start)) > 2){

		## work out read start and end CpGs
		readRange1<-aggregate(start ~ read_name,sam1,  range)
		readRange2<-aggregate(start ~ read_name,sam2,  range)

		## order
		readRange1<-readRange1[order(readRange1$start[,1], readRange1$start[,2]),]
		readRange2<-readRange2[order(readRange2$start[,1], readRange2$start[,2]),]	

		mat1<-mat1[readRange1$read_name,]
		mat2<-mat2[readRange2$read_name,]
		

		## filter out positions with < 10 reads
		mat1<-mat1[,which(colSums(!is.na(mat1)) > 10)]
		mat2<-mat2[,which(colSums(!is.na(mat2)) > 10)]
		## then filter reads left with no positions
		mat1<-mat1[which(rowSums(!is.na(mat1)) > 0),]
		mat2<-mat2[which(rowSums(!is.na(mat2)) > 0),]
		
		## merge into single matrix for calculating phasing statistics
		commonCols<-intersect(colnames(mat1), colnames(mat2))
		mat1<-mat1[,commonCols]
		mat2<-mat2[,commonCols]
		
		mat<-rbind(mat1,mat2)
		ldDat<-calcLDMat(mat)
		colnames(ldDat$Dprime)<-colnames(mat)
		rownames(ldDat$Dprime)<-colnames(mat)
		
		write.csv(ldDat$Dprime, paste0("MethPhasingStatsSignificantRegions", chr, "_", floor(start), "_", floor(stop), ".csv"))
		fishPIndex<-subjectHits(findOverlaps(GRanges(chr, IRanges(colnames(mat))), GRanges(fishP$seqnames.combined., IRanges(fishP$ranges.combined..start, fishP$ranges.combined..end))))

		fishPIndex[fishP$fishP[fishPIndex]< 0.05/nrow(fishP)]
		ll<-LDheatmap(ldDat$Dprime, as.numeric(colnames(mat)), LDmeasure="D", flip = TRUE, SNP.name= fishP$ranges.combined..start[fishPIndex[fishP$fishP[fishPIndex]< 0.05/nrow(fishP)]])
		llGenes <- LDheatmap.addGenes(ll, chr=chr, genome="hg38")
		

		pdf(paste0("Plots/LDHeatmapMethCallsSignificantRegions", chr, "_", floor(start), "_", floor(stop), ".pdf"), width = 10, height = 4)	
		LDheatmap.addScatterplot(llGenes, P=-log10(fishP$fishP[fishPIndex]))
		dev.off()
	

	}
}