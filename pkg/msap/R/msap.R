# msap - Statistical analysis for Methilattion-Sensitive Amplification Polimorphism data
# version: 0.1-2
# Author: Andrés Pérez-Figueroa (anpefi@uvigo.es)




msap <- function(datafile, name=datafile, uninformative=TRUE, nDec=4){
	#loading required packages
	#suppressPackageStartupMessages(require(ade4, warn.conflicts=FALSE))
	#suppressPackageStartupMessages(require(scrime, warn.conflicts=FALSE))
	#suppressPackageStartupMessages(require(pegas, warn.conflicts=FALSE))
	#suppressPackageStartupMessages(require(cba, warn.conflicts=FALSE))
	
	cat("\nmsap - Statistical analysis for Methilation-Sensitive Amplification Polimorphism data\n")

	#Read datafile
	cat("\nReading ", datafile,"\n")
	data <- read.csv(datafile, header=TRUE)
	#data file should have these columns:
	# 1. Label for factor levels (one-factor only)
	# 2. Arbitrary label
	# 3. Enzyme labels: HPA or MSP
	# 4+. band presence/absence: 1/0
	
	#get the loci names
	locus <- rownames(data)[4:length(data[1,])]
	#Set the first columns
	#sorting
	data <- data[with(data, order(data[,1],data[,2],data[,3])),]
	
	

	groups <- data[data[,3]=="HPA",1]

	ntt <- length(levels(groups))
	dataMSP <- data[data[,3]=="MSP",4:length(data[1,])]
	dataHPA <- data[data[,3]=="HPA",4:length(data[1,])]

	dataMIX <- dataHPA*10+dataMSP
	Met <- apply(dataMIX, 2, methStatusEval, type=uninformative)
	MSL <- which(Met)
	NML <- which(!Met)
	cat("Number of Methylation-Susceptible Loci (MSL): ",length(Met[MSL]),"\n")
	cat("Number of No Methylated Loci (NML): ",length(Met[NML]),"\n\n")
	
	
	
	if (uninformative){
	#This transformation assumes that HPA-/MSP- pattern is uninformative as it could represent full methylation of cytosines in the target or that target is missing by mutation. 
		#Herrera & Bazaga 2010
		
		
		dataMSL <- dataMIX[,MSL]
		dataNML <- dataMIX[,NML]
		dataMSL[dataMSL==0] <- NA
		dataMSL[dataMSL==1] <- 2
		dataMSL[dataMSL==11] <- 1
		dataMSL[dataMSL==10] <- 2
		dataNML <- ifelse(dataNML==0, 1, 2)
		
		matM <- as.matrix(dataMSL)
		matN <- as.matrix(dataNML)
	}
	else
	{
		#This transformation assumes that HPA-/MSP- pattern represents full methylation of cytosines in the target. Thus, all bands are scored as methylation.
		#Lu et al. 2008; Gupta et al. 2012

		dataMIXb <- ifelse(dataMIX==11, 2, 1)
		matM <- as.matrix(dataMIXb)
	}
	MSL.I <- apply(matM, 2, shannon)
	NML.I <- apply(matN, 2, shannon)
	cat("\nShannon's Diversity Index\n")
	cat("MSL: I = ", mean(MSL.I, na.rm=T),"  (SD: ", sd(MSL.I, na.rm=T),")\n")
	cat("NML: I = ", mean(NML.I, na.rm=T),"  (SD: ", sd(NML.I, na.rm=T),")\n")
	wt<-wilcox.test(MSL.I,NML.I)
	pval<-ifelse(wt$p.value<0.0001, "P < 0.0001", paste("P = ",wt$p.value))
	cat(wt$method,": ", names(wt$statistic)," = ",wt$statistic[[1]]," (",pval,")\n")
	png(filename=paste(name,"Shannon.png", sep='-'), width=680, height=680)
	boxplot(MSL.I, NML.I, names=c("MSL","NSL"), ylab="Shannon's I")
	
	dev.off()
	
	
	### MSL 
	cat("\nAnalysis of MSL\n")
	repMet(dataMIX[,MSL], groups, nDec)
	
	
	DM<-smc(matM, dist=TRUE) #Simple Matching Coefficient
	DM <- as.dist(DM)
	DM <- lingoes(DM)

	#PCoA
	pcol <- dudi.pco(DM, scannf = FALSE, nf = 2)
	2
	var <- pcol$eig / sum(pcol$eig) *100
	var
	var1 <- round(var[1], digits=1)
	var2 <- round(var[2], digits=1)
	spcoo<-split(pcol$li, groups)
	maxX <- max(pcol$li$A1)
	minX <- min(pcol$li$A1)
	maxY <- max(pcol$li$A2)
	minY <- min(pcol$li$A2)
	png(filename=paste(name,"MSL.png",sep='-'), width=680, height=680)
	par(bty = 'n')
	plot(0,0, main=paste(name," (MSL)"), type = "n",xlab=paste("C1 (",var1,"%)"),ylab=paste("C2 (",var2,"%)"), xlim=c(minX,maxX+0.1), ylim=c(minY,maxY+0.1), frame=TRUE, cex=1.5)
	bgcolors<-c("black", "green", "yellow", "red", "blue", "orange", "pink")
	symbs <- c(21,22,23,24,25,26,27)
	for(i in 1:ntt){
		points(spcoo[[i]], pch=symbs[i], col="black", bg=bgcolors[i])
	}	
	s.class(pcol$li, groups, cpoint=0, add.plot=TRUE)
	dev.off()

	#AMOVA
	cat("\nPerforming AMOVA\n")
	assign("DM", DM, envir=globalenv())
	assign("groups", groups, envir=globalenv())
	amv<-amova(DM ~ groups, nperm=10000)
	cat("AMOVA TABLE \td.f. \tSSD \t\tMSD \t\tVariance\n")
	phiST <- as.numeric(amv$varcomp[1,1]/(amv$varcomp[1,1]+amv$varcomp[2,1]))
	cat("among groups\t",amv$tab[1,3],"\t",format(amv$tab[1,1], digits=nDec),"\t",format(amv$tab[1,2], digits=nDec),"\t",format(amv$varcomp[1,1], digits=nDec),"\n")
	cat("within groups\t",amv$tab[2,3],"\t",format(amv$tab[2,1], digits=nDec),"\t",format(amv$tab[2,2], digits=nDec),"\t",format(amv$varcomp[2,1], digits=nDec),"\n")
	cat("Total        \t",amv$tab[3,3],"\t",format(amv$tab[3,1], digits=nDec),"\t",format(amv$tab[3,2], digits=nDec), " \n")
	cat("\n")
	pval<- amv$varcomp[1,2]
	pval<- ifelse(pval<0.0001, "(P<0.0001)", paste("(P=",format(pval,digits=nDec),")"))
	cat("Phi_ST = ", format(phiST, digits=nDec), " ", pval,"\n")

	cat("\nPairwise Phi_ST\n------------------------------------\n")

	for (i in 1:(ntt-1) ){
		for(j in seq(i+1,ntt ) ) {
			p1 <- which(groups==levels(groups)[i])
			p2 <- which(groups==levels(groups)[j])
			assign("ttos", groups[c(p1,p2)], envir=globalenv())
			assign("M", subset(DM,c(p1,p2)), envir=globalenv())
			a <- amova(M ~ ttos, nperm=10000)
			pval<-as.numeric(a$varcomp$P.value[1])

			pval<-ifelse(pval<0.0001, "(P<0.0001)", paste("(P=",format(pval, digits=nDec),")"))
			phiST <- as.numeric(a$varcomp[1,1]/(a$varcomp[1,1]+a$varcomp[2,1]))
			
			cat(levels(groups)[i], " - ",levels(groups)[j],": ",format(phiST,trim=T,digits=nDec),"\t ",pval," \n")
		}
	}
	
	
	
	#### NML
	
	cat("\nAnalysis of NML\n")
	DM<-smc(matN, dist=TRUE) #Simple Matching Coefficient
	DM <- as.dist(DM)
	DM <- lingoes(DM)

	#PCoA
	pcol <- dudi.pco(DM, scannf = FALSE, nf = 2)
	2
	var <- pcol$eig / sum(pcol$eig) *100
	var
	var1 <- round(var[1], digits=1)
	var2 <- round(var[2], digits=1)
	spcoo<-split(pcol$li, groups)
	maxX <- max(pcol$li$A1)
	minX <- min(pcol$li$A1)
	maxY <- max(pcol$li$A2)
	minY <- min(pcol$li$A2)
	png(filename=paste(name,"nML.png",sep='-'), width=680, height=680)
	par(bty = 'n')
	plot(0,0, main=paste(name," (NML)"), type = "n",xlab=paste("C1 (",var1,"%)"),ylab=paste("C2 (",var2,"%)"), xlim=c(minX,maxX+0.1), ylim=c(minY,maxY+0.1), frame=TRUE, cex=1.5)
	bgcolors<-c("black", "green", "yellow", "red", "blue", "orange", "pink")
	symbs <- c(21,22,23,24,25,26,27)
	for(i in 1:ntt){
		points(spcoo[[i]], pch=symbs[i], col="black", bg=bgcolors[i])
	}	
	s.class(pcol$li, groups, cpoint=0, add.plot=TRUE)
	dev.off()

	#AMOVA
	cat("\nPerforming AMOVA\n")
	assign("DM", DM, envir=globalenv())
	assign("groups", groups, envir=globalenv())
	amv<-amova(DM ~ groups, nperm=10000)
	cat("AMOVA TABLE \td.f. \tSSD \t\tMSD \t\tVariance\n")
	phiST <- as.numeric(amv$varcomp[1,1]/(amv$varcomp[1,1]+amv$varcomp[2,1]))
	cat("among groups\t",amv$tab[1,3],"\t",format(amv$tab[1,1], digits=nDec),"\t",format(amv$tab[1,2], digits=nDec),"\t",format(amv$varcomp[1,1], digits=nDec),"\n")
	cat("within groups\t",amv$tab[2,3],"\t",format(amv$tab[2,1], digits=nDec),"\t",format(amv$tab[2,2], digits=nDec),"\t",format(amv$varcomp[2,1], digits=nDec),"\n")
	cat("Total        \t",amv$tab[3,3],"\t",format(amv$tab[3,1], digits=nDec),"\t",format(amv$tab[3,2], digits=nDec), " \n")
	cat("\n")
	pval<- amv$varcomp[1,2]
	pval<- ifelse(pval<0.0001, "(P<0.0001)", paste("(P=",format(pval,digits=nDec),")"))
	cat("Phi_ST = ", format(phiST, digits=nDec), " ", pval,"\n")

	cat("\nPairwise Phi_ST\n------------------------------------\n")

	for (i in 1:(ntt-1) ){
		for(j in seq(i+1,ntt ) ) {
			p1 <- which(groups==levels(groups)[i])
			p2 <- which(groups==levels(groups)[j])
			assign("ttos", groups[c(p1,p2)], envir=globalenv())
			assign("M", subset(DM,c(p1,p2)), envir=globalenv())
			a <- amova(M ~ ttos, nperm=10000)
			pval<-as.numeric(a$varcomp$P.value[1])

			pval<-ifelse(pval<0.0001, "(P<0.0001)", paste("(P=",format(pval, digits=nDec),")"))
			phiST <- as.numeric(a$varcomp[1,1]/(a$varcomp[1,1]+a$varcomp[2,1]))
			
			cat(levels(groups)[i], " - ",levels(groups)[j],": ",format(phiST,trim=T,digits=nDec),"\t ",pval," \n")
		}
	}
	
	
}


repMet <-function(dataM,groups,nDec){
	cat("Report of methylation levels \n")
	dataM[dataM==11]<-"u"
	dataM[dataM==10]<-"h"
	dataM[dataM==1]<-"i"
	dataM[dataM==0]<-"f"
	dataM <- split(dataM,groups)
	
	ns<-names(dataM)
	res<-matrix(ncol=length(ns), nrow=4)
	colnames(res)<-ns
	rownames(res)<-c("HPA+/MSP+ (Unmethylated)", "HPA+/MSP- (Hemimethylated)","HPA-/MSP+ (Internal cytosine methylation)","HPA-/MSP- (Full methylation or absence of target)")
	for(x in seq_along(dataM)){
	res[1,x]<- as.numeric(table(as.matrix(dataM[[x]]))["u"])
	res[2,x]<- as.numeric(table(as.matrix(dataM[[x]]))["h"])
	res[3,x]<- as.numeric(table(as.matrix(dataM[[x]]))["i"])
	res[4,x]<- as.numeric(table(as.matrix(dataM[[x]]))["f"])
	suma<- sum(res[,x])
	res[1,x] <- res[1,x]/suma
	res[2,x] <- res[2,x]/suma
	res[3,x] <- res[3,x]/suma
	res[4,x] <- res[4,x]/suma
	}
	print(res, digits=nDec)
}

shannon <-function(x){
	#Calculate Shannon's Diversity Index given a list of binomial alleles
	t <- table(x)
	if (length(t)<2) return(NaN)
	P=t[1]/sum(t)
	Q=t[2]/sum(t)
	S=-(P*log(P)+Q*log(Q))
	return(S)
}

methStatusEval <- function(x, error=0.05, type=TRUE){
	#Evaluate for classification as 'methylation-susceptible locus' or 'non methylated"
	#Uses mixed data: 11, 10, 1, 0 are the four possible states for both enzymes
	threshold= ((error + error) - (2*error*error))
	#A locus is classified as "methylation-susceptible" if the frecuency of methylation evidence (10, 1 or 0) is above threshold
	if(type) fMet <- 1-(length(x[which(x==11|x==0)])/length(x))
	else fMet <- 1-(length(x[which(x==11)])/length(x))
	return(fMet>threshold)
}

