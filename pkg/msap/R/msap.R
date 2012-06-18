# msap - Statistical analysis for Methilattion-Sensitive Amplification Polimorphism data
# version: 0.1-2
# Author: Andrés Pérez-Figueroa (anpefi@uvigo.es)




msap <- function(datafile, name=datafile, uninformative=TRUE, nDec=4, meth=TRUE, rm.redundant=TRUE, rm.monomorphic=TRUE, do.pcoa=TRUE, do.shannon=TRUE, do.AMOVA=TRUE, do.pairwisePhiST=TRUE){
	
	cat("\nmsap - Statistical analysis for Methilation-Sensitive Amplification Polimorphism data\n")

	#Read datafile
	cat("\nReading ", datafile,"\n")
	data <- read.csv(datafile, header=TRUE)
	
	#get the loci names
	locus <- rownames(data)[4:length(data[1,])]
	cat("Number of loci: ",length(locus),"\n")
	#sorting
	data <- data[with(data, order(data[,1],data[,2],data[,3])),]
	
	
	if(meth){ #meth=TRUE -> MSAP
		groups <- data[data[,3]=="HPA",1]
		ntt <- length(levels(groups))
		cat("Number of samples/individuals: ",length(data[,1])/2,"\n")
		cat("Number of groups/populations: ",ntt,"\n")
		dataMSP <- data[data[,3]=="MSP",4:length(data[1,])]
		dataHPA <- data[data[,3]=="HPA",4:length(data[1,])]
		dataMIX <- dataHPA*10+dataMSP
		Met <- apply(dataMIX, 2, methStatusEval, type=uninformative)
		MSL <- which(Met)
		NML <- which(!Met)
		MSL.nloci <- length(Met[MSL])
		NML.nloci <- length(Met[NML])
		cat("Number of Methylation-Susceptible Loci (MSL): ",MSL.nloci,"\n")
		cat("Number of No Methylated Loci (NML): ",NML.nloci,"\n\n")
		
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
			matM <- as.matrix(dataMIXb[,MSL])
			matN <- as.matrix(dataMIXb[,NML])
		}
		
		if(rm.monomorphic){
			Poly <- apply(matM, 2, polymorphic)
			matM <- matM[,Poly]
			Poly <- apply(matN, 2, polymorphic)
			matN <- matN[,Poly]
			MSL.ploci <- length(matM[1,])
			NML.ploci <- length(matN[1,])
			cat("Number of polymorphic MSL: ",MSL.ploci," (",format(MSL.ploci/MSL.nloci*100,digits=1),"% of total MSL)\n")
			cat("Number of polymorphic NML: ",NML.ploci," (",format(NML.ploci/NML.nloci*100,digits=1),"% of total NML)\n\n")
		}
		
		MSL.I <- apply(matM, 2, shannon)
		NML.I <- apply(matN, 2, shannon)
		cat("\nShannon's Diversity Index \n")
		cat("MSL: I = ", mean(MSL.I, na.rm=T),"  (SD: ", sd(MSL.I, na.rm=T),")\n")
		cat("NML: I = ", mean(NML.I, na.rm=T),"  (SD: ", sd(NML.I, na.rm=T),")\n")
		wt<-wilcox.test(MSL.I,NML.I)
		pval<-ifelse(wt$p.value<0.0001, "P < 0.0001", paste("P = ",wt$p.value))
		cat(wt$method,": ", names(wt$statistic)," = ",wt$statistic[[1]]," (",pval,")\n")
		png(filename=paste(name,"Shannon.png", sep='-'), width=680, height=680)
		boxplot(MSL.I, NML.I, names=c("MSL","NSL"), ylab="Shannon's I")
	
		dev.off()
	
	
		### MSL 
		cat("\n\n*****************************\nAnalysis of MSL\n")
		repMet(dataMIX[,MSL], groups, nDec)
		
		#shannon for every group
	
		DM<-lingoes(as.dist(smc(matM, dist=TRUE))) #Simple Matching Coefficient, from scrime
		#DM <- as.dist(DM)
		#DM <- lingoes(DM)

		#PCoA
		pcoa(DM, groups, name, "MSL")

		#AMOVA
		cat("\nPerforming AMOVA\n")
		diffAmova(DM, groups, nDec)
	
	} #end if meth
	else{  #meth=FALSE -> AFLP
		groups <- data[,1] #Here, all rows are indepedent samples
		ntt <- length(levels(groups))
		matN <- as.matrix(data[,4:length(data[1,])])
		matN[matN==1]<-2
		matN[matN==0]<-1 
		cat("Number of samples/individuals: ",length(data[,1]),"\n")
		cat("Number of groups/populations: ",ntt,"\n")
	}
	#### NML
		
	if(meth) cat("\n\n*****************************\nAnalysis of NML\n")
	DM<-smc(matN, dist=TRUE) #Simple Matching Coefficient, from scrime
	DM <- as.dist(DM)
	DM <- lingoes(DM)

	#PCoA
	pcoa(DM, groups, name, "NML")

	#AMOVA
	cat("\nPerforming AMOVA\n")
	diffAmova(DM, groups, nDec)
	
	
	
	
	 
}







