# msap - Statistical analysis for Methilattion-Sensitive Amplification Polimorphism data
# version: 0.1-2
# Author: Andrés Pérez-Figueroa (anpefi@uvigo.es)




msap <- function(datafile, name=datafile, uninformative=TRUE, nDec=4, meth=TRUE, rm.redundant=TRUE, rm.monomorphic=TRUE, do.pcoa=TRUE, do.shannon=TRUE, do.amova=TRUE, do.pairwisePhiST=TRUE, do.cluster=TRUE, use.groups=NULL){
	
	cat("\nmsap - Statistical analysis for Methilation-Sensitive Amplification Polimorphism data\n")

	#Read datafile
	cat("\nReading ", datafile,"\n")
	data <- read.csv(datafile, header=TRUE)
	
	#get the loci names
	locus <- colnames(data)[4:length(data[1,])]
	cat("Number of loci: ",length(locus),"\n")
	#sorting
	data <- data[with(data, order(data[,1],data[,2],data[,3])),]
	
	
	#If user just want to use some of the groups in the datafile
	#Should be provided with the same name that in file
	if(is.vector(use.groups)){
		if(is.character(use.groups)){
			data <- data[which(is.element(data[,1],use.groups)),]
			cat("\nUsing a subset of the data:\n")
			print(levels(factor(data[,1]))) 
		}
		else cat("\nWARNING: use.groups argument wrong. Using all groups.\n\n")
	}
	
	
	if(meth){ #meth=TRUE -> MSAP
		inds <- as.character(data[data[,3]=="HPA",2])
		groups <- factor(data[data[,3]=="HPA",1])  #Needed to refactor here 
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
		
		PolyM <- apply(matM, 2, polymorphic)
		PolyN <- apply(matN, 2, polymorphic)
		MSL.ploci <- length(which(PolyM))
		NML.ploci <- length(which(PolyN))
		cat("Number of polymorphic MSL: ",MSL.ploci," (",format(MSL.ploci/MSL.nloci*100,digits=1),"% of total MSL)\n")
		cat("Number of polymorphic NML: ",NML.ploci," (",format(NML.ploci/NML.nloci*100,digits=1),"% of total NML)\n\n")
		NML.nloci <- NML.ploci
		MSL.nloci <- MSL.ploci
		if(rm.monomorphic){
			matM <- matM[,PolyM]
			if(NML.nloci>0) matN <- matN[,PolyN]
			#Note: 
		}
		
		MSL.I <- apply(matM, 2, shannon)
		if(NML.nloci>0) NML.I <- apply(matN, 2, shannon)
		cat("\nShannon's Diversity Index \n")
		cat("MSL: I = ", mean(MSL.I, na.rm=T),"  (SD: ", sd(MSL.I, na.rm=T),")\n")
		if(NML.nloci>0){
		cat("NML: I = ", mean(NML.I, na.rm=T),"  (SD: ", sd(NML.I, na.rm=T),")\n")
		wt<-wilcox.test(MSL.I,NML.I)
		pval<-ifelse(wt$p.value<0.0001, "P < 0.0001", paste("P = ",wt$p.value))
		cat(wt$method,": ", names(wt$statistic)," = ",wt$statistic[[1]]," (",pval,")\n")
		#png(filename=paste(name,"Shannon.png", sep='-'), width=680, height=680)
		#boxplot(MSL.I, NML.I, names=c("MSL","NSL"), ylab="Shannon's I")
		#dev.off()
		}
	
		### MSL 
		cat("\n\n*****************************\nAnalysis of MSL\n")
		repMet(dataMIX[,MSL], groups, nDec)
		
		#shannon for every group
	
		DM<-lingoes(as.dist(smc(matM, dist=TRUE))) #Simple Matching Coefficient, from scrime
		
		if(do.cluster){
			# cluster
			DM_copy <- DM
			attr(DM_copy, "Labels") <- inds
			np <- table(groups)[]
			MSL.cluster <-nj(DM_copy) 
			darksch <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#A6761D","#666666","#E64AB02")
			tipCol<-rep(darksch[1:length(np)],np)
			ecol<-unlist(lapply(MSL.cluster$edge[,2], edgeCol, length(MSL.cluster$tip.label), tipCol))
			png(filename=paste(name,"MSL-NJ.png", sep='-'))
			plot.phylo(MSL.cluster, tip.color=rep(darksch[1:length(np)],np), use.edge.length=T, edge.color=ecol, edge.width=3, show.tip.label=T, main="MSL")
			legend("bottomright", as.character(levels(groups)), col=darksch, lwd=3)
			dev.off()
		}
		
		#export 
		#assign("MSLmat", DM, envir=globalenv())
		#assign("pops", groups, envir=globalenv())
		#PCoA
		if(do.pcoa)pcoa(DM, groups, name, "MSL")
		assign("DM.MSL", DM, envir=globalenv())
		assign("pops", groups, envir=globalenv())
		if(do.amova){
			#AMOVA
			cat("\nPerforming AMOVA\n")
			diffAmova(DM, groups, nDec, do.pairwisePhiST)
		}
	
	} #end if meth
	else{  #meth=FALSE -> AFLP
		groups <- factor(data[,1]) #Here, all rows are indepedent samples
		ntt <- length(levels(groups))
		matN <- as.matrix(data[,4:length(data[1,])])
		matN[matN==1]<-2
		matN[matN==0]<-1 
		cat("Number of samples/individuals: ",length(data[,1]),"\n")
		cat("Number of groups/populations: ",ntt,"\n")
	}
	#### NML
	if(NML.nloci>0){	
		if(meth) cat("\n\n*****************************\nAnalysis of NML\n")
		DM<-lingoes(as.dist(smc(matN, dist=TRUE))) #smc(scrime), lingoes(ade4)
		if(do.cluster){
			# cluster
			DM_copy <- DM
			attr(DM_copy, "Labels") <- inds
			np <- table(groups)[]
			MSL.cluster <-nj(DM_copy) 
			darksch <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#A6761D","#666666","#E64AB02")
			tipCol<-rep(darksch[1:length(np)],np)
			ecol<-unlist(lapply(MSL.cluster$edge[,2], edgeCol, length(MSL.cluster$tip.label), tipCol))
			png(filename=paste(name,"NML-NJ.png", sep='-'))
			plot.phylo(MSL.cluster, tip.color=rep(darksch[1:length(np)],np), use.edge.length=T, edge.color=ecol, edge.width=3, show.tip.label=T, main="NML")
			legend("bottomright", as.character(levels(groups)), col=darksch, lwd=3)
			dev.off()
		}
		#PCoA
		if(do.pcoa) pcoa(DM, groups, name, "NML")
		assign("DM.MNL", DM, envir=globalenv())
		if(do.amova){
			#AMOVA
			cat("\nPerforming AMOVA\n")
			diffAmova(DM, groups, nDec, do.pairwisePhiST)
		}
	
	}
	
	 
}


edgeCol <- function(x, n, pal){
if(x > n) col <- "black"
else col <-pal[x]
return(col)
}





