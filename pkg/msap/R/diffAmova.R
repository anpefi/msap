#diffAmova.R
#package: msap
#Author: Andrés Pérez-Figueroa (anpefi@uvigo.es)

#This function should reduce redundant markers, i.e. markers that shown the same score across samples

diffAmova <- function(DM, groups, nDec){
	
	ntt <- length(levels(groups))

	assign("DM", DM, envir=globalenv())
	assign("groups", groups, envir=globalenv())

	amv<-amova(DM ~ groups, nperm=10000) #from pegas

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
