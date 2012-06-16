


pcoa <- function(DM, groups,name,surname){
	ntt <- length(levels(groups))
	#PCoA
	pcol <- dudi.pco(DM, scannf = FALSE, nf = 2) #from ade4
	2  #This is needed
	var <- pcol$eig / sum(pcol$eig) *100
	var
	var1 <- round(var[1], digits=1)
	var2 <- round(var[2], digits=1)
	spcoo<-split(pcol$li, groups)
	maxX <- max(pcol$li$A1)
	minX <- min(pcol$li$A1)
	maxY <- max(pcol$li$A2)
	minY <- min(pcol$li$A2)
	fullname <- paste(name,surname, sep='-')
	png(filename=paste(fullname,"png",sep='.'), width=680, height=680)
	par(bty = 'n')
	plot(0,0, main=paste(name,surname, sep=": "), type = "n",xlab=paste("C1 (",var1,"%)"),ylab=paste("C2 (",var2,"%)"), xlim=c(minX,maxX+0.1), ylim=c(minY,maxY+0.1), frame=TRUE, cex=1.5)
	bgcolors<-c("black", "green", "yellow", "red", "blue", "orange", "pink")
	symbs <- c(21,22,23,24,25,26,27) #What if we have > 7 groups???
	for(i in 1:ntt){
		points(spcoo[[i]], pch=symbs[i], col="black", bg=bgcolors[i])
	}	
	s.class(pcol$li, groups, cpoint=0, add.plot=TRUE)
	dev.off()

}
