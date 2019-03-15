repMet <-function(dataM,groups,nDec, enz1, enz2, pattern=c(1,2,2,NA)){
	cat("Report of methylation levels \n")
	dataM[dataM==11]<-"u"
	dataM[dataM==10]<-"h"
	dataM[dataM==1]<-"i"
	dataM[dataM==0]<-"f"
	dataM <- split(dataM,groups)
	
	if(identical(pattern,c(1,2,2,NA)) && enz1=="HPA" && enz2=="MSP") #Standard MSAP HPA/MSP
	        metLabels <- c("(Unmethylated)","(Hemimethylated)","(Internal C methylation)","(Uninformative)")
	else{
	        metLabels <- c("(Unmethylated)","(Methylated)")[pattern]
	        metLabels[is.na(metLabels)] <- "(Uninformative)"
	}
	
	ns<-names(dataM)
	res<-matrix(ncol=length(ns), nrow=4)
	colnames(res)<-ns
	rownames(res)<-c(
	        paste0(enz1,"+/",enz2,"+ ",metLabels[1]),
	        paste0(enz1,"+/",enz2,"- ", metLabels[2]),
	        paste0(enz1,"-/",enz2,"+ ", metLabels[3]),
	        paste0(enz1,"-/",enz2,"- ", metLabels[4])
	        )
	for(x in seq_along(dataM)){
	res[1,x]<- as.numeric(table(as.matrix(dataM[[x]]))["u"])  #Can I improve this?
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
  return(dataM)
}