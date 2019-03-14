methStatusEval <- function(x, error=0.05, pattern=c(1,2,2,NA)){
	#Evaluate for classification as 'methylation-susceptible locus' or 'non methylated"
	#Uses mixed data: 11, 10, 1, 0 are the four possible states for both enzymes
	threshold= ((error + error) - (2*error*error))
	#A locus is classified as "methylation-susceptible" if the frecuency of methylation evidence (2 in pattern) is above threshold
	
	#define the patterns methylated (2)
	counts <- as.numeric(table(factor(x,levels = c(11,10,1,0)) ) )
	fmet <- sum(counts[pattern==2],na.rm = T) / sum(counts[!is.na(pattern)],na.rm = T)
	
	return(fmet > threshold)
}
