#polymorphic.R
# Andrés Pérez-Figueroa

# From a binary matrix returns a list with indices of those polymorphic columns.
# Polymorphic threshold is now user defined (By default at least two ocurrences for every state)
polymorphic <- function(col, threshold=1){
	tcol <- table(col)
	if(length(tcol)>2){
		cat("ERROR! data matrix in 'polymorphic' is not binary\n")
	} 
	if(length(tcol)==1) return(FALSE)
	if(length(tcol)==2){
		if(min(tcol)>threshold-1) return (TRUE)
		else return (FALSE)
	}
}
