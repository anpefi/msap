\name{polymorphic}
\alias{polymorphic}

\title{
	Checks polymorphic state
}
\description{
	Checks if a list of binary alleles is polymorphic (with at least a defined number ocurrences of both states)
}
\usage{
polymorphic(col, threshold)

}
\arguments{
  	\item{col}{
  		A list of binary alleles
	}
	\item{threshold}{
  		An integer defining the minimum threshold of occurrences in each state to define the locus as polymorphic (default 1) 
	}
      
}

\author{
	Andres Perez-Figueroa (\email{anpefi@uvigo.es})
}

