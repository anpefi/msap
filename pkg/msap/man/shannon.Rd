\name{shannon}
\alias{shannon}

\title{
Shannon Index of phenotypic diversity
}
\description{
Calculate Shannon Index of phenotypic diversity given a list of binary alleles.   
}
\usage{
shannon(x)

}
\arguments{
  \item{x}{
  A list of binary values or alleles. 
}
  
}
\details{
  
The Shannon index of phenotypic diversity, S, derived from the Shannon-Weaver index (Shannon, 1948):
\deqn{ S=-\sum_{i=1}^{n} p_{i}\log_{2}p_i }{S=sum(p(i)log2(p(i)))}

where \eqn{p_{i}} is the frequency of the band presence at the ith marker within the population. This index gives more
weight to the presence than to the absence of bands. This has no real biological support, although it might account for the occurrence of homoplasic absences of bands (Bonin et al., 2007)


  
}
\value{
	The value of the Shannon Index of phenotypic diversity
}


\references{
Bonin, A., Ehrich, D. and Manel S. 2007. Statistical analysis of amplified fragment length polymorphism data: a toolbox for molecular ecologists and evolucionists. \emph{Molecular Ecology} 16, 3737-3758.

Shannon, C. E. 1948. A mathematical theory of comunication. \emph{Bell System Technical Journal} 27, 379-423.
}
\author{
Andres Perez-Figueroa (\email{anpefi@uvigo.es})
}

