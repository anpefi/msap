\name{repMet}
\alias{repMet}

\title{
Report methylation levels
}
\description{
Prints a report of the relative number or different HPA-MSP band patterns in all groups  
}
\usage{
repMet(dataM,groups,nDec)

}
\arguments{
  \item{dataM}{
  A matrix of mixed band patterns, where each row is a different sample. Values of 11, 10, 1 and 0 represent patterns HPA+/MSP+,HPA+/MSP-,HPA-/MSP+,HPA-/MSP-, respectively.
}
   \item{groups}{
   	A list with the label for the group for every sample (row) in 'dataM'
   }
   
   \item{nDec}{
  number of digits of precision for floating point output.
}
}

\author{
Andres Perez-Figueroa (\email{anpefi@uvigo.es})
}

