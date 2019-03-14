\name{repMet}
\alias{repMet}

\title{
        Report methylation levels
}
\description{
        Prints a report of the relative number or different HPA-MSP band patterns in all groups  
}
\usage{
        repMet(dataM,groups,nDec,enz1="HPA",enz2="MSP")

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
        \item{enz1}{
		String giving the label used in the datafile for the first enzyme. By default, as it considers standard MSAP with HPA/MSP, the value is "HPA". Note that this is case-sensitive and it should match the label given in the datafile.
	}
        \item{enz2}{
		String giving the label used in the datafile for the second enzyme. By default, as it considers standard MSAP with HPA/MSP, the value is "MSP". Note that this is case-sensitive and it should match the label given in the datafile.
	}
}

\author{
Andres Perez-Figueroa (\email{anpefi@uvigo.es})
}

