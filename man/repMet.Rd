\name{repMet}
\alias{repMet}

\title{
        Report methylation levels
}
\description{
        Prints a report of the relative number or different HPA-MSP band patterns in all groups  
}
\usage{
        repMet(dataM,groups,nDec,enz1="HPA",enz2="MSP", pattern=c(1,2,2,NA))

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
	\item{pattern}{
	        vector of methylation states (1 - unmethylated, 2- methylated, NA - missing data) for the four combinations of both enzymes enz1+/enz2+, enz1+/enz2-, enz1-/enz2+, enz1-/enz2-. By default (1,2,2,NA) corresponding with standard MSAP using HPA/MSP considering -/- as uninformative status. Note that the value of this argument could be modified if using the deprecated arguments no.bands or uninformative.
	}
}

\author{
Andres Perez-Figueroa (\email{anpefi@uvigo.es})
}

