\name{msap}
\alias{msap}

\title{
	Diversity analysis of MSAP data
}
\description{
	It performs all the analysis of MSAP data in order to asses both epigenetic and genetic diversity.   
}
\usage{
	msap(datafile, name=datafile, uninformative=TRUE, nDec=4, meth=TRUE, rm.redundant=TRUE, rm.monomorphic=TRUE, do.pcoa=TRUE, do.shannon=TRUE, do.amova=TRUE, do.pairwisePhiST=TRUE, do.cluster=TRUE, use.groups=NULL, do.mantel=FALSE, np.mantel=1000)
}

\arguments{
  	\item{datafile}{
  		String containing the url of the csv file with the data. Required.
	}
  	\item{name}{
  		a name for the dataset to be included in the output files. By default, the name of the given datafile is used.  
	}
  	\item{uninformative}{
  		A logical value determining how to deal with HPA-/MSP- pattern. 'FALSE' assumes that HPA-/MSP- (no band for both isoschizomers) pattern represents full methylation of cytosines in the target, while 'TRUE' (default value) consider that pattern as uninformative as could be caused by a missing target (mutation). See 'Details' below.
	}
  
  	\item{nDec}{
  		number of digits of precision for floating point output.
	}
	\item{meth}{
  	Logical value switching between MSAP ('TRUE') and standard AFLP ('FALSE') analysis. The difference lies in that for AFLP (meth=FALSE) the 'enzyme' column is ignored and every row in data represent an independent sample, without combination of data. 
	}
	\item{rm.redundant}{
		Not implemented yet.
	}
	\item{rm.monomorphic}{
		Logical value switching between the removal ('TRUE', by default) of monomorphic fragments (defined as those with only one state or just one ocurrence of the second state across the whole dataset) after data transformation.
	}
	\item{do.pcoa}{
		Option switcher for doing a Principal Coordinate Analysis for variation between groups. TRUE by default.
	}
	\item{do.shannon}{
		Option switcher for Shannon's Diversity Index comparison between MSL and NML.
	}
	\item{do.amova}{
		Option switcher for doing an AMOVA for differentiation between groups. TRUE by default.
	}
	\item{do.pairwisePhiST}{
		Logical value switching between the calculation of the pairwise Phi_st between pairs of groups/populations ('TRUE' by default) or skip it ('FALSE').
	}
	\item{do.cluster}{
		Calculates and plots a Neighbour-Joining tree ('TRUE' by default) or skip it ('FALSE').
	}
	\item{use.groups}{
		Gives the groups/populations/treatments of the datafile to be analysed. By default all groups are considered into de the analysis. To provide a subset of the groups a vector should be passed with the names of groups to be included. For example, in a datafile with 5 groups (Control, pop1, pop2, pop3 and pop4) we are interested only in Control and pops 1 and 3. Then, msap should be called with 'use.groups=c('Control','pop1','pop3')'.
	}
	\item{do.mantel}{
		Performs a Mantel test to obtain correlation between MSL and NML ('TRUE') or skip it ('FALSE' by default).
	}
	\item{np.mantel}{
		Gives the number of permutations for the above Mantel test (1000 by default) or skip it ('FALSE').
	}

}
\details{
  	This function is the main interface of the msap package. The only required argument is a string with the name (uri) of the data file to be analysed.
   
  	Data file should be a .csv file with markers as columns and two rows by sample, one for each isoschizomer reaction. The first row should include the markers name/references. The first column should provide the label for the group where the sample is included, with the aim to make comparisons between different gruops. Second column is reserved for an arbitrary label (i.e. to name the sample). Third column should identify the isochizomer with 'HPA' or 'MSP'.
  
}

\author{
	Andres Perez-Figueroa (\email{anpefi@uvigo.es})
}
\seealso{
	\code{\link{rmsap-package}}
}
\examples{
	#Perform all default analysis
	\dontrun{msap("MyDataFile", name="MyAnalysis")}
	#Perform analysis assumes that HPA-/MSP- represents full methylation 
	\dontrun{msap("MyDataFile", name="MyAnalysis", uninformative=FALSE)}
	#The same that above but skipping the PCoA results 
	\dontrun{msap("MyDataFile", name="MyAnalysis", uninformative=FALSE, do.pcoa=FALSE)}
	#Using only some of the populations
	\dontrun{msap("MyDataFile", name="MyAnalysis", use.groups=c('Control','pop1','pop3'))}
}

