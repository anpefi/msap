\name{msap}
\alias{msap}

\title{
Diversity analysis of MSAP data
}
\description{
It performs all the analysis of MSAP data in order to asses both epigenetic and genetic diversity.   
}
\usage{
msap(datafile, name=datafile, uninformative=TRUE, nDec=4)
}

\arguments{
  \item{datafile}{
  String containing the url of the csv file with the data. Required.
}
  \item{name}{
  a name for the dataset to be included in the output files. By default, the name of the given datafile is used.  
}
  \item{uninformative}{
  A logical value determining how to deal with HPA-/MSP- pattern. 'FALSE' assumes that HPA-/MSP- (no band for both isoschizomers) pattern represents full methylation of cytosines in the target, while 'TRUE' (default value) consider that pattern as uninformative as could be caused by a missing target (mutation). See 'Details' below
}
  
  \item{nDec}{
  number of digits of precision for floating point output.
}

}
\details{
  This function is the main interface of the msap package. The only required argument is a string with the name (uri) of the data file to be analysed.
   
  Data file should be a .csv file with markers as columns and two rows by sample, one for each isoschizomer reaction. The first row should include the markers name/references. The first column should provide the label for the group where the sample is included, with the aim to make comparisons between different gruops. Second column is reserved for an arbitrary label (i.e. to name the sample). Third column should identify the isochizomer with 'HPA' or 'MSP'



  
}

\author{
Andres Perez-Figueroa (\email{anpefi@uvigo.es})
}
\seealso{
\code{\link{rmsap-package}}
}

