\name{pcoa}
\alias{pcoa}

\title{
Principal Coordinate Analyses
}
\description{
This function uses ade4's dudi.pco and s.class to obtain a PCoA and represent its results in a 2-D figure.
}
\usage{
pcoa(DM, groups,name,surname)

}
\arguments{
  \item{DM}{
  A Euclidean Distance Matrix
}
   \item{groups}{
   	A factor with the label for the group for every sample (row) in 'dataM'
   }
   
  
\item{name}{
  Prefix for filename and title
}

\item{surname}{
  "MNL" or "MSL"
}
}

\author{
Andres Perez-Figueroa (\email{anpefi@uvigo.es})
}

