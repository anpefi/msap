\name{pcoa}
\alias{pcoa}

\title{
Principal Coordinate Analyses
}
\description{
This function uses ade4's dudi.pco and s.class to obtain a PCoA and represent its results in a 2-D figure. It also saves a .csv file with the coordinates for all individuals and another .csv file with the eigenvalues for all axis.
}
\usage{
pcoa(DM, groups, inds, name, surname)

}
\arguments{
  \item{DM}{
  A Euclidean Distance Matrix
}
   \item{groups}{
   	A factor with the label for the group for every sample (row) in 'dataM'
   }
   \item{inds}{
   	Label list for individuals in 'dataM'
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

