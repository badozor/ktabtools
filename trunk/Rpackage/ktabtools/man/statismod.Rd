\name{statismod}
\alias{statismod}
\title{
Statis
}
\description{
function statis from R package ade4 with addition of the computation of the trajectories for Individuals. 
}
\usage{
statismod(X, scannf = TRUE, nf = 3, tol = 1e-07)
}
\arguments{
  \item{X}{
an object of class 'ktab'
}
  \item{scannf}{
a logical value indicating whether the number of kept axes for the compromise should be asked.
}
  \item{nf}{
if ‘scannf’ FALSE, an integer indicating the number of kept axes for the compromise.
}
  \item{tol}{
a tolerance threshold to test whether the distance matrix is Euclidean : an eigenvalue is considered positive if it is larger than ‘-tol*lambda1’ where ‘lambda1’ is the largest eigenvalue.
}
}
\value{
statismod returns a object 'statis' with trajectories of Lavit and their modification (see more details in the documentation of ADE-4):
  \item{La}{
Lavit's trajectories
}
  \item{Lamod}{
modified Lavit's trajectories
}
}
\references{
}
\seealso{
}
\examples{
}
