\name{RV.mcoa}
\alias{RV.mcoa}
\title{
Computation of RV coefficient based on object 'mcoa'
}
\description{
The function computes RV coefficients between reference (Synthetic scores) and individual 
representation of each table or RV coefficient among all the individual 
representation of each table.
}
\usage{
RV.mcoa(m,nfs=1:m$nf,option=c("reference"),distance=FALSE,...)
}
\arguments{
  \item{m}{
an object of class 'mcoa'.
}
  \item{nfs}{
a numeric vector indicating the kept axes
}
  \item{otpion}{
a character for the type of RV computation. If option is equal to "reference", the function computes RV 
coefficients between reference (Synthetic scores) and individual representation of each table (by default). 
If option is equal to "all", the function computes RV coefficient among all the individual representation of each table
}
  \item{distance}{
a logical value indicating if the RV distances are returned.
}
  \item{\dots}{
further arguments passed to or from other methods
}
}
\seealso{
mcoa
}
\references{
Chessel, D. and Hanafi, M. (1996) Analyses de la co-inertie de K nuages de points, _Revue de Statistique Appliquée_, *44*, 35-60.\cr
Escoufier, Y. 1973. Le traitement des variables vectorielles. _Biometrics_, *29*, 750-760.\cr
}
\examples{
require(ade4)
data(friday87)
w1 <- data.frame(scale(friday87$fau, scal = FALSE))
w2 <- ktab.data.frame(w1, friday87$fau.blo, tabnames = friday87$tab.names)
mcoa1 <- mcoa(w2, "lambda1", scan = FALSE)
mcoa1
# RV for individual representation and reference
RV.mcoa(mcoa1)
# RV among individual representation
RV.mcoa(mcoa1,option="all")
# table clustering
drv1 <- RV.mcoa(mcoa1,option="all",distance=TRUE)
hc1 <- hclust(drv1,method="ward.D")
plot(hc1,hang=-1)
}

