\name{mbpls}
\alias{mbpls}
\alias{mbpls.default}
\alias{mbpls.ktab}
\alias{print.mbpls}
\title{
Multi-Block Partial Least Squares (MBPLS)
}
\description{
Performs a multi-block partail least squares, using list of K data.frames (Xb) associated with a data.frame (Y) 
or object 'ktab' (Xb) associated with object 'dudi' (Y).
}
\usage{
mbpls(X, ...)
\method{mbpls}{default}(X, Y, option = c("inertia", "lambda1", "uniform"), deflation = c("super", "block"), 
	nf = 3, tol = 1e-06, max.iter = 100, center.x=TRUE,scale.x=TRUE,center.y=TRUE,scale.y=TRUE,...)
\method{mbpls}{ktab}(X, Y, option = c("inertia", "lambda1", "uniform", "internal"), deflation = c("super", "block"), 
	nf = 3, tol = 1e-06, max.iter = 100, ...)
}
\arguments{
  \item{X}{
a list of data.frame or object 'ktab'
}
  \item{Y}{
a data.frame or object 'dudi'
}
  \item{option}{
a character describing the weighting of the table (inertia, lambda1 or uniform, see details)
}
  \item{deflation}{
a character describing deflation method ("super" or "block", see details)
}
  \item{nf}{
number of components
}
  \item{tol}{
a positive real, the tolerance used in the iterative algorithm.
}
  \item{max.iter}{
integer, the maximum number of iterations.
}
  \item{center.x}{
either a logical value or a numeric vector of length equal to the number of columns of ‘X’
}
  \item{scale.x}{
either a logical value or a numeric vector of length equal to the number of columns of ‘X’.
}
  \item{center.y}{
either a logical value or a numeric vector of length equal to the number of columns of ‘Y’
}
  \item{scale.y}{
either a logical value or a numeric vector of length equal to the number of columns of ‘Y’.
}
  \item{\dots}{
further arguments passed to or from other methods
}
}
\details{
Two methods of deflation are proposed in this function. The first (deflation="super") is the deflation mode proposed by Westerhuis and Coenegracht (1997) 
using the super score 'Tt'for the deflation step. The super score are orthogonal and the block scores 'Tb' are slightly correlated (more details in Westerhuis et al. 1998). 
The second method (deflation="block") is suggested by Wangen and Kowalski and it used the block scores 'Tb' for the deflation of the data block Xb. This deflation strategy 
forces the blocks scores to be orthogonal ('Tb') and the super scores 'Tt' becomes correlated.
}
\note{
The version based on 'ktab' and 'dudi' object is in development and not stabilized!!!
}
\value{
mpls and all the functions that use it return a list with the following components: 
\item{Tt}{an object 'matrix' containing the super scores}
\item{Qt}{an object 'matrix' containing the weights of the variables in Y}
\item{Wt}{an object 'matrix' containing the super weights} 
\item{Pb}{an object 'matrix' containing the loadings of variables in block Xb}
\item{U}{an object 'matrix' containing score of Y}
\item{Wb}{an object 'matrix' containing the weights of the variables in block Xb}    
\item{Tb}{an object 'matrix' containing all the block score of Xb}
\item{iter}{the number of iteration}
\item{Y}{an object 'data.frame' corresponding to the response data} 
\item{X}{a list of 'data.frame' corresponding to descriptor data (all Xb)}
\item{InerY}{Inertia of Y}    
\item{InerX}{Inertia of the block Xb}
\item{ExpVarX}{Explained variances (inertia) in block Xb by component} 
\item{ExpVarY}{Explained variances (inertia) in Y by component}    
\item{blo}{a vector of variable number by blocks}
\item{TC}{a numeric vector with the factors for columns}
\item{TL}{a numeric vector with the factors for rows}
\item{call}{original call} 
}
\references{
González I., Lé Cao, K-A. and Déjean S. (2011) mixOmics: Omics Data Integration Project. URL: http://www.math.univ-toulouse.fr/~biostat/mixOmics/
Lé Cao, K-A., González I. and Déjean S. (2009) intergrOmics: an R package to unravel relationships between two omics data sets. Bioinformatics. 25(21):2855-2856.
Li W, Zhang S, Liu CC, Zhou XJ.(2012)Identifying multi-layer gene regulatory modules from multi-dimensional genomic data.Bioinformatics, 28(19):2458-66.
Westerhuis et al. (1998) Analysis of multiblock and hierarchical pca and pls models, journal of chemometrics, 12, 301-321.
}
\seealso{
}
\examples{
# example1: Comparison of pls and mbpls for two tables (Y and X)
require(ade4)
require(mixOmics)
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
pls1 <- pls(X, Y, mode = "classic",ncomp=3)
mb1 <- mbpls(list(X),Y,option="uniform",deflation="super")
par(mfrow=c(2,2))
plot(pls1$variates$Y[,1],mb1$U[,1])
abline(0,1,col="red")
plot(pls1$variates$Y[,2],mb1$U[,2])
abline(0,1,col="red")
plot(pls1$variates$Y[,3],mb1$U[,3])
abline(0,1,col="red")
}

