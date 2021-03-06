\name{RV}
\alias{RV}
\alias{bootRV}
\alias{RVrtest}
\alias{RVdiff.rtest}
\alias{multiRV.rtest}
\alias{kRVplot}
\title{
Set of functions related to of RV coefficients
}
\description{
Pairwise Analyse of the similarity among several tables
}
\usage{
RV(df1,df2,center=TRUE)
bootRV(df1, df2, nrepet = 99,center=TRUE){
RVdiff.rtest(df1,df2,df3,df4=NULL, nrepet = 99,center=TRUE)
RVrtest(df1, df2, nrepet = 99,center=TRUE) 
multiRV.rtest(X,nrepet=99,method=NULL,...)
kRVplot(X, whichinrow = NULL, whichincol = NULL, gap = 4, nclass = 10, coeff = 1, nrepet = 99, col = jet.colors(10), digits = 2, ...)
}
\arguments{
  \item{X}{
an object of class 'ktab'.
}
  \item{df1,df2,df3,df4}{
data frames with the same rows.
}
  \item{center}{
a logical value (TRUE by default).

}
  \item{nrepet}{
df1, df2: two data frames with the same rows.
}
  \item{whichinrow}{
a vector of integers to select the graphs in rows (if NULL all the graphs are computed)
}
  \item{whichincol}{
a vector of integers to select the graphs in columns (if NULL all the graphs are computed).
}
  \item{gap}{
an integer to determinate the space between two graphs.
}
  \item{nclass}{
a number of intervals for the histogram.
}
  \item{coeff}{
an integer to fit the magnitude of the graph.
}
  \item{nrepet}{
the number of permutations
}
  \item{col}{
a list of colors such as that generated by 'rainbow','heat.colors', 'topo.colors', 'terrain.colors' or similar functions. 
}
  \item{digits}{
minimal number of _significant_ digits (by default, 2).
}
  \item{\dots}{
further arguments passed to or from other methods
}
}
\details{
}
\note{
}
\seealso{
}
\examples{
require(ade4)
require(MASS)
# k-table simulation
tab0 <- as.data.frame(expand.grid(1:10,1:5))
row.names(tab0)<- paste("s",1:50,sep="")
tab1 <- cbind(tab0, mvrnorm(50,mu=rnorm(5,0,1),Sigma=diag(abs(rnorm(5,0,1)))))
tab2 <- cbind(tab0, mvrnorm(50,mu=rnorm(5,0,1),Sigma=diag(abs(rnorm(5,0,1)))))
tab3 <- cbind(tab0, mvrnorm(50,mu=rnorm(5,0,1),Sigma=diag(abs(rnorm(5,0,1)))))
listab <- list(tab1,tab2,tab3)
names(listab) <- c("tab1","tab2","tab3")
listpca <- list()
for(i in 1:3)
listpca[[i]] <- dudi.pca(listab[[i]],scannf=FALSE,nf=2) 
# construction of the K-table
ktab1 <- ktab.list.dudi(listpca,tabnames= names(listab))
# graphical representation 
kRVplot(ktab1)
multiRV.rtest(ktab1)
RVdiff.rtest(listpca[[1]]$tab,listpca[[2]]$tab,listpca[[3]]$tab)
RV(listpca[[1]]$tab,listpca[[2]]$tab)
bootRV(listpca[[1]]$tab,listpca[[2]]$tab)
}

