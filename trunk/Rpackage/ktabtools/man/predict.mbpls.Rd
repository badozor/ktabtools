\name{predict.mbpls}
\alias{predict.mbpls}
\alias{predict.smbpls}
\title{
Prediction for object (s)mbpls
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.mbpls(x,data=NULL,nf=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
}
}
\note{
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
pred1 <- predict(mb1)
pred2 <- predict(mb1,list(X))
}
