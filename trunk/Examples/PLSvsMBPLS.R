#### comparison of pls and mbpls for two tables (Y and X)
#### 2012-11-21 by pbady

require(ade4)
require(mixOmics)

source("mbpls.R")
source("mbpls.default.R")

data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological

pls1 <- pls(X, Y, mode = "classic",ncomp=3)
mb1 <- mbpls(list(X),Y,option="uniform",deflation="super")

par(mfrow=c(2,2))
plot(pls1$variates$Y[,1],smb1$U[,1])
abline(0,1,col="red")
plot(pls1$variates$Y[,2],smb1$U[,2])
abline(0,1,col="red")
plot(pls1$variates$Y[,3],smb1$U[,3])
abline(0,1,col="red")

