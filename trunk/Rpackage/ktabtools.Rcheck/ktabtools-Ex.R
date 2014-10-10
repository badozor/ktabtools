pkgname <- "ktabtools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ktabtools')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("kRVplot")
### * kRVplot

flush(stderr()); flush(stdout())

### Name: kRVplot
### Title: Multiple plot representation of RV coefficients
### Aliases: kRVplot

### ** Examples

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



cleanEx()
nameEx("mbpls")
### * mbpls

flush(stderr()); flush(stdout())

### Name: mbpls
### Title: Multi-Block Partial Least Squares (MBPLS)
### Aliases: mbpls mbpls.default mbpls.ktab print.mbpls

### ** Examples

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



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("predict.mbpls")
### * predict.mbpls

flush(stderr()); flush(stdout())

### Name: predict.mbpls
### Title: Prediction for object (s)mbpls
### Aliases: predict.mbpls predict.smbpls

### ** Examples

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



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
