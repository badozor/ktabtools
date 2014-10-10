### extract of RV.rtest from R package ade4
RV <- function(df1,df2,center=TRUE){
	if(nrow(df1) != nrow(df2))
		stop("non convenient dimension!")
	X <- as.matrix(df1)
	Y <- as.matrix(df2)
	if(center){    
		X <- scale(X, scale = FALSE)
		Y <- scale(Y, scale = FALSE)
	}else{
		X <- X
		Y <- Y
	}
	X <- X/(sum(svd(X)$d^4)^0.25)
	Y <- Y/(sum(svd(Y)$d^4)^0.25)
	X <- as.matrix(X)
	Y <- as.matrix(Y)
	sum(svd(t(X) %*% Y)$d^2)			
}

