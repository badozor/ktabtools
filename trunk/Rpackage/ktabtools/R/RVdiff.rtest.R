
#-------------------------------------
# additional function
#-------------------------------------
RVdiff.rtest <- function (df1,df2,df3,df4=NULL, nrepet = 99,center=TRUE){
	if (!is.data.frame(df1)) 
		stop("data.frame expected")
	if (!is.data.frame(df2)) 
		stop("data.frame expected")
	if (!is.data.frame(df3)) 
		stop("data.frame expected")
	l1 <- nrow(df1)
	if (nrow(df2) != l1) 
		stop("Row numbers are different")
	if (nrow(df3) != l1) 
		stop("Row numbers are different")
	if (any(row.names(df2) != row.names(df1))) 
		stop("row names are different")
	if(is.null(df4))
		df4 <- df3
	if (!is.data.frame(df4)) 
		stop("data.frame expected")	
	if (nrow(df4) != l1) 
		stop("Row numbers are different")	
	if(center){    
		X <- scale(df1, scale = FALSE)
		Y <- scale(df2, scale = FALSE)
		Z <- scale(df3, scale = FALSE)
		W <- scale(df4, scale = FALSE)
	}else{
		X <- df1
		Y <- df2
		Z <- df3
		W <- df4
	}
	X <- X/(sum(svd(X)$d^4)^0.25)
	Y <- Y/(sum(svd(Y)$d^4)^0.25)
	Z <- Z/(sum(svd(Z)$d^4)^0.25)
	W <- W/(sum(svd(W)$d^4)^0.25)	
	X <- as.matrix(X)
	Y <- as.matrix(Y)
	Z <- as.matrix(Z)
	W <- as.matrix(W)	
	# absolute difference between the two observed RVs
	obs <- abs(sum(svd(t(X) %*% Z)$d^2)-sum(svd(t(Y) %*% W)$d^2))
	if (nrepet == 0) 
		return(obs)
	perm <- matrix(0, nrow = nrepet, ncol = 1)
	# empricial distribution of the absolute difference between two 'random' RVs
	perm <- apply(perm, 1, function(x) abs(sum(svd(t(X[sample(l1),]) %*% Z[sample(l1),])$d^2)-sum(svd(t(Y[sample(l1),]) %*% W[sample(l1),])$d^2)))
	w <- as.randtest(obs = obs, sim = perm, call = match.call(),alter="greater")
	return(w)
}
