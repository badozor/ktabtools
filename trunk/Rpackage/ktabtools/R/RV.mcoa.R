RV.mcoa <- function(m,nfs=1:m$nf,option=c("reference"),distance=FALSE,...){
	if (!inherits(m, "mcoa"))
		stop("non convenient data")
	blo <- sort(unique(m$TL[, 1]))
	nblo <- length(blo)
	res <- NULL
	if(option=="reference"){   
		for(i in 1:nblo){
			X <- scale(m$SynVar[,nfs],...)
			Y <- scale(m$Tl1[m$TL[,1]==levels(m$TL[,1])[i],nfs],...)
			X <- X/(sum(svd(X)$d^4)^0.25)
			Y <- Y/(sum(svd(Y)$d^4)^0.25)
			X <- as.matrix(X)
			Y <- as.matrix(Y)
			w <- sum(svd(t(X) %*% Y)$d^2)
			res <- c(res,w)
		}
		names(res)<- row.names(m$cov2)
	}else if(option=="all"){
		res <- matrix(NA,nc=nblo,nr=nblo)
		for(i in 2:nblo)
			for(j in 1:(i-1)){
			x <- scale(as.matrix(m$Tl1[m$TL[,1]==levels(m$TL[,1])[i],nfs]),...)
			y <- scale(as.matrix(m$Tl1[m$TL[,1]==levels(m$TL[,1])[j],nfs]),...)
			x <- x/(sum(svd(x)$d^4)^0.25)
			y <- y/(sum(svd(y)$d^4)^0.25)
			x <- as.matrix(x)
			y <- as.matrix(y)
			res[j,i] <- res[i,j] <- sum(svd(t(x) %*% y)$d^2)
			}
		diag(res) <- rep(1,nrow(res))
		rownames(res) <- colnames(res) <- rownames(m$cov2)
	}else stop("non convenient 'option'!")
	if(distance){
		res <- as.dist(sqrt(2*(1-res)))
	}
	return(res)
}
