#--------------------------------------------------
# P.BADY                25/11/2003 17:30:44
# + caluls RV pour l'ACOM => 14/02/2005 12:26:07
#--------------------------------------------------
Iner <- function(x) sum(svd(x)$d^2)
eigen1 <- function(x) svd(x)$d[1]^2
kdist2list <- function(kd,...){
	if (!inherits(kd, "kdist")) 
		stop("object 'kdist' expected")
	size <- attributes(kd)$size
	labels <- attributes(kd)$labels
	k <-length(attributes(kd)$names)
	res <- lapply(kd[1:k],function(x)vec2dist(m=x,size=size,labels=labels))
	return(res)
}
kobs2mat <- function(kd,labels=NULL,...){
	if (!inherits(kd, "corkdist")) 
		stop("object 'corkdist' expected")
	size <- max(attributes(kd)$design)
	res <- vec2dist(unlist(lapply(kd,function(x) x$obs)),size=size)
	res <- dist2mat(res)
	diag(res) <- rep(1,size)
	if(!is.null(labels)){
		attributes(res)$dimnames[[1]] <- labels
		attributes(res)$dimnames[[2]] <- labels
	}
	return(res)
}
kpval2mat <- function(kd,labels=NULL,...){
	if (!inherits(kd, "corkdist")) 
		stop("object 'corkdist' expected")
	size <- max(attributes(kd)$design)
	res <- vec2dist(unlist(lapply(kd,function(x) x$pvalue)),size=size)
	res <- dist2mat(res)
	diag(res) <- rep(1,size)
	if(!is.null(labels)){
		attributes(res)$dimnames[[1]] <- labels
		attributes(res)$dimnames[[2]] <- labels
	}
	return(res)
}
diag.RV <- function(RV,...){
	if(!inherits(RV,"data.frame"))
		stop("non convient data")
	dnames <- names(RV)
	res <- list()
	eig1 <- eigen(RV, sym = TRUE)
	res$RV.eig <- eig1$values
	if (any(eig1$vectors[, 1] < 0)) 
		eig1$vectors[, 1] <- -eig1$vectors[, 1]
	tabw <- eig1$vectors[, 1]
	res$RVw <- tabw
	w <- t(t(eig1$vectors) * sqrt(eig1$values))
	w <- as.data.frame(w)
	row.names(w) <- dnames
	names(w) <- paste("S", 1:ncol(w), sep = "")
	res$RV.coo <- w[, 1:min(4, ncol(w))]
	return(res)
}
### to do: add RV between inital tables and SynVar
RV.mcoa <- function(m,...){
	if (!inherits(m, "mcoa")) 
		stop("non convenient data")
	blo <- sort(unique(m$TL[, 1]))
	nblo <- length(blo)
	res <- NULL
	for(i in 1:nblo){
		X <- scale(m$SynVar, scale = FALSE)
		Y <- scale(m$Tl1[m$TL[,1]==i,], scale = FALSE)
		X <- X/(sum(svd(X)$d^4)^0.25)
		Y <- Y/(sum(svd(Y)$d^4)^0.25)
		X <- as.matrix(X)
		Y <- as.matrix(Y)
		w <- sum(svd(t(X) %*% Y)$d^2)
		res <- c(res,w)
	}
	names(res)<- row.names(m$cov2)
	return(res)
}
### 2013-01-10 by pbady
RV <- function(X,Y,center=TRUE){
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
### function in test (not fixed!!!)
RVselect <- function(df1,r=ncol(df1), ...){
	listnames <- list()
	listRV <- list()
	auxi <- apply(df1,2,function(x) RV(df1,as.data.frame(x)))	
	listRV[[1]] <- res0[auxi==max(auxi)][1]
	listnames[[1]] <- names(listRV[[1]])
	id <- which(colnames(df1)==listnames[[1]])
	for(i in 2:r){
		auxi <- apply(df1,2,function(x) RV(df1,as.data.frame(cbind(x,df1[,id]))))
		listRV[[i]] <- auxi[auxi==max(auxi)][1]
		listnames[[i]] <- c(listnames[[i-1]],names(listRV[[i]]))
		id <- c(id,which(colnames(df1)==names(listRV[[i]])))
	}
	return(list(RVs=unlist(listRV),varnames=listnames,id=id))
}
