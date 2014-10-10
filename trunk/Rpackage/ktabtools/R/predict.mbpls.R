predict.mbpls <-
function (x,data=NULL,nf=NULL,...) 
{
print("function in test!")
print("problem, if we use weighting of table and new datasets")
print("ok for uniform ponderation!")
	Iner <- function(x) sum(svd(x)$d^2)
	eigen1 <- function(x) svd(x)$d[1]^2
	option <- x$option
	deflation <- x$deflation
	bloc <- x$blo
	nbloc <- length(bloc)
	indicablo <- x$TC[,1]
	if(is.null(nf))
		nf <- ncol(x$U)
	if (is.null(data)){
	        X <- x$X
	}else{
	if (!inherits(data, "list")) 
	        stop("object 'list' expected")
    	if (any(unlist(lapply(data, function(x) !inherits(x, "data.frame"))))) 
        	stop("list of 'data.frame' object expected")
	if(length(data)!=length(bloc))
		stop("stop non convenient dimension for data")
    	InerX <- Eig <- rep(NA, nbloc)
### problème si ponderation du tableau different de 1 pour un nouveau jeu de donnees	
	for (i in 1:nbloc){
		mi <- attributes(x$X[[i]])$"scaled:center"
		sdi <- attributes(x$X[[i]])$"scaled:scale"
        	X[[i]] <- scale(data[[i]],mi,sdi)
        	w <- svd(X[[i]])$d^2
        	Eig[i] <- w[1]
        	InerX[i] <- sum(w)
    	}
    	if (option == "lambda1") {
        	for (i in 1:nbloc) tabw <- c(tabw, 1/Eig[i])
    	}else if (option == "inertia") {
        	for (i in 1:nbloc) tabw <- c(tabw, 1/InerX[i])
    	}else if (option == "uniform"){
        	tabw <- rep(1, nbloc)
    	}else stop("Unknown option")
    		for (i in 1:nbloc) {
        	X[[i]] <- X[[i]] * sqrt(tabw[i])
   	}
	}
	util1 <- klistutil(X)
	Y <- matrix(0,ncol=ncol(x$Y),nrow=nrow(X[[1]]))
	Xnew <- list()
	for (k in 1:nbloc) {
		Xnew[[k]] <- matrix(0,ncol=ncol(X[[k]]),nrow=nrow(X[[k]]))
	}
    Tt <- list()
    Tb <- list()
    Pb <- list()
    Wt <- list()
	wb <- x$Wb
for(h in 1:nf){		    
    	u <- x$U[,h]
	q <- x$Qt[,h]
	tb <- list()
	pb <- list()
        for (k in 1:nbloc) {
         	tb[[k]] <- as.matrix(X[[k]]) %*% wb[indicablo==k,h]
        }
            T <- do.call("cbind", tb)
            wt <- t(T) %*% (u/drop(crossprod(u)))
            wt <- wt/sqrt(drop(crossprod(wt)))
            tt <- T %*% (wt/drop(crossprod(wt)))
        if (deflation == "block") {
            pb <- list()
            for (k in 1:nbloc) {
                pb[[k]] <- t(X[[k]]) %*% (tb[[k]]/drop(crossprod(tb[[k]])))
		auxi <- as.matrix(tb[[k]]) %*%t(as.matrix(pb[[k]]))
                X[[k]] <- as.matrix(X[[k]]) - auxi
		Xnew[[k]] <- Xnew[[k]] + auxi
            }
            Y <- Y + as.matrix(tt) %*% t(q)
        }
        else if (deflation == "super") {
            pb <- list()
            for (k in 1:nbloc) {
                pb[[k]] <- t(X[[k]]) %*% (tt/drop(crossprod(tt)))
		auxi <- as.matrix(tt) %*%t(as.matrix(pb[[k]]))
          	X[[k]] <- as.matrix(X[[k]]) - auxi
		Xnew[[k]] <- Xnew[[k]] + auxi
            }
            Y <- Y + as.matrix(tt) %*% t(q)
        }
        else stop("non convenient deflation method!")
        Tt[[h]] <- tt
        Pb[[h]] <- unlist(pb)
        Tb[[h]] <- unlist(tb)
        Wt[[h]] <- wt
	}
    Tt <- do.call("cbind", Tt)
    Tb <- do.call("cbind", Tb)
	colnames(Tt) <- colnames(Tb) <- paste("CS",1:nf,sep="")
	rownames(Tt) <- rownames(Y)
	rownames(Tb) <- util1$rownames
    Wt <- do.call("cbind", Wt)
	colnames(Wt) <- paste("W",1:nf,sep="")
	rownames(Wt) <- util1$tabnames
    Pb <- do.call("cbind", Pb)
	colnames(Pb) <- paste("LD",1:nf,sep="")
	rownames(Pb) <- util1$colnames
	return(list(Y=Y,X=Xnew,Tt=tt,Tb=Tb,Pb=Pb,Wt=Wt,nf))
}
