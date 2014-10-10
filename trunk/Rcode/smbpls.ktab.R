smbpls.ktab <- 
function(X,Y,keepX=NULL,keepY=NULL,nf=3,option = c("inertia", "lambda1", "uniform", "internal"),deflation=c("super","block"),tol=1e-6,max.iter=100,...)
{
	print("In construction!!! problem with weighting!!!")	
	if (!inherits(X, "ktab")) 
		stop("object 'ktab' expected")
	if (!inherits(Y, "dudi")) 
		stop("object 'dudi' expected")
	lw <- X$lw
	if(all(lw != Y$lw))
		stop("non convient arguments!")
	if(is.null(keepX))
		keepX <- matrix(rep(X$blo,nf),ncol=nf,nrow=length(X$blo))
	if(is.null(keepY))
		keepY <- rep(ncol(Y$tab),nf)
	Iner <- function(x) sum(svd(x)$d^2)
	nlig <- length(lw)
	cw <- X$cw
	ncol <- length(cw)
	nbloc <- length(X$blo)
	bloc <- X$blo
	indicablo <- X$TC[, 1]
	Xsepan <- sepan(X, nf = 4)
	rank.fac <- factor(rep(1:nbloc, Xsepan$rank))
	tabw <- NULL
	auxinames <- ktab.util.names(X)
	if (option == "lambda1") {
		for (i in 1:nbloc) tabw <- c(tabw, 1/Xsepan$Eig[rank.fac == 
									i][1])
	}
	else if (option == "inertia") {
		for (i in 1:nbloc) tabw <- c(tabw, 1/sum(Xsepan$Eig[rank.fac == 
											i]))
	}
	else if (option == "uniform") {
		tabw <- rep(1, nbloc)
	}
	else if (option == "internal") 
		tabw <- X$tabw
	else stop("Unknown option")
	# ponderation des tableaux
	for (i in 1:nbloc){
		X[[i]] <- X[[i]] * sqrt(tabw[i])
	}
	# ponderation ligne et colonne X
	InerX <- rep(NA,length(bloc))
	Xtemp <- unclass(X)[1:nbloc]
	for (i in 1:nbloc){
		Xtemp[[i]] <- Xtemp[[i]] * sqrt(X$lw)
		Xtemp[[i]] <- t(t(Xtemp[[i]]) * sqrt(X$cw[X$TC[,1]==i]))
		InerX[i] <- Iner(Xtemp[[i]])
	}	
	# ponderation ligne et colonne Y	
	Ytemp <- Y$tab
	Ytemp  <- Ytemp  * sqrt(Y$lw)
	Ytemp  <- t(t(Ytemp) * sqrt(Y$cw))
	InerY <- Iner(Ytemp)
	# th and u by components
	Tt <- list() 
	Tb <- list() 
	Qt <- list()
	Pb <- list()
	U <- list()
	Wb <- list()
	Wt <- list()
	ExpVarX <- matrix(rep(X$blo,nf),ncol=nf,nrow=length(X$blo))
	ExpVarY <- rep(NA,nf)
	hiter <- NULL
	for(h in 1:nf){
		u <- Ytemp[,1]
		iter <- 1
		tt.old <- 0
		# loop until convergence of tt
		repeat{
			wb <- list()
			tb <- list()			
			for(k in 1:nbloc){
				nx <- ncol(Xtemp[[k]])-keepX[k,h]
				w <- t(Xtemp[[k]])%*% u
				if(nx !=0){
					w = ifelse(abs(w) > abs(w[order(abs(w))][nx]),(abs(w) - abs(w[order(abs(w))][nx])) * sign(w), 0)	
				}
				w <- w/drop(crossprod(u))
				w <- w/sqrt(drop(crossprod(w)))
				wb[[k]] <- w
				tb[[k]] <- as.matrix(Xtemp[[k]])%*% w
			}
			T <- do.call("cbind",tb)
			wt <- t(T)%*% (u/drop(crossprod(u)))
			wt <- wt/sqrt(drop(crossprod(wt)))
			tt <- T%*% (wt/drop(crossprod(wt)))
			q <- t(Ytemp)%*% tt
			ny <- ncol(Ytemp)-keepY[h]	
			if (ny != 0) {
				q = ifelse(abs(q) > abs(q[order(abs(q))][ny]),(abs(q) - abs(q[order(abs(q))][ny])) * sign(q), 0)
			}
			q <- q/drop(crossprod(tt))
			u <- as.matrix(Ytemp)%*% (q/drop(crossprod(q)))
			if (crossprod(tt - tt.old) < tol) 
				break
			if (iter == max.iter) {
				warning(paste("Maximum number of iterations reached for dimension",h), call. = FALSE)
				break
			}
			tt.old <- tt
			iter <- iter + 1
		}
		hiter <- c(hiter,iter)
	# Deflation Wangen and Kowalski (1988)
	if(deflation=="block"){	
			pb <- list()	
			for(k in 1:nbloc){
				pb[[k]] <- t(Xtemp[[k]])%*%(tb[[k]]/drop(crossprod(tb[[k]])))
				Xtemp[[k]] <- as.matrix(Xtemp[[k]])-as.matrix(tb[[k]])%*%t(as.matrix(pb[[k]]))
				ExpVarX[k,h] <- 1-Iner(Xtemp[[k]])/InerX[[k]]
			}
			Ytemp <- as.matrix(Ytemp)-as.matrix(tt)%*%t(q)
			ExpVarY[h] <- (1-Iner(Ytemp))/InerY
		}else if(deflation=="super"){
			pb <- list()	
			for(k in 1:nbloc){
				pb[[k]] <- t(Xtemp[[k]])%*%(tt/drop(crossprod(tt)))
				Xtemp[[k]] <- as.matrix(Xtemp[[k]])-as.matrix(tt)%*%t(as.matrix(pb[[k]]))
				ExpVarX[k,h] <- 1-Iner(Xtemp[[k]])/InerX[[k]]
			}
			Ytemp <- as.matrix(Ytemp)-as.matrix(tt)%*%t(q)
			ExpVarY[h] <- 1-Iner(Ytemp)/InerY
		}else
			stop("non convenient deflation method!")
		# stockage 	
		Tt[[h]] <- tt
		Qt[[h]] <- q
		Pb[[h]] <- unlist(pb)
		Tb[[h]] <- unlist(tb)
		U[[h]] <- u
		Wb[[h]] <- unlist(wb)
		Wt[[h]] <- wt
	}
	Tt <- do.call("cbind",Tt) 
	Qt <- do.call("cbind",Qt)
	U <- do.call("cbind",U)
	Wt <- do.call("cbind",Wt)
	Wb <- do.call("cbind",Wb)
	Pb <- do.call("cbind",Pb)
	Tb <- do.call("cbind",Tb)
	res <- list(Tt,Qt,Wt,Pb,U,Wb,Tb,hiter,X$TC,X$TL)
	names(res)<-c("Tt","Qt","Wt","Pb","U","Wb","Tb","iter","TC","TL")
	res$Y <- Y$tab
	res$X <- X[1:nbloc]
	res$InerY <- InerY
	res$InerX <- InerX
	res$ExpVarX <- ExpVarX
	res$ExpVarY <- ExpVarY
	res$blo <- bloc
	res$Xnames <- auxinames
	res$keepX <- keepX
	res$keepY <- keepY
	res$class <- "smbpls"
	res$call <- match.call()	
	return(res)
} 


# ========================= END ! =================================
