#===========================================================
# test global de concordance entre plusieurs tableaux
# d'un K-tableau (avec operator, cf STATIS)
# P.BADY                25/11/2003 17:30:44
#===========================================================

multiRV.rtest <- function(X,nrepet=99,method=NULL,...){
# based on the ADE-4 software
# verification
	if (!inherits(X, "ktab")) stop("object 'ktab' expected")
	if (is.null(method)) {
		cat("1 = Sum of RV-Coefficients\n")
		cat("2 = First Eigenvalue of interstructure\n")
		cat("3 = Squared Sum of interstructure Eigenvalues\n")
		cat("Select an integer (1-3): ")
		method <- as.integer(readLines(n = 1))
	}
	# Data preparartion
	lw <- X$lw
	nlig <- length(lw)
	cw <- X$cw
	ncol <- length(cw)
	blo <- X$blo
	ntab <- length(X$blo)
	indicablo <- X$TC[, 1]
	Xk <- list()
	lwsqrt <- sqrt(lw)
	for (k in 1:ntab) {
		Xk[[k]] <- as.matrix( X[[k]] )
	}
# functions 
	RVs <- function(Xk,cw,lwsqrt,ntab,indicablo){
		sep <- list()
		for (k in 1:ntab) {
			ak <- sqrt(cw[indicablo == k])
			wk <- as.matrix(Xk[[k]]) * lwsqrt
			wk <- t(t(wk) * ak)
			wk <- wk %*% t(wk)
			sep[[k]] <- wk
		}
		sep <- matrix(unlist(sep), nlig * nlig, ntab)
		RV <- t(sep) %*% sep
		ak <- sqrt(diag(RV))
		RV <- sweep(RV, 1, ak, "/")
		RV <- sweep(RV, 2, ak, "/")
		return(RV)
	}
	fun1 <-function(RV){
		# i!=j
		w <- RV - diag(diag(RV),ntab,ntab)
		w <- sum(w)/2
		return(w)
	}
	fun2 <-function(RV){ 
		w <- eigen(RV, sym = TRUE)$values[1]
		return(w)
	}
	fun3 <-function(RV){ 
		w <- sum(eigen(RV, sym = TRUE)$values^2)
		return(w)
	}
# method
	if (any(method == 1)){
		fun=fun1
	}else if(any(method==2)){
		fun=fun2
	}else if(any(method == 3)){
		fun=fun3
	}else stop("Non convenient method")
# computation
	obs <- fun(RVs(Xk,cw,lwsqrt,ntab,indicablo))
	if (any(nrepet == 0)) return(obs)
# permutations
	perm <- matrix(0, nrow = nrepet, ncol = 1)
	for (j in 1:nrepet){
		pXk <-list()
		for(i in 1:ntab){
			vec <-sample(1:nlig)
			pXk[[i]] <- Xk[[i]][vec,]
			plwsqrt <- lwsqrt[vec]
		}
		perm[j] <- fun(RVs(pXk,cw,plwsqrt,ntab,indicablo))
	}
# results
	w <- as.rtest(obs = obs, sim = perm, call = match.call() )
	return(w)
}
