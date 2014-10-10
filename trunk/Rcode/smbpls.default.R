smbpls.default <-
function(X, Y, keepX = NULL, keepY = NULL, nf = 3, option = c("inertia","lambda1", "uniform", "internal"), deflation = c("super","block"), tol = 1e-06, max.iter = 100, ...) 
{
    if (!inherits(Y, "data.frame")) 
        stop("object 'ktab' expected")
    if (!inherits(X, "list")) 
        stop("object 'list' expected")
    if (any(unlist(lapply(X, function(x) !inherits(x, "data.frame"))))) 
        stop("list of 'data.frame' object expected")
	util1 <- klistutil(X)
    Iner <- function(x) sum(svd(x)$d^2)
    eigen1 <- function(x) svd(x)$d[1]^2
    nlig <- nrow(Y)
    bloc <- util1$bloc
    nbloc <- length(bloc)
    indicablo <- util1$TC[,1]
    if (is.null(keepX)) 
        keepX <- matrix(rep(bloc, nf), ncol = nf, nrow = length(bloc))
    if (is.null(keepY)) 
        keepY <- rep(ncol(Y), nf)
    InerX <- Eig <- rep(NA, nbloc)
    for (i in 1:nbloc) {
        X[[i]] <- scale(X[[i]])
        w <- svd(X[[i]])$d^2
        Eig[i] <- w[1]
        InerX[i] <- sum(w)
    }
    if (option == "lambda1") {
        for (i in 1:nbloc) tabw <- c(tabw, 1/Eig[1])
    }
    else if (option == "inertia") {
        for (i in 1:nbloc) tabw <- c(tabw, 1/InerX[i])
    }
    else if (option == "uniform") {
        tabw <- rep(1, nbloc)
    }
    else stop("Unknown option")
    for (i in 1:nbloc) {
        X[[i]] <- X[[i]] * sqrt(tabw[i])
    }
    InerX <- rep(NA, nbloc)
    Xtemp <- X
    for (i in 1:nbloc) {
        InerX[i] <- Iner(Xtemp[[i]])
    }
    Y <- as.matrix(Y)
    Y <- scale(Y)
    Ytemp <- Y
    InerY <- Iner(Ytemp)
    Tt <- list()
    Tb <- list()
    Qt <- list()
    Pb <- list()
    U <- list()
    Wb <- list()
    Wt <- list()
    ExpVarX <- matrix(rep(bloc, nf), ncol = nf, nrow = length(bloc))
    ExpVarY <- rep(NA, nf)
    hiter <- NULL
    for (h in 1:nf) {
        u <- Ytemp[, 1]
        iter <- 1
        tt.old <- 0
        repeat {
            wb <- list()
            tb <- list()
            for (k in 1:nbloc) {
                nx <- ncol(Xtemp[[k]]) - keepX[k, h]
                w <- t(Xtemp[[k]]) %*% u
                if (nx != 0) {
                  w = ifelse(abs(w) > abs(w[order(abs(w))][nx]), 
                    (abs(w) - abs(w[order(abs(w))][nx])) * sign(w), 
                    0)
                }
                w <- w/drop(crossprod(u))
                w <- w/sqrt(drop(crossprod(w)))
                wb[[k]] <- w
                tb[[k]] <- as.matrix(Xtemp[[k]]) %*% w
            }
            T <- do.call("cbind", tb)
            wt <- t(T) %*% (u/drop(crossprod(u)))
            wt <- wt/sqrt(drop(crossprod(wt)))
            tt <- T %*% (wt/drop(crossprod(wt)))
            q <- t(Ytemp) %*% tt
            ny <- ncol(Ytemp) - keepY[h]
            if (ny != 0) {
                q = ifelse(abs(q) > abs(q[order(abs(q))][ny]), 
                  (abs(q) - abs(q[order(abs(q))][ny])) * sign(q), 
                  0)
            }
            q <- q/drop(crossprod(tt))
            u <- as.matrix(Ytemp) %*% (q/drop(crossprod(q)))
            if (crossprod(tt - tt.old) < tol) 
                break
            if (iter == max.iter) {
                warning(paste("Maximum number of iterations reached for dimension", 
                  h), call. = FALSE)
                break
            }
            tt.old <- tt
            iter <- iter + 1
        }
        hiter <- c(hiter, iter)
        if (deflation == "block") {
            pb <- list()
            for (k in 1:nbloc) {
                pb[[k]] <- t(Xtemp[[k]]) %*% (tb[[k]]/drop(crossprod(tb[[k]])))
                Xtemp[[k]] <- as.matrix(Xtemp[[k]]) - as.matrix(tb[[k]]) %*% 
                  t(as.matrix(pb[[k]]))
                ExpVarX[k, h] <- 1 - Iner(Xtemp[[k]])/InerX[[k]]
            }
            Ytemp <- as.matrix(Ytemp) - as.matrix(tt) %*% t(q)
            ExpVarY[h] <- (1 - Iner(Ytemp))/InerY
        }
        else if (deflation == "super") {
            pb <- list()
            for (k in 1:nbloc) {
                pb[[k]] <- t(Xtemp[[k]]) %*% (tt/drop(crossprod(tt)))
                Xtemp[[k]] <- as.matrix(Xtemp[[k]]) - as.matrix(tt) %*% 
                  t(as.matrix(pb[[k]]))
                ExpVarX[k, h] <- 1 - Iner(Xtemp[[k]])/InerX[[k]]
            }
            Ytemp <- as.matrix(Ytemp) - as.matrix(tt) %*% t(q)
            ExpVarY[h] <- 1 - Iner(Ytemp)/InerY
        }
        else stop("non convenient deflation method!")
        Tt[[h]] <- tt
        Qt[[h]] <- q
        Pb[[h]] <- unlist(pb)
        Tb[[h]] <- unlist(tb)
        U[[h]] <- u
        Wb[[h]] <- unlist(wb)
        Wt[[h]] <- wt
    }
    Tt <- do.call("cbind", Tt)
	U <- do.call("cbind", U)
    Tb <- do.call("cbind", Tb)
	colnames(Tt) <- colnames(Tb) <- colnames(U) <- paste("CS",1:nf,sep="")
	rownames(U) <- rownames(Y)
	rownames(Tt) <- rownames(Y)
	rownames(Tb) <- util1$rownames
    Qt <- do.call("cbind", Qt)
    Wt <- do.call("cbind", Wt)
    Wb <- do.call("cbind", Wb)
	colnames(Qt) <- colnames(Wt) <- colnames(Wb) <- paste("W",1:nf,sep="")
	rownames(Qt) <- colnames(Y)
	rownames(Wb) <- util1$colnames
	rownames(Wt) <- util1$tabnames
    Pb <- do.call("cbind", Pb)
	colnames(Pb) <- paste("LD",1:nf,sep="")
	rownames(Pb) <- util1$colnames	
	res <- list(Tt, Qt, Wt, Pb, U, Wb, Tb, hiter)
    names(res) <- c("Tt", "Qt", "Wt", "Pb", "U", "Wb", "Tb","iter")
    res$Y <- Y
    res$X <- X
    res$InerY <- InerY
	names(InerX) <- util1$tabnames	
    res$InerX <- InerX
	rownames(ExpVarX) <- util1$tabnames
	colnames(ExpVarX) <- paste("CS",1:nf,sep="")	
    res$ExpVarX <- ExpVarX
	names(ExpVarY) <- paste("CS",1:nf,sep="")
    res$ExpVarY <- ExpVarY
    res$blo <- bloc
	res$TL <- util1$TL
	res$TC <- util1$TC
	rownames(keepX) <- util1$tabnames
	colnames(keepX) <- paste("CS",1:nf,sep="")		
    res$keepX <- keepX
	names(keepY) <- paste("CS",1:nf,sep="")
    res$keepY <- keepY
    res$call <- match.call()
	res$class <- "smbpls"
    return(res)
}
