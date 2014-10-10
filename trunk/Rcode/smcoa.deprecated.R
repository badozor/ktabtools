# ==============================================================================
# Title: smcoa-0.1.R
# Description: sparse Multiple co-inertia analysis (smcoa)
# Author:  Pierre Bady <pierre.bady@unil.ch>
# Date : 2012-11-07
# Version: 0.1
# Revision: 2012-11-07
# Comments: RAS
# License: GPL version 2 or newer
# Copyright (C) 2003-2013  Pierre Bady
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# ==============================================================================
### based on svd method used in sPCA from R package mixOmics
ssvd1 <- function(x,keepX = ncol(x),iter.max = 500, tol = 1e-06,...){
		n <- nrow(x)
		p <- ncol(x)
		X.temp <- x
       svd.X = svd(X.temp)
        u.new = svd.X$u[, 1]
        v.new = svd.X$d[1] * svd.X$v[, 1]
        v.stab = FALSE
        u.stab = FALSE
        iter = 0
        nx = p - keepX
        while ((v.stab == FALSE) || (u.stab == FALSE)) {
            iter = iter + 1
            u.old = u.new
            v.temp = t(X.temp) %*% u.old
            v.old = v.new
            if (nx != 0) {
                v.new = ifelse(abs(v.temp) > abs(v.temp[order(abs(v.temp))][nx]), 
                  (abs(v.temp) - abs(v.temp[order(abs(v.temp))][nx])) * 
                    sign(v.temp), 0)
            }
            u.new = as.vector(X.temp %*% v.new)
            u.new = u.new/sqrt(drop(crossprod(u.new)))
           if (crossprod(u.new - u.old) < tol) {
                u.stab = TRUE
            }
            if (crossprod(v.new - v.old) < tol) {
                v.stab = TRUE
            }
            if ((is.na(v.stab)) | (is.na(u.stab)) | (iter >= 
                iter.max)) {
                v.stab = TRUE
                u.stab = TRUE
            }
        }
	return(list(d=svd.X$d,u=svd.X$u,v=svd.X$v,u.new=u.new,v.new=v.new,keepX=keepX))
}

ssvdk <- function(x,bloc=ncol(x),keepX = bloc,iter.max = 500, tol = 1e-06,...){
nblo <- length(bloc)
indicablo <- rep.int(1:nblo,bloc)
n <- nrow(x)
p <- ncol(x)
X.temp <- x
       svd.X = svd(X.temp)
        u.new = svd.X$u[, 1]
        v.new = svd.X$d[1] * svd.X$v[, 1]
        v.stab = FALSE
        u.stab = FALSE
        iter = 0
        while ((v.stab == FALSE) || (u.stab == FALSE)) {
            iter = iter + 1
            u.old = u.new
            v.temp = t(X.temp) %*% u.old
            v.old = v.new
### adaptation for ktab
	v.new <- v.temp
	for(i in 1:nblo){
	nx <- bloc[i] - keepX[i]
            if (nx != 0) {
		auxi <- v.temp[indicablo==i] 
                v.new[indicablo==i] = ifelse(abs(auxi) > abs(auxi[order(abs(auxi))][nx]), 
                  (abs(auxi) - abs(auxi[order(abs(auxi))][nx])) * sign(auxi), 0)


            }
            
}
#####
            u.new = as.vector(X.temp %*% v.new)
            u.new = u.new/sqrt(drop(crossprod(u.new)))
   

	if (crossprod(u.new - u.old) < tol) {
                u.stab = TRUE
            }
            if (crossprod(v.new - v.old) < tol) {
                v.stab = TRUE
            }
            if ((is.na(v.stab)) | (is.na(u.stab)) | (iter >= iter.max)) {
                v.stab = TRUE
                u.stab = TRUE
            }
        }
	return(list(d=svd.X$d,u=svd.X$u,v=svd.X$v,u.new=u.new,v.new=v.new,keepX=keepX,indic=indicablo,bloc=bloc))
}
#sk <- ssvdk(X,bloc=c(2,4),keepX=c(1,1))
### pas top :/ #####
# probleme pour conserver l'orthogonalité => à vérifier
# => à revoir!!!!
#
smcoa <- function (X, keepX = X$blo,option = c("inertia", "lambda1", "uniform", "internal"),scannf = TRUE, nf = 3, tol = 1e-07,iter.max = 500) 
{
    if (!inherits(X, "ktab")) 
        stop("object 'ktab' expected")
    option <- option[1]
    if (option == "internal") {
        if (is.null(X$tabw)) {
            warning("Internal weights not found: uniform weigths are used")
            option <- "uniform"
        }
    }
    lw <- X$lw
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
    for (i in 1:nbloc) X[[i]] <- X[[i]] * sqrt(tabw[i])
    Xsepan <- sepan(X, nf = 4)
    normaliserparbloc <- function(scorcol) {
        for (i in 1:nbloc) {
            w1 <- scorcol[indicablo == i]
            w2 <- sqrt(sum(w1 * w1))
            if (w2 > tol) 
                w1 <- w1/w2
            scorcol[indicablo == i] <- w1
        }
        return(scorcol)
    }
    recalculer <- function(tab, scorcol) {
        for (k in 1:nbloc) {
            soustabk <- tab[, indicablo == k]
            uk <- scorcol[indicablo == k]
            soustabk.hat <- t(apply(soustabk, 1, function(x) sum(x * 
                uk) * uk))
            soustabk <- soustabk - soustabk.hat
            tab[, indicablo == k] <- soustabk
        }
        return(tab)
    }
    tab <- as.matrix(X[[1]])
    for (i in 2:nbloc) {
        tab <- cbind(tab, X[[i]])
    }
    names(tab) <- auxinames$col
    tab <- tab * sqrt(lw)
    tab <- t(t(tab) * sqrt(cw))
    compogene <- list()
    uknorme <- list()
    valsing <- NULL
    nfprovi <- min(c(20, nlig, ncol))
    for (i in 1:nfprovi) {
#### svd part with variable selection ####
        af <- svd(tab)
        u.new = af$u[, 1]
        v.new = af$d[1] * af$v[, 1]
        v.stab = FALSE
        u.stab = FALSE
        iter = 0
        while ((v.stab == FALSE) || (u.stab == FALSE)) {
            iter = iter + 1
            u.old = u.new
            v.temp = t(tab) %*% u.old
            v.old = v.new
### adaptation for ktab
	v.new <- v.temp
	for(j in 1:nbloc){
	nx <- bloc[j] - keepX[j]
            if (nx != 0) {
		auxi <- v.temp[indicablo==j] 
                v.new[indicablo==j] = ifelse(abs(auxi) > abs(auxi[order(abs(auxi))][nx]), 
                  (abs(auxi) - abs(auxi[order(abs(auxi))][nx])) * sign(auxi), 0)
            }
	}
        u.new = as.vector(tab %*% v.new)
        u.new = u.new/sqrt(drop(crossprod(u.new)))
	if (crossprod(u.new - u.old) < tol) {
                u.stab = TRUE
            }
            if (crossprod(v.new - v.old) < tol) {
                v.stab = TRUE
            }
            if ((is.na(v.stab)) | (is.na(u.stab)) | (iter >= iter.max)) {
                v.stab = TRUE
                u.stab = TRUE
            }
        }
	compogene[[i]] <- u.new
	w <- normaliserparbloc(v.new)
### recalcule du tableau ###
        tab <- recalculer(tab, w)
        w <- w/sqrt(cw)
        uknorme[[i]] <- w
        w <- af$d[1]
        valsing <- c(valsing, w)
    }
    pseudoeig <- valsing^2
    if (scannf) {
        barplot(pseudoeig)
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0) 
        nf <- 2
    acom <- list()
    acom$pseudoeig <- pseudoeig
    w <- matrix(0, nbloc, nf)
    for (i in 1:nbloc) {
        w1 <- Xsepan$Eig[rank.fac == i]
        r0 <- Xsepan$rank[i]
        if (r0 > nf) 
            r0 <- nf
        w[i, 1:r0] <- w1[1:r0]
    }
    w <- data.frame(w)
    row.names(w) <- Xsepan$tab.names
    names(w) <- paste("lam", 1:nf, sep = "")
    acom$lambda <- w
    w <- matrix(0, nlig, nf)
    for (j in 1:nf) w[, j] <- compogene[[j]]
    w <- data.frame(w)
    names(w) <- paste("SynVar", 1:nf, sep = "")
    row.names(w) <- row.names(X)
    acom$SynVar <- w
    w <- matrix(0, ncol, nf)
    for (j in 1:nf) w[, j] <- uknorme[[j]]
    w <- data.frame(w)
    names(w) <- paste("Axis", 1:nf, sep = "")
    row.names(w) <- auxinames$col
    acom$axis <- w
    w <- matrix(0, nlig * nbloc, nf)
    covar <- matrix(0, nbloc, nf)
    i1 <- 0
    i2 <- 0
# recalcul de SynVar pour conserver l'orthogonalite
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + nlig
        urk <- as.matrix(acom$axis[indicablo == k, ])
        tab <- as.matrix(X[[k]])
        urk <- urk * cw[indicablo == k]
        urk <- tab %*% urk
        w[i1:i2, ] <- urk
        urk <- urk * acom$SynVar * lw
        covar[k, ] <- apply(urk, 2, sum)
    }
#########
    w <- data.frame(w, row.names = auxinames$row)
    names(w) <- paste("Axis", 1:nf, sep = "")
    acom$Tli <- w
    covar <- data.frame(covar)
    row.names(covar) <- tab.names(X)
    names(covar) <- paste("cov2", 1:nf, sep = "")
    acom$cov2 <- covar^2
    w <- matrix(0, nlig * nbloc, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + nlig
        tab <- acom$Tli[i1:i2, ]
        tab <- scalewt(tab, wt = lw, center = FALSE, scale = TRUE)
        w[i1:i2, ] <- tab
    }
    w <- data.frame(w, row.names = auxinames$row)
    names(w) <- paste("Axis", 1:nf, sep = "")
    acom$Tl1 <- w
    w <- matrix(0, ncol, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + ncol(X[[k]])
        urk <- as.matrix(acom$SynVar)
        tab <- as.matrix(X[[k]])
        urk <- urk * lw
        w[i1:i2, ] <- t(tab) %*% urk
    }
    w <- data.frame(w, row.names = auxinames$col)
    names(w) <- paste("SV", 1:nf, sep = "")
    acom$Tco <- w
    acom$TSco <- w
    for(i in 1:nf){
        auxi <- acom$axis[,i]==0
        acom$TSco[auxi,i] <- rep(0,sum(auxi))
        }
    var.names <- NULL
    w <- matrix(0, nbloc * 4, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + 4
        urk <- as.matrix(acom$axis[indicablo == k, ])
        tab <- as.matrix(Xsepan$C1[indicablo == k, ])
        urk <- urk * cw[indicablo == k]
        tab <- t(tab) %*% urk
        for (i in 1:min(nf, 4)) {
            if (tab[i, i] < 0) {
                for (j in 1:nf) tab[i, j] <- -tab[i, j]
            }
        }
        w[i1:i2, ] <- tab
        var.names <- c(var.names, paste(Xsepan$tab.names[k], 
            ".a", 1:4, sep = ""))
    }
    w <- data.frame(w, row.names = auxinames$tab)
    names(w) <- paste("Axis", 1:nf, sep = "")
    acom$Tax <- w
    acom$nf <- nf
    acom$TL <- X$TL
    acom$TC <- X$TC
    acom$T4 <- X$T4
    acom$keepX <- keepX
    class(acom) <- "mcoa"
    acom$call <- match.call()
    return(acom)
}
plot.smcoa <- function(x,...){
}
kplot.smcoa <- function(x,...){
}
summary.smcoa <- function(x,...){
}
# ============================== END! ==========================================
