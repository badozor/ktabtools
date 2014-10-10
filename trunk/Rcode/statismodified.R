# ==============================================================================
# Title: statismod-0.1.R
# Description: statis with trajectories and additional utility functions
# Author:  Pierre Bady <pierre.bady@unil.ch>
# Date : 30/09/2003 17:24:42
# Version: 0.1
# Revision: 07/09/2010
# Comments: RAS
# License: GPL version 2 or newer
# Copyright (C) 2003-2010  Pierre Bady
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

# from statis of the R package ade4 (Chesselet al. 2004)
# addition of lavit and modified lavit trajectories
# (pedagogically interesting, but MCOA and MFA are 
# more efficent for the representation of tables onto the compromise)

 statismod <- function (X, scannf = TRUE, nf = 3, tol = 1e-07){
    if (!inherits(X, "ktab"))
        stop("object 'ktab' expected")
    lw <- X$lw
    nlig <- length(lw)
    cw <- X$cw
    ncol <- length(cw)
    blo <- X$blo
    ntab <- length(X$blo)
    indicablo <- X$TC[, 1]
    tab.names <- tab.names(X)
    auxinames <- ktab.util.names(X)
    statis <- list()
    sep <- list()
    lwsqrt <- sqrt(lw)
    for (k in 1:ntab) {
        ak <- sqrt(cw[indicablo == k])
# integration du poids des lignes par lwsqrt
        wk <- as.matrix(X[[k]]) * lwsqrt
        wk <- t(t(wk) * ak)
        wk <- wk %*% t(wk)
        sep[[k]] <- wk
    }
# sep contient les WkD=XkQkXktD
    sep <- matrix(unlist(sep), nlig * nlig, ntab)
    RV <- t(sep) %*% sep
    ak <- sqrt(diag(RV))
    RV <- sweep(RV, 1, ak, "/")
    RV <- sweep(RV, 2, ak, "/")
    dimnames(RV) <- list(tab.names, tab.names)
    statis$RV <- RV
    eig1 <- eigen(RV, sym = TRUE)
    statis$RV.eig <- eig1$values
    if (any(eig1$vectors[, 1] < 0))
        eig1$vectors[, 1] <- -eig1$vectors[, 1]
    tabw <- eig1$vectors[, 1]
    statis$RV.tabw <- tabw
    w <- t(t(eig1$vectors) * sqrt(eig1$values))
    w <- as.data.frame(w)
    row.names(w) <- tab.names
    names(w) <- paste("S", 1:ncol(w), sep = "")
    statis$RV.coo <- w[, 1:min(4, ncol(w))]
    sep <- t(t(sep)/ak)
# remplacer sep par sep1 apres la diagonalisation du compromis	
    sep1 <- lapply(split(sep,col(sep)),function(x) matrix(x,nlig,nlig))
    C.ro <- apply(t(sep) * tabw, 2, sum)
    C.ro <- matrix(unlist(C.ro), nlig, nlig)
    eig1 <- eigen(C.ro, sym = TRUE)
    eig <- eig1$values
    rank <- sum((eig/eig[1]) > tol)
    if (scannf) {
        barplot(eig[1:rank])
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0)
        nf <- 2
    if (nf > rank)
        nf <- rank
    statis$C.eig <- eig[1:rank]
    statis$C.nf <- nf
    statis$C.rank <- rank
    wref <- eig1$vectors[, 1:nf]
# vec2 == wref
    vec1 <- wref
    wref <- wref/lwsqrt
    w <- data.frame(t(t(wref) * sqrt(eig[1:nf])))
    row.names(w) <- row.names(X)
    names(w) <- paste("C", 1:nf, sep = "")
    statis$C.li <- w
    w <- as.matrix(X[[1]])
    for (k in 2:ntab) {
        w <- cbind(w, as.matrix(X[[k]]))
    }
    w <- w * lw
    w <- t(w) %*% wref
    w <- data.frame(w, row.names = auxinames$col)
    names(w) <- paste("C", 1:nf, sep = "")
    statis$C.Co <- w
    sepanL1 <- sepan(X, nf = 4)$L1
    w <- matrix(0, ntab * 4, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:ntab) {
        i1 <- i2 + 1
        i2 <- i2 + 4
        tab <- as.matrix(sepanL1[X$TL[, 1] == k, ])
        tab <- t(tab * lw) %*% wref
        for (i in 1:min(nf, 4)) {
            if (tab[i, i] < 0) {
                for (j in 1:nf) tab[i, j] <- -tab[i, j]
            }
        }
        w[i1:i2, ] <- tab
    }
    w <- data.frame(w, row.names = auxinames$tab)
    names(w) <- paste("C", 1:nf, sep = "")
    statis$C.T4 <- w
    w <- as.matrix(statis$C.li) * lwsqrt
    w <- w %*% t(w)
    w <- w/sqrt(sum(w * w))
    w <- as.vector(unlist(w))
    sep <- sep * unlist(w)
    w <- apply(sep, 2, sum)
    statis$cos2 <- w
    statis$tab.names <- tab.names
    statis$TL <- X$TL
    statis$TC <- X$TC
    statis$T4 <- X$T4
# traj. de lavit
	vec2 <- vec1/lwsqrt
	vec2 <- t(t(vec2) / sqrt(eig[1:nf]))
	L1 <- NULL
	for(k in 1:ntab){
		w <- sep1[[k]]
		w <- (as.matrix(w))%*%as.matrix(vec2)
		L1 <- rbind(L1,w)
	}
# traj. de lavit modifiees, see Thema51 (Chessel et al.)
	L1mod <- NULL
	for (k in 1:ntab) {
		w <- sep1[[k]]
		w <- (as.matrix(w))%*%as.matrix(vec2)       
		tab <- as.matrix(sepanL1[X$TL[, 1] == k, ])
		tab <- t(tab * lw) %*% vec2
		ak <- apply(tab, 2, function(x) sqrt(sum(x*x)))
		w <- sweep(w,2,ak,"/")
		L1mod <- rbind(L1mod,w)	
	}
# ajouter les trajectoires de Place

# ajouter les trajectoires basees sur l'analyse procuste

#----------------
w <- C.ro
w <- (as.matrix(w))%*%as.matrix(vec2)
w <- as.data.frame(w)
rownames(w) <- rownames(statis$C.li)
colnames(w) <- colnames(statis$C.li)
statis$C.la <- w
L1 <- as.data.frame(L1)
colnames(L1) <- colnames(statis$C.li)
rownames(L1) <- auxinames$row
L1mod <- as.data.frame(L1mod)
colnames(L1mod) <- colnames(statis$C.li)
rownames(L1mod) <- auxinames$row
statis$La <- L1
statis$Lamod <- L1mod
#statis$sep1 <- sep1
#statis$C.ro <- C.ro
#==============================
    class(statis) <- "statis"
    return(statis)
}
# some tests:
# statis3 <- statismod(ktab1,nf=2,scannf=FALSE)
# s.match(statis3$La[statis3$TL[,1]==1,],statis3$C.la)
# summary(statis3$La)
# summary(statis3$C.la)
# summary(statis3$C.li)
# summary(statis3$La[statis3$TL[,1]==1,])
# summary(statis3$La[statis3$TL[,1]==2,])
# summary(statis3$La[statis3$TL[,1]==3,])
# s.class(statis3$Lamod,fac=statis3$TL[,2])
# points(statis3$C.la,pch=21,bg="red")
# ww <- statis3$La[statis3$TL[,1]==1,]*statis3$RV.tabw[1]+
# statis3$La[statis3$TL[,1]==2,]*statis3$RV.tabw[2]+
# statis3$La[statis3$TL[,1]==3,]*statis3$RV.tabw[3]
# summary(ww)
# hcc <-hclust(dist(statis3$C.li),method="ward")
# ddc <- as.dendrogram(hcc)
# ddc <- reorder(ddc, statis3$C.li[,1])
# colInd <- order.dendrogram(ddc)
# plot(ddc)

######################
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
kRVplot <- function(X,whichinrow = NULL, whichincol = NULL,gap=4,nclass = 10,coeff = 1,nrepet=99,col=jet.colors(10),digits=2,...){
    "hist.simul.util" <- function(sim, obs, nclass, coeff, title = "") {
        r0 <- c(sim, obs)
        h0 <- hist(sim, plot = FALSE, nclass = nclass, xlim = xlim0)
        y0 <- max(h0$counts)
        l0 <- max(sim) - min(sim)
        w0 <- l0/(log(length(sim), base = 2) + 1)
        w0 <- w0 * coeff
        xlim0 <- range(r0) + c(-w0, w0)
        hist(sim, plot = TRUE, nclass = nclass, xlim = xlim0, 
            main = title, col = grey(0.9))
        lines(c(obs, obs), c(y0/2, 0))
        points(obs, y0/2, pch = 18, cex = 2)
    }
	if (!inherits(X, "ktab")) 
        stop("object 'ktab' expected")
    lw <- X$lw
    nlig <- length(lw)
    cw <- X$cw
    ncol <- length(cw)
    ntab <- length(X$blo)
    indicablo <- X$TC[, 1]
    labels <- tab.names(X)
    sep <- list()
    lwsqrt <- sqrt(lw)
    for (k in 1:ntab) {
        ak <- sqrt(cw[indicablo == k])
        wk <- as.matrix(X[[k]]) * lwsqrt
        wk <- t(t(wk) * ak)
        wk <- wk %*% t(wk)
        sep[[k]] <- wk
    }
	if(is.null(whichinrow)) 
		whichinrow <- 1:ntab
	if(is.null(whichincol))
		whichincol <- 1:ntab
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    oma <- c(2, 2, 1, 1)
    par(mfrow = c(length(whichinrow), length(whichincol)), mar = rep(gap/2,4), oma = oma)
   for (i in whichinrow) {
        for (j in whichincol) {
            if (i == j) {
                plot.default(0, 0, type = "n", asp = 1, xlab = "", 
                  ylab = "", xaxt = "n", yaxt = "n", xlim = c(0, 
                    1), ylim = c(0, 1), xaxs = "i", yaxs = "i", 
                  frame.plot = FALSE)
                l.wid <- strwidth(labels, "user")
                cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                text(0.5, 0.5, labels[i], cex = cex.labels, font = 1)
            }
            else if (i > j) {
				rv1 <- RV.rtest(as.data.frame(sep[[i]]),as.data.frame(sep[[j]]),nrepet=nrepet)
				sim <- rv1$sim 
				obs <- rv1$obs
				titre <- paste("RV = ",round(obs,digits)," (p = ",rv1$pvalue,")",sep="")
                hist.simul.util(sim, obs, title = titre, nclass = nclass, 
                  coeff = coeff)
            }
            else if (j > i) {
				image(sep[[i]]%*%t(sep)[[j]],col=col,...)
            }
        }
    }
}


################################ END ! #####################"
