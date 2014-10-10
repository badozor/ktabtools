statismod <-
function (X, scannf = TRUE, nf = 3, tol = 1e-07) 
{
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
        wk <- as.matrix(X[[k]]) * lwsqrt
        wk <- t(t(wk) * ak)
        wk <- wk %*% t(wk)
        sep[[k]] <- wk
    }
    sep <- matrix(unlist(sep), nlig * nlig, ntab)
    RV <- t(sep) %*% sep
    ak <- sqrt(diag(RV))
    RV <- sweep(RV, 1, ak, "/")
    RV <- sweep(RV, 2, ak, "/")
    dimnames(RV) <- list(tab.names, tab.names)
    statis$RV <- RV
    eig1 <- eigen(RV, symmetric = TRUE)
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
    sep1 <- lapply(split(sep, col(sep)), function(x) matrix(x, 
        nlig, nlig))
    C.ro <- apply(t(sep) * tabw, 2, sum)
    C.ro <- matrix(unlist(C.ro), nlig, nlig)
    eig1 <- eigen(C.ro, symmetric = TRUE)
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
    vec2 <- vec1/lwsqrt
    vec2 <- t(t(vec2)/sqrt(eig[1:nf]))
    L1 <- NULL
    for (k in 1:ntab) {
        w <- sep1[[k]]
        w <- (as.matrix(w)) %*% as.matrix(vec2)
        L1 <- rbind(L1, w)
    }
    L1mod <- NULL
    for (k in 1:ntab) {
        w <- sep1[[k]]
        w <- (as.matrix(w)) %*% as.matrix(vec2)
        tab <- as.matrix(sepanL1[X$TL[, 1] == k, ])
        tab <- t(tab * lw) %*% vec2
        ak <- apply(tab, 2, function(x) sqrt(sum(x * x)))
        w <- sweep(w, 2, ak, "/")
        L1mod <- rbind(L1mod, w)
    }
    w <- C.ro
    w <- (as.matrix(w)) %*% as.matrix(vec2)
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
    class(statis) <- "statis"
    return(statis)
}
