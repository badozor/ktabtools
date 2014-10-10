kRVplot.R <-
function (X, whichinrow = NULL, whichincol = NULL, gap = 4, nclass = 10, 
    coeff = 1, nrepet = 99, col = jet.colors(10), digits = 2, 
    center=TRUE,...) 
{
    "hist.simul.util" <- function(sim, obs, nclass, coeff, title = "") {
        r0 <- c(sim, obs)
        h0 <- hist(sim, plot = FALSE, nclass = nclass)
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
### slight modification of the original function RV.rtest from R package
RVrtest <- function (df1, df2, nrepet = 99,center=FALSE) 
{
    if (!is.data.frame(df1)) 
        stop("data.frame expected")
    if (!is.data.frame(df2)) 
        stop("data.frame expected")
    l1 <- nrow(df1)
    if (nrow(df2) != l1) 
        stop("Row numbers are different")
    if (any(row.names(df2) != row.names(df1))) 
        stop("row names are different")
if(center){    
X <- scale(df1, scale = FALSE)
    Y <- scale(df2, scale = FALSE)
}else{
X <- df1
Y <- df2
}
    X <- X/(sum(svd(X)$d^4)^0.25)
    Y <- Y/(sum(svd(Y)$d^4)^0.25)
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    obs <- sum(svd(t(X) %*% Y)$d^2)
    if (nrepet == 0) 
        return(obs)
    perm <- matrix(0, nrow = nrepet, ncol = 1)
    perm <- apply(perm, 1, function(x) sum(svd(t(X) %*% Y[sample(l1), 
        ])$d^2))
    w <- as.rtest(obs = obs, sim = perm, call = match.call())
    return(w)
}
###
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
    if (is.null(whichinrow)) 
        whichinrow <- 1:ntab
    if (is.null(whichincol)) 
        whichincol <- 1:ntab
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    oma <- c(2, 2, 1, 1)
    par(mfrow = c(length(whichinrow), length(whichincol)), mar = rep(gap/2, 
        4), oma = oma)
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
                rv1 <- RVrtest(as.data.frame(sep[[i]]), as.data.frame(sep[[j]]), 
                  nrepet = nrepet,center=center)
                sim <- rv1$sim
                obs <- rv1$obs
                titre <- paste("RV = ", round(obs, digits), " (p = ", 
                  rv1$pvalue, ")", sep = "")
                hist.simul.util(sim, obs, title = titre, nclass = nclass, 
                  coeff = coeff)
            }
            else if (j > i) {
                image(sep[[i]] %*% t(sep)[[j]], col = col, ...)
            }
        }
    }
}


