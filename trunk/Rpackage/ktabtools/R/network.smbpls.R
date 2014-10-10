network.smbpls <-
function (object, comp = 1, X.names = NULL, Y.names = NULL, keep.var = TRUE, 
    threshold = 0.1, ...) 
{
    nbloc <- length(object$blo)
    bloc <- object$blo
    indicablo <- object$TC[,1]
	keepX <- object$keepX
    simMat <- list()
    if (is.null(keepX) & is.null(keepX)) {
        for (k in 1:nbloc){
	w <- object$Wb[object$TC[,1] == k, 1:comp] %*% t(object$Qt[, 1:comp])
	colnames(w) <- rownames(object$Qt)
	rownames(w) <- rownames(object$Wb[object$TC[,1] == k,])
	simMat[[k]] <- w	
	}
    }
    else if (is.null(keepX) & !is.null(keepX)) {
        .NotYetImplemented()
    }
    else if (!is.null(keepX) & is.null(keepX)) {
        .NotYetImplemented()
    }
    simMat <- as.matrix(do.call("rbind", simMat))
    mixOmics:::network.default(simMat, threshold = threshold,...)
}
