klistutil <- function (x){
    bloc <- unlist(lapply(x, ncol))
    nbloc <- length(bloc)
	TC <- matrix(0,ncol=2,nrow=sum(bloc))
    TC[,1] <- rep.int(1:nbloc,bloc)
	w <- NULL
	for(i in 1:nbloc)
		w <- c(w,1:bloc[i])
	TC[,2] <- w
	nligne <- unlist(lapply(x, nrow))
	TL <- matrix(0,ncol=2,nrow=sum(nligne))
    TL[,1] <- rep.int(1:nbloc, nligne)
	w <- NULL
	for(i in 1:nbloc)
		w <- c(w,1:nligne[i])
	TL[,2] <- w
	w <- unlist(lapply(x,colnames))
    if (any(duplicated(w))) 
        w <- paste(w, as.character(TC[, 1]), sep = ".")
	rownames(TC) <- w
	w <- unlist(lapply(x,rownames))
	if (any(duplicated(w))) 
        w <- paste(w, as.character(TL[, 1]), sep = ".")	
	rownames(TL) <- w
	colnames(TL) <- colnames(TC) <- c("I1","I2")
	if(is.null(names(x))){
		tabnames <- paste("X",1:nbloc,sep=".")
	}else{
		tabnames <- names(x)
	}
    return(list(bloc=bloc,TL = TL, TC = TC, tabnames = tabnames,rownames= rownames(TL),colnames=rownames(TC)))
}