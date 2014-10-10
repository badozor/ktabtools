# ==============================================================================
# Title: mbplsutil.R
# Description:
# Author:  Pierre Bady <pierre.bady@unil.ch>
# Date : 
# Version: 
# Revision: Nov 16, 2012 4:22:41 PM
# Comments: RAS
# License: GPL version 2 or newer
# Copyright (C) 2012  Pierre Bady
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

print.mbpls <- print.smbpls <- function(...){
	.NotYetImplemented()
}
summary.mbpls <- summary.smbpls <- function(...){
	.NotYetImplemented()
}
predict.mbpls <- predict.smbpls <- function(...){
	.NotYetImplemented()
}
plot.mbpls <- plot.smbpls <- function(...){
	.NotYetImplemented()
}
kplot.mbpls <- kplot.smbpls <- function(...){
	.NotYetImplemented()
}
valid.mbpls <- valid.smbpls <- function(...){
	.NotYetImplemented()
}
network.mbpls <- network.smbpls <- function(object, comp = 1, X.names = NULL, Y.names = NULL, keep.var = TRUE,threshold = 0.1, ...) 		
{
	nbloc <- length(object$blo)
	bloc <- object$blo
	indicablo <- object$TC[, 1]
	simMat <- list()
	if(is.null(keepX) & is.null(keepX)){
		for(k in 1:nbloc)
			simMat[[k]] <- object$Wb[object$TC[,1]==k,1:comp] %*% t(object$Qt[,1:comp])
	}else if(is.null(keepX) & !is.null(keepX)){
		.NotYetImplemented()
	}else if(!is.null(keepX) & is.null(keepX)){
		.NotYetImplemented()
	}
	simMat <- as.matrix(do.call("rbind",simMat))
	mixOmics:::network.default(simMat, threshold=threshold,...)
}

# ========================== END! ==================================================
