# ==============================================================================
# Title: bootRV.R
# Description:
# Author:  Pierre Bady <pierre.bady@unil.ch>
# Date : 
# Version: 
# Revision: Jan 8, 2013 10:00:24 AM
# Comments: RAS
# License: GPL version 2 or newer
# Copyright (C) 2013  Pierre Bady
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

### slight modification of the original function RV.rtest from R package ade4
# by default center=TRUE as in the function RV.rtest
bootRV <- function (df1, df2, nrepet = 99,center=TRUE){
	if (!is.data.frame(df1)) 
		stop("data.frame expected")
	if (!is.data.frame(df2)) 
		stop("data.frame expected")
	l1 <- nrow(df1)
	if (nrow(df2) != l1) 
		stop("Row numbers are different")
	if (any(row.names(df2) != row.names(df1))) 
		stop("row names are different")		
	RV <- function(X,Y,center=TRUE){
		if(center){    
			X <- scale(X, scale = FALSE)
			Y <- scale(Y, scale = FALSE)
		}else{
			X <- X
			Y <- Y
		}
		X <- X/(sum(svd(X)$d^4)^0.25)
		Y <- Y/(sum(svd(Y)$d^4)^0.25)
		X <- as.matrix(X)
		Y <- as.matrix(Y)
		sum(svd(t(X) %*% Y)$d^2)			
	}
	obs <- RV(df1,df2,center=center)
	if (nrepet == 0) 
		return(obs)
	perm <- matrix(0, nrow = nrepet, ncol = 1)
	perm <- apply(perm, 1, function(x) {vec1 <-sample(l1,l1,replace=TRUE); RV(df1[vec1,],df2[vec1,],center=center)})
	w <- list(obs = obs, sim = perm,mean= mean(perm),se=sd(perm),bias=mean(perm)-obs ,call = match.call())
	return(w)
}
	
		
# ========================= END! ===============================================
