# ==============================================================================
# Title: RVrtest.R
# Description:
# Author:  Pierre Bady <pierre.bady@unil.ch>
# Date : 
# Version: 
# Revision: Jan 8, 2013 10:03:19 AM
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
RVrtest <- function (df1, df2, nrepet = 99,center=TRUE) 
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


# ========================= END! ===============================================
