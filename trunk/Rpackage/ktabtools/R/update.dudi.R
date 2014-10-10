update.dudi <- function(dudi,...,evaluate=TRUE){
	if (!inherits(dudi, "dudi")) 
		stop("Object of class 'dudi' expected")
	call <- dudi$call
	if (is.null(call)) 
		stop("need an object with call component")
	extras <- match.call(expand.dots = FALSE)$...
	if (length(extras) > 0) {
		existing <- !is.na(match(names(extras), names(call)))
		for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
		if (any(!existing)) {
			call <- c(as.list(call), extras[!existing])
			call <- as.call(call)
		}
	}
	if (evaluate) 
		eval(call, parent.frame())
	else call
}	