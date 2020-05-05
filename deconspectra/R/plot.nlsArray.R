plot.nlsArray <- function(x, col=NULL, normalized = FALSE, ...) {
	if (is.null(col)){col = rainbow(length(names(x$values)))}

	if (normalized == T){ x.values = data.frame(mapply(function(x){ if(max(abs(x)) != 0){ x/max(abs(x)) } else{x} }, x$values))}
	else { x.values = x$values }

	maxY = max(x.values)
	minY = min(c(0,min(x.values)))

	plot(x$time, ylim = c(minY, maxY), unlist(x.values[1]), col = col[1], ...)
	for (i in 2:length(names(x.values))){
		points(x$time, unlist(x.values[i]), col = col[i])
	}
}

