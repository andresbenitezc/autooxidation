####
# Author: Andres Benitez
####
#R CMD INSTALL /Library/Frameworks/R.framework/Versions/3.2/Resources/library deconspectra
# detach("package:deconspectra", unload = TRUE)

deconSpectra <- function(S, W, R, M = NULL, baseline = TRUE, lambda = TRUE, start = NULL, ...){

	if (is.null(M)){ M = as.numeric(gsub('[A-Z]+','',names(S))) }
	if (length(S) != length(M) || !is.numeric(M)) { stop('Time series length does not match time values vector length') } 
	if (length(S[,1]) != length(W)) { stop('spectrum lengths do not match wavelength vector length') }	
	if (length(S[,1]) != length(R[,1])) { stop('spectrum lengths do not match reference spectrum length') }	
	
	out <- structure(list(fits = NULL, values = NULL, species = c(names(R)), wavelength = W, time = M), class = 'nlsArray')
	if (!is.null(start)){
		out$values <- start
	}else{
		out$values <- data.frame(as.list(c(1, rep(0, length(R) - 1))))
		names(out$values) <- c(names(R))
	}
	if (baseline == T && is.null(out$values[['base']]) ) { out$values[['base']] = 0 }
	if (lambda == T && is.null(out$values[['lambda']]) ) { out$values[['lambda']] = 0 }

	k = 0	
	for (i in 1:length(S)) {
	
		if (i == 1){k = 1}
		else {k = i-1}
	
		trace = data.frame(R, absorb = S[[i]], wavelength = W)
		start.fom = paste(letters[1:length(names(R))], names(R), sep = '*')
		
		decon.start = NULL
		decon.start[c(letters[1:length(names(R))])] = out$values[k,1:(length(names(R)))]
		decon.start = as.list(decon.start)
		decon.form = paste(letters[1:length(names(R))], names(R), sep = '*')
		decon.form = paste(decon.form, collapse = ' + ')
		decon.form = paste('absorb', decon.form, sep = ' ~ ')
		decon.lower = rep(-1, length(R))
		if (baseline == T){
			decon.start['z'] = out$values[['base']][k]
			decon.form = paste(decon.form, 'z', sep = ' + ')
			decon.lower = c(decon.lower, - 10)
		}
		if (lambda == T){ 
			decon.start['y'] = out$values[['lambda']][k]
			decon.form = paste(decon.form, 'y*10^11/wavelength^4', sep = ' + ')
			decon.lower = c(decon.lower, - 10)
		}		
		decon.solution <- nlsLM(as.formula(decon.form), data = trace, start = decon.start, lower = decon.lower,
			control = list(maxiter = 1000))
		out$fits[[i]] = decon.solution
		out$values[i,] <- coef(decon.solution)
	}
	out
}
