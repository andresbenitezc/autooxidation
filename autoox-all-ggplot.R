rm(list=ls())
#Read data for growth rates and kinetics.
grates = read.table('GrowthRates_Norm3-2.txt',header = TRUE)
kinetics  = read.table('TetX2_kinetics3-2.txt',header = TRUE)
runs = c('B1','B2','B3','B4','B5','B6','B7','B8')

#Set up variables and constants.
bgr = grates[,10]
min = grates[,1]
#Kinetics constants.
kcat = kinetics[,3]
tetx = kinetics[,5]/kinetics[,3]
knadph = kinetics[,2]
kmin = kinetics[,1]
grmax = 1.0
#Growth rate data.
x1 = grates[,1]
y2 = grates[,2:9]
y1 = asin((grates[,2:9])^0.5)
y1[is.nan(y1)] = asin(1)
nadph = 150
alpha = 4

#Fit the no-TetX growth curve using a hill equation; function is the growth rate dependent on minocycline inside the cell.
MinGrowth = function(pargr,min){asin((1-min^pargr[2]/(pargr[1]+min^pargr[2]))^0.5)}
MinGrowth_asin2 = function(pargr,min){1-min^pargr[2]/(pargr[1]+min^pargr[2])}

#Set up equations for calculating minocycline concentration inside the cell.
TetX_Vmax = function(D,kcat,tetx,knadph,nadph){(kcat*tetx/D)*nadph/(knadph+nadph)}
TetX_Vmax_2 = function(D,kcat,knadph,nadph){(kcat/D)*nadph/(knadph+nadph)}
TetX_kM = function(kmin,knadph,nadph){kmin*nadph/(knadph+nadph)}
MIN_Free = function(a,b,c){(-b-(b^2-4*a*c)^0.5)/(2*a)}
MIN_Free2 = function(vmax,km,min0,alpha){(-(min0-km/alpha-vmax)-((min0-km/alpha-vmax)^2+4/alpha*km*min0)^0.5)/(-2/alpha)}

#Fit all data with dr, nadph, and tetx as variables. tetx can vary between mutants.
fmin <- function(p,kcat,tetx,nadph,knadph,kmin,min0){
	vmax <- TetX_Vmax(p,kcat,tetx,knadph,nadph)
	km <- TetX_kM(kmin,knadph,nadph)
	km2 <- km/alpha
	a <- -1/alpha
	b <- outer(km2+vmax,min0,function(x,y){y-x})
	c <- outer(km,min0)
	MIN_Free(a,b,c)
}
fsum <- function(p){
	sum((t(y1[,-c]) - MinGrowth(p[1:2],fmin(p[3],kcat[-c],tetx[-c],nadph,knadph[-c],kmin[-c],x1)))^2)
}
fmin2 <- function(d,kcat_f,nadph_f,knadph_f,kmin_f,min0){
	vmax <- TetX_Vmax_2(d,kcat_f,knadph_f,nadph_f)
	km <- TetX_kM(kmin_f,knadph_f,nadph_f)
	km2 <- km/alpha
	a <- -1/alpha
	b <- min0-km2-vmax
  	c <- km*min0
  	MIN_Free(a,b,c)
}
fsum2 <- function(p){
	kmin2 <- p[1]
  	kcat2 <- p[1]*p[2]
  	gr2 <- y1[,i]
  	sum((t(gr2) - MinGrowth(param[1:2],fmin2(param[3],kcat2,nadph,knadph[i],kmin2,x1)))^2)
}
fmin3 <- function(d,kcat_f,nadph_f,knadph_f,kmin_f,min0){
	vmax <- TetX_Vmax_2(d,kcat_f,knadph_f,nadph_f)
	km <- TetX_kM(kmin_f,knadph_f,nadph_f)
	km2 <- km/alpha
	a <- -1/alpha
	b <- min0-km2-vmax
  	c <- km*min0
  	MIN_Free(a,b,c)
}
fsum3 <- function(p){
  	kmin2 <- p
  	kcat2 <- kact_m[i]*p
  	gr2 <- y1[,i]
  	sum((t(gr2) - MinGrowth(param[1:2],fmin3(param[3],kcat2,nadph,knadph[i],kmin2,x1)))^2)
}
fmin4 <- function(p,kcat,nadph,knadph,kmin,min0){
	vmax <- TetX_Vmax_2(p,kcat,knadph,nadph)
	km <- TetX_kM(kmin,knadph,nadph)
	km2 <- km/alpha
	a <- -1/alpha
	b <- outer(km2+vmax,min0,function(x,y){y-x})
	c <- outer(km,min0)
	MIN_Free(a,b,c)
}
fsum4 <- function(p){
	sum((t(y1) - MinGrowth(p[1:2],fmin4(p[3],kcat,nadph,knadph,kmin,x1)))^2)
}	

parameters <- NULL
for(c in c(1:8)){	
	#Model fitting using non-linear least squares.
	param = NULL
	out = optim(fsum,p = c(20,2.3,0.001), hessian = 1);
	print(out$convergence)
	param = out$par
	print(param[1:2]) 
	print(param[3])
	p <- param
	parameters <- rbind(parameters,p)
	run <- runs[c]
	assign(paste('out', run, sep = '_'),out)


	kmin_p <- NULL
	kact_p <- NULL
	index <- NULL
	minsum <- NULL
	n = 50
	for (k in 1:n){
		for (j in 1:n){
      			kmin_pt <- NULL
      			kact_pt <- NULL
      			minsum_t <- NULL
      			for (i in 1:length(kcat)){
        			out2 = optim(fsum2,p = c((k)*2+10,(j)*0.001+0.001),hessian = 1);
        			print(c(out2$convergence,k,j))
        			kmin_pt <- rbind(kmin_pt,out2$par[1])
        			kact_pt <- rbind(kact_pt,out2$par[2])
        			minsum_t <- rbind(minsum_t,out2$value)
      			}
      			index <- rbind(index,c((k)*2+10,(j)*0.001))
      			kmin_p <- cbind(kmin_p,kmin_pt)
      			kact_p <- cbind(kact_p,kact_pt)
      			minsum <- cbind(minsum,minsum_t)
    		}
	}

	varB <- NULL
	for (i in 1:(length(x1)-4)){
		varb <- (colSums(x1[i:(i+4)]*y1[i:(i+4),]) - 1/5*(colSums(y1[i:(i+4),])*sum(x1[i:(i+4)])))/(sum(x1[i:(i+4)]^2) - 1/5*(sum(x1[i:(i+4)]))^2)
    		varB <- rbind(varB,varb)
	}
	varA <- NULL
	varC <- NULL
	y3 <- NULL
	for (i in 2:8){
		yt <- colMeans(y2[i-1:i+1,])
		y3 <- rbind(y3,yt)
	}
	for (l in 1:length(kmin)){
		xp <- x1[which.min(varB[,l]):(which.min(varB[,l])+4)]
		#xp <- x1[(which.min((y2[3:6,l]-0.5)^2)):(which.min((y2[3:6,l]-0.5)^2)+4)]
		varC[l] <- min(varB[,l])
		#varC[l] <- varB[which.min((y2[3:6,l]-0.5)^2)+2,l]
    		yp <- t(MinGrowth(p[1:2],fmin4(p[3],kact_p[l,]*kmin_p[l,],nadph,knadph[l],kmin_p[l,],xp)))
    		vara <- (colSums(xp*yp) - 1/5*(colSums(yp)*sum(xp)))/(sum(xp^2) - 1/5*(sum(xp)^2))
    		varA <- cbind(varA,vara)
	}
	colnames(varA) <- colnames(y1)
	kmin_ps <- NULL
	kact_ps <- NULL
	fitbeta <- NULL
	for(i in 1:length(kcat)){
    		fitbeta[i] <- which.min((varA[,i] - varC[i])^2)
    		kmin_ps <- rbind(kmin_ps,kmin_p[i,fitbeta[i]])
    		kact_ps <- rbind(kact_ps,kact_p[i,fitbeta[i]])
	}

	tiff(file = paste(Sys.Date(),run,'solutions',sep = '_'), 
		width = 16.0, height = 8.7, units = 'cm', 
		res = 800, antialias = 'default')

	mutnames = names(grates)
	par(mfrow = c(2,4), oma = c(4,4,1,1),mar = c(0,0,1,0),xpd = NA)
	lbl = c('y','n','n','n','b','x','x','x')
	for (l in 1:length(kcat)){
 		plot(grates[,1],grates[,l+1], xlim = (c(0,grates[length(x1),1])), ylim = c(0,1),
         		labels = F, frame.plot = T, xlab = NA, ylab = NA, xaxp = c(0,80,4), pch = 20, cex = 0.66);
	    for(i in 1:length(kmin_p[l,])){
	      curve(MinGrowth_asin2(p[1:2],fmin3(p[3],kact_p[l,i]*kmin_p[l,i],nadph,knadph[l],kmin_p[l,i],x)),
        	    0, max(x1), add = T, lwd = 0.1)
	    }

	    curve(MinGrowth_asin2(p[1:2],fmin3(p[3],kact_ps[l]*kmin_ps[l],nadph,knadph[l],kmin_ps[l],x)),
        	  0, max(x1), add = T, lwd = 0.3, col = 2)
	    mtext(paste(mutnames[l+1],'  '), outer = F, side = 3, adj = 1, line = -1.5, cex = 0.6)
	    if (lbl[l] == 'x'){
	      axis(side = 1, tick = F, cin = 6/72)
	    } else if (lbl[l] == 'y'){
	      axis(side = 2, tick = F, at = c(0,0.5,1), cin = 6/72)
	    } else if (lbl[l] == 'b'){
	      axis(side = 1, tick = F)
	      axis(side = 2, tick = F, at = c(0,0.5,1), cin = 6/72)
	    }
	}
	mtext(paste('A)'), line = 1.2, outer = T, cin = 6/72, adj = 0)
	mtext('Minocycline (ug/ml)', outer = T, side = 1, line = 2.2, cin = 6/72)
	mtext('Growth Rate', outer = T, side = 2, line = 2.4)
	dev.off()

	tiff(file = paste(Sys.Date(),run,'ratios',sep = '_'), 
		width = 16.0, height = 8.7, units = 'cm', 
		res = 800, antialias = 'default')
	mutnames = names(grates)
	par(mfrow = c(2,4), oma = c(4,4,1,1),mar = c(0,0,1,0),xpd = NA)
	lbl = c('y','n','n','n','b','x','x','x')
	for (l in 1:length(kcat)){	
		plot(grates[,1],grates[,l+1], xlim = (c(0,grates[length(x1),1])), ylim = c(0,1),
			labels = F, frame.plot = T, xlab = NA, ylab = NA, xaxp = c(0,80,4), pch = 20, cex = 0.66, col = 3)
		for(i in 1:1000){
			curve(MinGrowth_asin2(p[1:2],fmin3(p[3],kact_ps[l]*kmin_ps[l]*(i/100),nadph,knadph[l],kmin_ps[l]*(i/100),x)),
		    		0, max(x1), add = T, lwd = 0.1)
		}
		curve(MinGrowth_asin2(p[1:2],fmin3(p[3],kact_ps[l]*kmin_ps[l],nadph,knadph[l],kmin_ps[l],x)),
			  0, max(x1), add = T, lwd = 0.3, col = 2)
		mtext(paste(mutnames[l+1],'  '), outer = F, side = 3, adj = 1, line = -1.5, cex = 0.6)
		if (lbl[l] == 'x'){
			axis(side = 1, tick = F, cin = 6/72)
		} else if (lbl[l] == 'y'){
			axis(side = 2, tick = F, at = c(0,0.5,1), cin = 6/72)
		} else if (lbl[l] == 'b'){
			axis(side = 1, tick = F)
			axis(side = 2, tick = F, at = c(0,0.5,1), cin = 6/72)
		}
	}
	mtext(paste('A)'), line = 1.2, outer = T, cin = 6/72, adj = 0)
	mtext('Minocycline (ug/ml)', outer = T, side = 1, line = 2.2, cin = 6/72)
	mtext('Growth Rate', outer = T, side = 2, line = 2.4)
	dev.off()

	kact_m <- kact_ps
	#kact_m <- apply(kact_p,1,mean)

	kmin_pd <- NULL
	kcat_pd <- NULL
	for(i in 1:8){
	  out3 = optim(fsum3,p = 50,hessian = 1);
	  print(c(out3$convergence))
	  kmin_pd <- rbind(kmin_pd,out3$par[1])
	  kcat_pd <- rbind(kcat_pd,out3$par[1]*kact_m[i])
	}

	tiff(file = paste(Sys.Date(),run,'solution-final',sep = '_'), 
		width = 16.0, height = 8.7, units = 'cm', 
		res = 800, antialias = 'default')

	mutnames = names(grates)
	par(mfrow = c(2,4), oma = c(4,4,1,1),mar = c(0,0,1,0),xpd = NA)
	lbl = c('y','n','n','n','b','x','x','x')
	for (l in 1:length(kcat)){
		plot(grates[,1],grates[,l+1], xlim = (c(0,grates[length(x1),1])), ylim = c(0,1),
			labels = F, frame.plot = T, xlab = NA, ylab = NA, xaxp = c(0,80,4), pch = 20, cex = 0.66);
		curve(MinGrowth_asin2(p[1:2],fmin3(p[3],kcat_pd[l],nadph,knadph[l],kmin_pd[l],x)),
			0, max(x1), add = T, lwd = 0.3, col = 2)
		mtext(paste(mutnames[l+1],'  '), outer = F, side = 3, adj = 1, line = -1.5, cex = 0.6)
		if (lbl[l] == 'x'){
		axis(side = 1, tick = F, cin = 6/72)
		} else if (lbl[l] == 'y'){
		axis(side = 2, tick = F, at = c(0,0.5,1), cin = 6/72)
		} else if (lbl[l] == 'b'){
		axis(side = 1, tick = F)
		axis(side = 2, tick = F, at = c(0,0.5,1), cin = 6/72)
		}
	}
	mtext(paste('A)'), line = 1.2, outer = T, cin = 6/72, adj = 0)
	mtext('Minocycline (ug/ml)', outer = T, side = 1, line = 2.2, cin = 6/72)
	mtext('Growth Rate', outer = T, side = 2, line = 2.4)
	dev.off()

	Results <- data.frame('KM MCN' = kmin_pd, 'Vmax' = kcat_pd, 'Activity' = kact_m)
	rownames(Results) <- mutnames[2:9]
	assign(paste('Results', run, sep = '_'),Results)

	write.table(Results, file = paste(Sys.Date(),run,'Results',sep = '_'), sep = ',', row.names = T, col.names = T)
}

print(Results_B1)
print(Results_B2)
print(Results_B3)
print(Results_B4)
print(Results_B5)
print(Results_B6)
print(Results_B7)
print(Results_B8)
print(parameters)

write.table(parameters, file = paste(Sys.Date(),'B','Parameters',sep = '_'), sep = ',', row.names = F, col.names = F)
