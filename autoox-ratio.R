####
# R script for calculating ratios of species from spectrums.
# This R script requires a csv file with the spectrums that have to be read and the reference curves for each species.
# The script also takes an imput for a file name to output the data table and plot to.
# Author: Andres Benitez
# Cite: Benitez Cardenas AS, Samuel PP, Olson JS. Current Challenges in the Development of Acellular Hemoglobin Oxygen Carriers by Protein Engineering. Shock. 2019;52(1S Suppl 1):28‐40. doi:10.1097/SHK.0000000000001053
####

library(scatterplot3d)
library(minpack.lm)
library(deconspectra)
library(signal)

file <- commandArgs(trailingOnly = TRUE)

data.raw <- read.csv(file[1], header = TRUE)

n = which(names(data.raw) == "X0.0")
abs_curves <- data.frame(data.raw[-(1:(n-1))])
wavelength <- data.raw[[1]]
ref_curves <- data.frame(data.raw[2:(n-1)])

abs_curves <- data.frame(apply(t(abs_curves),1,function(x){sgolayfilt(x, p = 7, n = 11)}))

#ref_curves <- data.frame(lapply(ref_curves, function(x){ x*ifelse(wavelength >= 500 & wavelength <= 600, 4, 1) }))
#abs_curves <- data.frame(lapply(abs_curves, function(x){ x*ifelse(wavelength >= 500 & wavelength <= 600, 4, 1) }))

nls.fit <- deconSpectra(abs_curves, wavelength, ref_curves)

all.time <- rep(nls.fit$time, each = length(nls.fit$wavelength))
all.wavelength <- rep(nls.fit$wavelength, times = length(nls.fit$time))
all.residuals <- c(residuals(nls.fit))
all.oxy <- as.vector(ref_curves$oxy %o% nls.fit$values$oxy)
all.fer <- as.vector(ref_curves$fer %o% nls.fit$values$fer)
all.hch <- as.vector(ref_curves$hch %o% nls.fit$values$hch)
all.fyl <- as.vector(ref_curves$fyl %o% nls.fit$values$fyl)
all.hem <- as.vector(ref_curves$hem %o% nls.fit$values$hem)
all.dox <- as.vector(ref_curves$dox %o% nls.fit$values$dox)

pdf(file = sub("csv", "pdf", file[1]), width = 8, height = 5)
plot(nls.fit, normalized = TRUE, xlab = "time (min)", ylab = "relative concentration")
plot(nls.fit, normalized = FALSE, xlab = "time (min)", ylab = "concentration relative to reference" )
scatterplot3d(all.wavelength, all.time, all.residuals, pch = ".")

dataout <- cbind(time = nls.fit$time, nls.fit$values)

write.csv(dataout, file = sub(".csv", " traces.csv", file[1]), row.names = FALSE)
