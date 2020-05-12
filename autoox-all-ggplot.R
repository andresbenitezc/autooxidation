  
####
# R script to create plots from output data 
# this R script does not require an input, it finds all of the matching files in the current 
# folder to plot and analyze.
# Author: Andres Benitez
####

files = list.files(path = "./", pattern = "autooxCurves-seqexp.csv")
raw.files = list.files(path = "./", pattern = "-raw.csv")

values <- data.frame(y0.mean = double(), y0.sd = double(),
	y1.mean = double(), y1.sd = double(), y2.mean = double(), y2.sd = double(), 
	k1.mean = double(), k1.sd = double(), k2.mean = double(), k2.sd = double(), 
	v0.mean = double(), v0.sd = double(), km.mean = double(), km.sd = double(), y1y2 = double())
labels <- names(values)

color = c("blue", "red", "green", "purple")

for (i in files){
	data = read.csv(i, header = TRUE)
	values = rbind(values, c(rbind(data$mean, data$stdev), data$mean[2]/data$mean[3]))
}

names(values) = labels
row.names(values) = sub("_autooxCurves-seqexp.csv", "", files)

values.n <- values

values.n$v0.sd = values.n$v0.sd/values.n$v0.mean[1]
values.n$v0.mean = values.n$v0.mean/values.n$v0.mean[1]

values.n$km.sd = values.n$km.sd/values.n$km.mean[1]
values.n$km.mean = values.n$km.mean/values.n$km.mean[1]

values.n$k1.sd = values.n$k1.sd/values.n$k1.mean[1]
values.n$k1.mean = values.n$k1.mean/values.n$k1.mean[1]

values.n$k2.sd = values.n$k2.sd/values.n$k2.mean[1]
values.n$k2.mean = values.n$k2.mean/values.n$k2.mean[1]

write.table(values, file = "autoox_all.csv", sep = ",", col.names = NA, row.names = TRUE)

data.total <- NULL
for ( i in raw.files){
	
	data.raw <- read.csv(i, header = TRUE) %>% gather("file", "HbOxy", -X) %>% rename("Time" = X) %>%
		group_by(file) %>% mutate(HbOxy = (HbOxy- min(HbOxy))/(HbOxy[1] - min(HbOxy))) %>% ungroup() %>%
		group_by(Time) %>% summarize(mean(HbOxy), sd(HbOxy)) %>% mutate(file = i)
	ifelse(is.null(data.total), data.total <- data.raw, data.total <- union(data.total, data.raw))
}

data.total <- data.total %>% mutate(file = sub("-raw.csv", "", file)) %>% 
	rename("Hb Oxy (normalized)" = `mean(HbOxy)`) %>%
	rename("SD" = `sd(HbOxy)`) %>% 
	rename("Time (min)" = Time) %>%
	rename("Hb Mutant" = file)

data.plot <- ggplot(data.total, aes(y = `Hb Oxy (normalized)`, x = `Time (min)`, color = `Hb Mutant`)) +
	geom_point() +
	geom_errorbar(aes(ymin = `Hb Oxy (normalized)` - SD, ymax = `Hb Oxy (normalized)` + SD), width = 50) +
	theme_bw()  
