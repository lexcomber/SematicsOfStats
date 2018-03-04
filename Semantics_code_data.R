##########################################################

# The varying semantics of statistical modelling: 
# the example of regression, LASSO, GWR and GW-LASSO

# Lex Comber
# March 2018
# email: a.comber@leeds.ac.uk

# The code describes all of the analyses undertaken in the paper 

##########################################################
library(GISTools)
library(rgdal)
library(raster)
library(reshape2)
library(tidyverse)
library(OpenStreetMap)
library(grid)
library(GWmodel)
library(gwrr)
library(ggthemes)
library(tmap)
library(glmnet)
library(repmis)

#### Part 1: Load Data 
source_data("https://github.com/lexcomber/SematicsOfStats/blob/master/data_in.RData?raw=True")

# get rid of anomly
## convert aspect to Northness and Eastness
# https://groups.google.com/forum/#!msg/Maxent/5VA0ATKJrPI/utSe0PaV2W4J
# aspect
a <- data.30$aspect
data.30$east <- sin(a)
data.30$north <- cos(a)
a <- data.90$aspect
data.90$east <- sin(a)
data.90$north <- cos(a)
a <- data.180$aspect
data.180$east <- sin(a)
data.180$north <- cos(a)
# flow direction
a <- ((data.30$flowdir/128) * 360) 
data.30$flowne <- sin(a)
a <- ((data.90$flowdir/128) * 360) 
data.90$flowne <- sin(a)
a <- ((data.180$flowdir/128) * 360) 
data.180$flowne <- sin(a)

data.30$Ks <-data.30$Ks.0.10cm.
data.90$Ks <-data.90$Ks.0.10cm.
data.180$Ks <-data.180$Ks.0.10cm.

#### Part 2: Descriptive statistics and plots
# see http://spatialreference.org/ref/?search=krasovsky+1940+Albers
tmp.proj <- CRS("+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +x_0=4000000 +y_0=0 +datum=WGS84 +units=m +no_defs")
source_data("https://github.com/lexcomber/SematicsOfStats/blob/master/LP_outline.RData?raw=True")

crs.val <- CRS("+proj=longlat +datum=WGS84 ")
out.ll <- spTransform(out, crs.val)
data.ll <- spTransform(data.sp, crs.val)
# Study Area Figure
ul <- as.vector(cbind(bbox(out.ll)[2,2], 
                      bbox(out.ll)[1,1]))
lr <- as.vector(cbind(bbox(out.ll)[2,1], 
                      bbox(out.ll)[1,2]))
MyMap <- openmap(ul,lr, type = "stamen-terrain")
### Figure 1
png("fig1.png", width = 10, height = 6, units = 'in', res = 300)
par(mar = c(0,0,0,0))
plot(MyMap, removeMargin=FALSE)
plot(spTransform(out.ll, osm()), add = TRUE, lwd = 2)
plot(spTransform(data.ll, osm()),bg = "#25252580", 
     add = TRUE, pch = 21, cex = 1.5)
scalebar(200000, type = "bar", xy = c(11329039, 4081403), 
         divs = 4, below = "km",
         label = c(0,100,200) )
dev.off()

### Figure 2: boxplots
## ggplot prep
index <- c(4:6,10:12, 19:21, 22:24)
tmp <- data.30@data[, index]
depths <- c("0 to 10", "10 to 20", "20 to 40")
depth.i <- rep(1, nrow(tmp))
m<- melt(tmp)
head(m)
depth.i <- rep(NA, nrow(m))
i <- grep("0.10", m$variable)
depth.i[i] <- depths[1]
i <- grep("0.20", m$variable)
depth.i[i] <- depths[2]
i <- grep("0.40", m$variable)
depth.i[i] <- depths[3]
m <- data.frame(m, depth = depth.i)
m$variable <- gsub(".0.10cm", "",m$variable )
m$variable <- gsub("0.10cm", "",m$variable )
m$variable <- gsub(".10.20cm", "",m$variable )
m$variable <- gsub("10.20cm", "",m$variable )
m$variable <- gsub(".20.40cm", "",m$variable )
m$variable <- gsub("20.40cm", "",m$variable )
m$variable <- gsub("10.20", "",m$variable )
m$variable <- gsub("20.40", "",m$variable )
m$variable <- gsub("SSWC.", "SSWC",m$variable )
m$variable <- gsub("Sand", "San",m$variable )
m$variable <- gsub("San", "Sand",m$variable )
m$variable <- gsub("Clay", "Cla",m$variable )
m$variable <- gsub("Cla", "Clay",m$variable )
m$variable <- gsub("FC.", "FC",m$variable )

# check 
table(m$variable)
png("fig2.png", width = 8, height = 3, units = 'in', res = 300)
ggplot(m, aes(
      y = as.numeric(value),
      x = as.factor(depth),
      fill = as.factor(depth))
  )+ 
  geom_boxplot(show.legend=FALSE) +
  scale_y_continuous(limits = quantile(m$value, c(0.25, 0.75)))+
  coord_flip() +
  facet_wrap(~variable, scales = "free",ncol = 4) +
  scale_fill_brewer(type="qual", palette = "Set2") +
  ylab("Variable Value") +
  xlab("Measurement Depth") +
  theme_bw()
dev.off()

### Figure 3: density plots
# function
hist.func.dem <- function(name = "aspect", tit = name, bw = 40) {
  data.tmp <- data.frame(data.30@data[, c(name)], data.90@data[, c(name)], 
                         data.180@data[, c(name)])
  names(data.tmp) <- c("30 m", "90 m", "180 m")
  data.tmp <- melt(data.tmp)
  p <- ggplot(data.tmp, aes(x=value)) + 
    geom_histogram(aes(y=..density..),
                   binwidth=bw,colour="white") +
    geom_density(alpha=.2, fill="#FF6666") +
    facet_wrap( ~ variable, ncol = 3) +
    #theme(axis.text.y = element_text(size=2))+
    theme(axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    #ylab("Probability Density") +
    #coord_flip() +
    xlab(tit)
  return(p)
}  

# c("Intercept", "Clay","Sand","FC","SSWC","TRI","Slope","East-ness","North-ness","Flow-NE")
p1 <- hist.func.dem("north", "North-ness", 0.2)
p2 <- hist.func.dem("east", "East-ness", 0.2)
p3 <- hist.func.dem("flowne", "Flow-NE", 0.2)
p4 <- hist.func.dem("slope", "Slope", bw = 5)
p5 <- hist.func.dem("tri", "Terrain Ruggedness Index", bw = 5)
### Figure 3
png("fig3.png", width = 8, height = 6, units = 'in', res = 300)
# open a new plot page
grid.newpage()
# set up the layout
pushViewport(viewport(layout=grid.layout(2,2)))
# plot using he print command
print(p4, vp=viewport(layout.pos.col = 1, layout.pos.row = 1, height = 5))
print(p3, vp=viewport(layout.pos.col = 1, layout.pos.row = 2, height = 5))
print(p1, vp=viewport(layout.pos.col = 2, layout.pos.row = 1, height = 5))
print(p5, vp=viewport(layout.pos.col = 2, layout.pos.row = 2, height = 5))
dev.off()

#### Part 3: Analysis - LM, LASSO, GWR and GWL
# see new_anal.R for initial model selection
# 3 models 
reg.mod.10 <- as.formula(Ks~Clay0.10cm+Sand0.10cm+FC.0.10cm.+SSWC.0.10cm.+tri+slope+east+north+flowne)
reg.mod.20 <- as.formula(Ks~Clay10.20+Sand10.20+FC.10.20cm.+SSWC.10.20cm.+tri+slope+east+north+flowne)              
reg.mod.40 <- as.formula(Ks~Clay20.40+Sand20.40+FC.20.40cm.+SSWC.20.40cm.+tri+slope+east+north+flowne)  

### 9 models of each type
# 3.1 lm
lm.10.r30 <- lm(reg.mod.10, data = data.30)                       
lm.10.r90 <- lm(reg.mod.10, data = data.90)                        
lm.10.r180 <- lm(reg.mod.10, data = data.180)  
lm.20.r30 <- lm(reg.mod.20, data = data.30)                       
lm.20.r90 <- lm(reg.mod.20, data = data.90)                        
lm.20.r180 <- lm(reg.mod.20, data = data.180)  
lm.40.r30 <- lm(reg.mod.40, data = data.30)                       
lm.40.r90 <- lm(reg.mod.40, data = data.90)                        
lm.40.r180 <- lm(reg.mod.40, data = data.180)  

# 3.2 lasso
x.10.names <- c("Clay0.10cm","Sand0.10cm","FC.0.10cm.","SSWC.0.10cm.","tri","slope","east","north","flowne")
x.20.names <- c("Clay10.20","Sand10.20","FC.10.20cm.","SSWC.10.20cm.","tri","slope","east","north","flowne")
x.40.names <- c("Clay20.40","Sand20.40","FC.20.40cm.","SSWC.20.40cm.","tri","slope","east","north","flowne")
y.names <- "Ks"
eln.func <- function(y.names, x.names, data = data.30@data) {
	x <- as.matrix(data[x.names])
	y <- as.matrix(data[y.names])
	cvfit = cv.glmnet(x, y)
	r2 <- cvfit$glmnet.fit$dev.ratio[which(cvfit$glmnet.fit$lambda == cvfit$lambda.min)]
	yhat <- as.vector(predict(cvfit,newx=x,s="lambda.min"))
	return(list(coef(cvfit, s = "lambda.min"), r2, yhat))
}
eln.10.r30 <- eln.func(y.names, x.10.names, data = data.30@data) 
eln.10.r90 <- eln.func(y.names, x.10.names, data = data.90@data) 
eln.10.r180 <- eln.func(y.names, x.10.names, data = data.180@data) 
eln.20.r30 <- eln.func(y.names, x.20.names, data = data.30@data) 
eln.20.r90 <- eln.func(y.names, x.20.names, data = data.90@data) 
eln.20.r180 <- eln.func(y.names, x.20.names, data = data.180@data) 
eln.40.r30 <- eln.func(y.names, x.40.names, data = data.30@data) 
eln.40.r90 <- eln.func(y.names, x.40.names, data = data.90@data) 
eln.40.r180 <- eln.func(y.names, x.40.names, data = data.180@data) 

# 3.3 GWR
bw.gwr.10.30 <- bw.gwr(as.formula(reg.mod.10), data = data.30, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwr.10.r30 <- gwr.basic(as.formula(reg.mod.10), data = data.30, bw = bw.gwr.10.30, 
                  kernel = "bisquare", adaptive = T)
bw.gwr.20.30 <- bw.gwr(as.formula(reg.mod.20), data = data.30, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwr.20.r30 <- gwr.basic(as.formula(reg.mod.20), data = data.30, bw = bw.gwr.20.30, 
                  kernel = "bisquare", adaptive = T)
bw.gwr.40.30 <- bw.gwr(as.formula(reg.mod.40), data = data.30, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwr.40.r30 <- gwr.basic(as.formula(reg.mod.40), data = data.30, bw = bw.gwr.40.30, 
                  kernel = "bisquare", adaptive = T)

bw.gwr.10.90 <- bw.gwr(as.formula(reg.mod.10), data = data.90, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwr.10.r90 <- gwr.basic(as.formula(reg.mod.10), data = data.90, bw = bw.gwr.10.90, 
                  kernel = "bisquare", adaptive = T)
bw.gwr.20.90 <- bw.gwr(as.formula(reg.mod.20), data = data.90, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwr.20.r90 <- gwr.basic(as.formula(reg.mod.20), data = data.90, bw = bw.gwr.20.90 , 
                  kernel = "bisquare", adaptive = T)
bw.gwr.40.90 <- bw.gwr(as.formula(reg.mod.40), data = data.90, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwr.40.r90 <- gwr.basic(as.formula(reg.mod.40), data = data.90, bw = bw.gwr.40.90, 
                  kernel = "bisquare", adaptive = T)

bw.gwr.10.180 <- bw.gwr(as.formula(reg.mod.10), data = data.180, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwr.10.r180 <- gwr.basic(as.formula(reg.mod.10), data = data.180, bw = bw.gwr.10.180, 
                  kernel = "bisquare", adaptive = T)
bw.gwr.20.180 <- bw.gwr(as.formula(reg.mod.20), data = data.180, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwr.20.r180 <- gwr.basic(as.formula(reg.mod.20), data = data.180, bw = bw.gwr.20.180, 
                  kernel = "bisquare", adaptive = T)
bw.gwr.40.180 <- bw.gwr(as.formula(reg.mod.40), data = data.180, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwr.40.r180 <- gwr.basic(as.formula(reg.mod.40), data = data.180, bw = bw.gwr.40.180, 
                  kernel = "bisquare", adaptive = T)

# 3.4 GWR-L
source("funcs.R") # needs to be modified to do adapt = F, aic, lamda.2se
dMat <- gw.dist(coordinates(data.30))
bw.gwl.10.r30 <- bw.lass(as.formula(reg.mod.10), data = data.30, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwl.10.r30 <- gwr.lass(reg.mod.10, data = data.30, regression.points = data.30, 
			bw = bw.gwl.10.r30, kernel = "bisquare", adaptive = T, p = 2, 
			theta = 0, longlat = F, 
			dMat = dMat, F123.test = F, cv = T, W.vect = NULL) 
bw.gwl.20.r30 <- bw.lass(as.formula(reg.mod.20), data = data.30, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwl.20.r30 <- gwr.lass(reg.mod.20, data = data.30, regression.points = data.30, 
			bw = bw.gwl.20.r30, kernel = "bisquare", adaptive = T, p = 2, 
			theta = 0, longlat = F, 
			dMat = dMat, F123.test = F, cv = T, W.vect = NULL) 
bw.gwl.40.r30 <- bw.lass(as.formula(reg.mod.40), data = data.30, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwl.40.r30 <- gwr.lass(reg.mod.40, data = data.30, regression.points = data.30, 
			bw = bw.gwl.40.r30, kernel = "bisquare", adaptive = T, p = 2, 
			theta = 0, longlat = F, 
			dMat = dMat, F123.test = F, cv = T, W.vect = NULL) 

bw.gwl.10.r90 <- bw.lass(as.formula(reg.mod.10), data = data.90, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwl.10.r90 <- gwr.lass(reg.mod.10, data = data.30, regression.points = data.90, 
			bw = bw.gwl.10.r90, kernel = "bisquare", adaptive = T, p = 2, 
			theta = 0, longlat = F, 
			dMat = dMat, F123.test = F, cv = T, W.vect = NULL) 
bw.gwl.20.r90 <- bw.lass(as.formula(reg.mod.20), data = data.90, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwl.20.r90 <- gwr.lass(reg.mod.20, data = data.30, regression.points = data.90, 
			bw = bw.gwl.20.r90, kernel = "bisquare", adaptive = T, p = 2, 
			theta = 0, longlat = F, 
			dMat = dMat, F123.test = F, cv = T, W.vect = NULL) 
bw.gwl.40.r90 <- bw.lass(as.formula(reg.mod.40), data = data.90, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwl.40.r90 <- gwr.lass(reg.mod.40, data = data.30, regression.points = data.90, 
			bw = bw.gwl.40.r90, kernel = "bisquare", adaptive = T, p = 2, 
			theta = 0, longlat = F, 
			dMat = dMat, F123.test = F, cv = T, W.vect = NULL) 

bw.gwl.10.r180 <- bw.lass(as.formula(reg.mod.10), data = data.180, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwl.10.r180 <- gwr.lass(reg.mod.10, data = data.30, regression.points = data.180, 
			bw = bw.gwl.10.r180, kernel = "bisquare", adaptive = T, p = 2, 
			theta = 0, longlat = F, 
			dMat = dMat, F123.test = F, cv = T, W.vect = NULL) 
bw.gwl.20.r180 <- bw.lass(as.formula(reg.mod.20), data = data.180, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwl.20.r180 <- gwr.lass(reg.mod.20, data = data.30, regression.points = data.180, 
			bw = bw.gwl.20.r180, kernel = "bisquare", adaptive = T, p = 2, 
			theta = 0, longlat = F, 
			dMat = dMat, F123.test = F, cv = T, W.vect = NULL) 
bw.gwl.40.r180 <- bw.lass(as.formula(reg.mod.40), data = data.180, 
                kernel = "gaussian", adaptive = T, approach = "CV") 
gwl.40.r180 <- gwr.lass(reg.mod.40, data = data.30, regression.points = data.180, 
			bw = bw.gwl.40.r180, kernel = "bisquare", adaptive = T, p = 2, 
			theta = 0, longlat = F, 
			dMat = dMat, F123.test = F, cv = T, W.vect = NULL) 

save.image("new_anal2.RData")

#### Part 4 Model summaries and visiualisation 
lm.list <- ls(pattern = "lm.")
las.list <- ls(pattern = "eln")[1:9]
gwr.list <- ls(pattern = "gwr.")[10:18]
bwgwr.list <- ls(pattern = "gwr.")[1:9]
gwl.list <- ls(pattern = "gwl.")[10:18]
bwgwl.list <- ls(pattern = "gwl.")[1:9]

vn <- c("Intercept", "Clay","Sand","FC","SSWC",
		"TRI","Slope","East-ness","North-ness","Flow-NE")

#### 4.1 Tables
## OLS
tab <- matrix(ncol = 18, nrow = 11)
lm.lis <- lm.list[c(2,3,1,5,6,4,8,9,7)]
for (i in 1: length(lm.lis)) {
    lm.i <- get(lm.lis[i])
    ii <- i *2
    tab[1:10,c((ii-1):ii)] <- summary(lm.i)$coefficients[, c(1,4)]
    tab[11,(ii-1)]<- unlist(summary(lm.i)[8])
}
cn <- c(paste0(lm.lis, " Coef"), paste0(lm.lis, " P"))
colnames(tab) <- cn
rownames(tab) <- c(vn, "R2")
write.csv(tab, "tab1.csv")

## Lasso
r2.function <- function(y, yhat) {
	R2 <- 1 - (sum((y-yhat )^2)/sum((y-mean(y))^2))
	#list(R2, RMSE)
	R2
}
y = data.30$Ks
tab <- matrix(nrow = 11, ncol = 9)
lm.lis <- las.list[c(2,3,1,5,6,4,8,9,7)]
for (i in 1: length(lm.lis)){
	l.i <- get(lm.lis[i])
	coef <- as.vector(l.i[[1]])
	tab[1:10,i] <- coef
   	tab[11, i] <- r2.function(y,l.i[[3]])
}
cn <- c(paste0(lm.lis, " Coef"))
colnames(tab) <- cn
rownames(tab) <- c(vn, "R2")
write.csv(tab, "tab2.csv")
round(tab,3)

## GWR
tab <- matrix(nrow = 12, ncol = 27)
lm.lis <- gwr.list[c(2,3,1,5,6,4,8,9,7)]
bw.lis <- bwgwr.list[c(2,3,1,5,6,4,8,9,7)]
for(i in 1:length(lm.lis)) {
	ii <- 3 * i
	gw.i <- get(lm.lis[i])
	tab[1:10,(ii-2):ii] <- (t(apply(gw.i$SDF@data[,1:10], 2, summary))[, c(2,3,5)])
	tab[11, (ii-1)] <- r2.function(y,gw.i$SDF$yhat)
	tab[12, (ii-1)] <- get(bw.lis[i])
}
cn <- c(paste0(lm.lis, " Q1"), paste0(lm.lis, " Q2"), paste0(lm.lis, " Q3"))
cn <- cn[c(1,10, 19, 2,11, 20, 3,12,21, 4,13,22, 5,14,23, 6,15,24, 7,16,25, 8,17, 26, 9,18, 27)]
colnames(tab) <- cn
rownames(tab) <- c(vn, "R2", "BW")
write.csv(t(tab), "tab3.csv")
round(t(tab),3)

## GWL
tab <- matrix(nrow = 12, ncol = 27)
lm.lis <- gwl.list[c(2,3,1,5,6,4,8,9,7)]
bw.lis <- bwgwl.list[c(2,3,1,5,6,4,8,9,7)]
for(i in 1:length(lm.lis)) {
	ii <- 3 * i
	gl.i <- get(lm.lis[i])
	tab[1:10,(ii-2):ii] <- (t(apply(gl.i$SDF@data[,1:10], 2, summary))[, c(2,3,5)])
	tab[11, (ii-1)] <- r2.function(y,gl.i$yhat)
	tab[12, (ii-1)] <- get(bw.lis[i])

}
cn <- c(paste0(lm.lis, " Q1"), paste0(lm.lis, " Q2"), paste0(lm.lis, " Q3"))
cn <- cn[c(1,10, 19, 2,11, 20, 3,12,21, 4,13,22, 5,14,23, 6,15,24, 7,16,25, 8,17, 26, 9,18, 27)]
colnames(tab) <- cn
rownames(tab) <- c(vn, "R2", "BW")
write.csv(t(tab), "tab4.csv")
round(t(tab),3)

#### 4.2 Plots
## Dot Plots of Coefficients 
cn <- c("Coef", "lower", "upper")

do.3val.plot <- function(df) {
	b1_sign <- c("#7b3294","#008837","grey")
	#b1_sign <- c("#7b3294","#008837")
	p <- df %>% 
	    #mutate(sign=ifelse(df$Coef>0,"pos","neg")) %>%
	    mutate(sign=ifelse(df$Coef==0,"zero", ifelse(df$Coef<0,"neg","pos"))) %>%
	      ggplot(aes(y=Coef, x=Variable)) +
	      geom_pointrange(aes(ymin=lower, ymax=upper, colour = factor(sign)),
	                size=0.6,                    
	                position=position_dodge(.9), 
	                show.legend = F) +
	      facet_wrap(facets = "method", nrow = 1) +        
	      geom_hline(yintercept=0, linetype="dashed")+
	      coord_flip() +
	      theme_bw() +
		  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
	      #scale_x_discrete(limits = order_expl$row )+
	      #scale_y_continuous(limits=c(-1.5,1.5))+
	      scale_colour_manual(values = b1_sign )+
	      labs(x="",y="")
      	  #theme_dark()
	return(p)
}
#do.3val.plot(df)

for(i in 1:length(lm.list)) {

	l.i <- get(lm.list[i])
	ls.i <- get(las.list[i])	
	gw.i <- get(gwr.list[i])
	gl.i <- get(gwl.list[i])

	df.1 <- data.frame(
		OLS = coef(l.i))
	df.1 <-cbind(df.1, confint(l.i, level=0.95))
	names(df.1)  <- cn
	df.1$method <- "OLS"
	df.1$Variable <- vn	#df.1 <- unname(df.4)
	rownames(df.1) <- NULL
	df.2 <- data.frame(
		Lasso = as.vector(ls.i[[1]]), 
		lower = rep(0, length(as.vector(ls.i[[1]]))),
		upper = rep(0, length(as.vector(ls.i[[1]]))))
	names(df.2)  <- cn
	df.2$method <- "Lasso"
	df.2$Variable <- vn
	rownames(df.2) <- NULL
	df.3 <- data.frame(t(apply(gw.i$SDF@data[,1:10], 2, summary))[, c(3,2,5)])
	names(df.3)  <- cn
	df.3$method <- "GWR"
	df.3$Variable <- vn	
	rownames(df.3) <- NULL	
	df.4 <- data.frame(t(apply(gl.i$SDF@data, 2, summary))[, c(3,2,5)])
	names(df.4)  <- cn
	df.4$method <- "GWL"
	df.4$Variable <- vn
	rownames(df.4) <- NULL
	df <- rbind(df.1,df.2,df.3,df.4)
	df$Variable <- factor(df$Variable, levels = rev(vn))
	df$method <- factor(df$method, levels = (c("OLS", "Lasso", "GWR", "GWL")))
	#df[,1:3] <- as_tibble(t(apply(df[,1:3], 1, function(x) scale(x, center = 0))))
	#df[is.na(df$Coef), 1:3] <- c(0,0,0)
	df <- df[df$Variable != "Intercept", ]
	p <- do.3val.plot(df)
	tit <- lm.list[i]
	tit <- gsub("lm.", "Plot", tit)
	assign(tit, p)
	#tab <- dcast(melt(df), Variable~variable+method, mean)
}
plot.list <- ls(pattern = "Plot")[c(2,3,1,5,6,4,8,9,7)]

png("fig4.png", width = 12, height = 6, units = 'in', res = 300)
# open a new plot page
grid.newpage()
# set up the layout
pushViewport(viewport(layout=grid.layout(3,3)))
# plot using he print command
print(get(plot.list[1]), vp=viewport(layout.pos.col = 1, layout.pos.row = 1, height = 5))
print(get(plot.list[2]), vp=viewport(layout.pos.col = 1, layout.pos.row = 2, height = 5))
print(get(plot.list[3]), vp=viewport(layout.pos.col = 1, layout.pos.row = 3, height = 5))
print(get(plot.list[4]), vp=viewport(layout.pos.col = 2, layout.pos.row = 1, height = 5))
print(get(plot.list[5]), vp=viewport(layout.pos.col = 2, layout.pos.row = 2, height = 5))
print(get(plot.list[6]), vp=viewport(layout.pos.col = 2, layout.pos.row = 3, height = 5))
print(get(plot.list[7]), vp=viewport(layout.pos.col = 3, layout.pos.row = 1, height = 5))
print(get(plot.list[8]), vp=viewport(layout.pos.col = 3, layout.pos.row = 2, height = 5))
print(get(plot.list[9]), vp=viewport(layout.pos.col = 3, layout.pos.row = 3, height = 5))
dev.off()

# do pred v obs
for(i in 1:length(lm.list)) {

	l.i <- get(lm.list[i])
	ls.i <- get(las.list[i])	
	gw.i <- get(gwr.list[i])
	gl.i <- get(gwl.list[i])
	y = data.30$Ks
	df.1 <- data.frame(
				y,
				yhat = as.vector(l.i$fitted.values)) 
	df.1$method <- "OLS"
	df.2 <- data.frame(
				y,
				yhat = ls.i[[3]])
	df.2$method <- "Lasso"
	df.3 <- data.frame(
				y,
				yhat = gw.i$SDF$yhat)
	df.3$method <- "GWR"
	df.4 <- data.frame(
				y,
				yhat = gl.i$yhat)
	df.4$method <- "GWL"
	df2 <- rbind(df.1, df.2,df.3,df.4)
	r2 <- data.frame(
				OLS = summary(l.i)$r.squared,
				Lasso = ls.i[[2]],
				GWR = r2.function(y,gw.i$SDF$yhat),
				GWL = r2.function(y,gl.i$yhat))
	tits <- paste0(c("OLS:", "Lasso:", "GWR:", "GWL:"), round(r2,3))
	df2$method <- factor(df2$method, levels = (c("OLS", "Lasso", "GWR", "GWL")))	
	levels(df2$method) <- tits
	
	p <- ggplot(data = df2, aes(y=y,x=yhat)) +
	  geom_point(alpha = 0.5, shape = 16) +
	  facet_wrap( ~ method, nrow = 1, scales = "fixed") +
	  geom_smooth(method='lm', col = "#A50F15") +
	  xlab("Predicted Ks") +
	  ylab("Observed Ks") +
	  theme_bw() +
	  theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
	  	axis.text.y = element_blank())
	     
	tit <- lm.list[i]
	tit <- gsub("lm.", "Scat", tit)
	assign(tit, p)
}
scatt.list <- ls(pattern = "Scat")[c(2,3,1,5,6,4,8,9,7)]


png("fig5.png", width = 12, height = 6, units = 'in', res = 300)
# open a new plot page
grid.newpage()
# set up the layout
pushViewport(viewport(layout=grid.layout(3,3)))
# plot using he print command
print(get(scatt.list[1]), vp=viewport(layout.pos.col = 1, layout.pos.row = 1, height = 5))
print(get(scatt.list[2]), vp=viewport(layout.pos.col = 1, layout.pos.row = 2, height = 5))
print(get(scatt.list[3]), vp=viewport(layout.pos.col = 1, layout.pos.row = 3, height = 5))
print(get(scatt.list[4]), vp=viewport(layout.pos.col = 2, layout.pos.row = 1, height = 5))
print(get(scatt.list[5]), vp=viewport(layout.pos.col = 2, layout.pos.row = 2, height = 5))
print(get(scatt.list[6]), vp=viewport(layout.pos.col = 2, layout.pos.row = 3, height = 5))
print(get(scatt.list[7]), vp=viewport(layout.pos.col = 3, layout.pos.row = 1, height = 5))
print(get(scatt.list[8]), vp=viewport(layout.pos.col = 3, layout.pos.row = 2, height = 5))
print(get(scatt.list[9]), vp=viewport(layout.pos.col = 3, layout.pos.row = 3, height = 5))
dev.off()

######################### END ###########################