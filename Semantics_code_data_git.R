#### Part 2: Descriptive statistics and plots

#RS imagery # ls5_20030924_30m – cloud free and nearest to 2007# ls5_20110930_30m – cloud free and nearest to 2012# .img files created and clipped using Norfolk_clip.shp# converted to IDRISI raster files in rs_data folder# 6 bands for each date
# from https://landsat.usgs.gov/what-are-band-designations-landsat-satellites
# Band 1 - Blue	0.45-0.52	30
# Band 2 - Green	0.52-0.60	30
# Band 3 - Red	0.63-0.69	30
# Band 4 - Near Infrared (NIR)	0.76-0.90	30
# Band 5  - Shortwave Infrared (SWIR) 1	1.55-1.75	30
# Band 6 - Thermal	10.40-12.50	120* (30)

##########################################################
library(GISTools)
library(rgdal)
library(raster)
library(reshape2)
library(tidyverse)
library(OpenStreetMap)
library(grid)
library(GWmodel)
library(tmap)
library(nlme)
library(repmis)

###### START PART 1 - Demonstration #######

### Data set up
source_data("https://github.com/lexcomber/regression_semantics/blob/master/rc.RData?raw=True")
new.proj <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.999601272 +x_0=400000 +y_0=-100000
+ +ellps=airy +towgs84=375,-111,431,0,0,0,0 +units=m +no_defs")

#y  <- as(readGDAL("2011_all/2011_b1.rst"), "RasterLayer"); proj4string(y)  <- new.proj
#x1 <- as(readGDAL("2011_all/2011_b2.rst"), "RasterLayer"); proj4string(x1) <- new.proj
#x2 <- as(readGDAL("2011_all/2011_b3.rst"), "RasterLayer"); proj4string(x2) <- new.proj
#x3 <- as(readGDAL("2011_all/2011_b4.rst"), "RasterLayer"); proj4string(x3) <- new.proj

#r <- stack(y,x1,x2,x3)
#plotRGB(stack(r), 1, 2, 3, stretch="lin")
#scalebar(10000, type = "bar", xy = c(616074.1, 291496.7), 
#	divs = 2, col = "white", below = "m")

#ex <- extent(604000, 619000, 290000, 305000)
#ex <- extent(601000, 604000, 290000, 293000)

#rc <- crop(r, ex) 
#names(rc) <- c("y", "x1", "x2","x3")

# save(rc, file = "rc.RData")

ex <- extent(601000, 602500, 290000, 291500)
rc <- crop(rc, ex) 

#png(filename = "F1.png", w = 8, h = 8, units = "in", res = 600)
plotRGB(stack(rc), 3,2,1, stretch="lin")
scalebar(500, type = "bar", xy = c(601874.1, 290096.7), divs = 4, col = "white", below = "m")
#dev.off()
df <- data.frame(getValues(rc))

##### Analyses
ncols <- 5
### 1. Standard: y, y*, x1 where y predicted using a simple linear regression model only.
rm.standard <- as.formula(y~x1+0)
lm.stand <- lm(rm.standard, data = df)
summary(lm.stand)
ry <- raster(rc); ry  <- setValues(ry,  df$y)	
ryh<- raster(rc); ryh <- setValues(ryh, fitted(lm.stand))	
rx1<- raster(rc); rx1 <- setValues(rx1, df$x1)
rxe<- raster(rc); rxe <- setValues(rxe, residuals(lm.stand))	
r.stand <- brick(ry,ryh,rx1, rxe)

# Plot Functions
plot.func.fixed <- function(r = r.stand, var = "layer.1", tit = "y", breaks = NULL) {
	p <- tm_shape(r)+
		tm_raster(col = var, alpha = 1, 
			breaks = breaks,
     		title = tit) +
     tm_layout(
     	legend.title.size=1,
     	legend.text.size = 0.6,
     	legend.bg.color = "white",
     	legend.bg.alpha = 1) 
     return(p)
}
plot.func <- function(r = r.stand, var = "layer.1", tit = "y") {
	p <- tm_shape(r)+
		tm_raster(col = var, alpha = 1, 
			style = "kmeans", n = 6, 
     		title = tit) +
     tm_layout(
     	legend.title.size=1,
     	legend.text.size = 0.6,
     	legend.bg.color = "white",
     	legend.bg.alpha = 1) 
     return(p)
}
## Fixed breaks for legend
b.error <- c(-3300, -200, -100, 100, 200, 500, 2500)
b.ystar <- c(500, 900, 1300, 1600, 2000, 3000, 4600)

tits <- c("y", "y*", expression(paste("x"[1])), expression(paste(epsilon)))
for (i in 1:dim(r.stand)[3]) {
	if (i == 2) p <- plot.func.fixed(r.stand, var = names(r.stand)[i], tit = tits[i], 
		breaks = b.ystar)
	if (i == 4) p <- plot.func.fixed(r.stand, var = names(r.stand)[i], tit = tits[i], 
		breaks = b.error)
	if( i==1 | i == 3 | i == 5) 
		p <- plot.func(r.stand, var = names(r.stand)[i], tit = tits[i]) 
	assign(paste0("p",i), p)
}
#png(filename = "F2a.png", w = 9*1.5, h = 2.3*1.2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(1,5)))
print(p1, vp=viewport(layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.col = 5))
#dev.off()

### 2. x issue: y, y*, x1 * x2: the number kx and choice cx 
rm.x <- as.formula(y~x1+x2+0)
lm.x <- lm(rm.x, data = df)
summary(lm.x)
ry <- raster(rc); ry  <- setValues(ry,  df$y)	
ryh<- raster(rc); ryh <- setValues(ryh, fitted(lm.x))	
rx1<- raster(rc); rx1 <- setValues(rx1, df$x1)	
rx2<- raster(rc); rx2 <- setValues(rx2, df$x2)	
rxe<- raster(rc); rxe <- setValues(rxe, residuals(lm.x))	
r.x <- brick(ry,ryh,rx1,rx2,rxe)

tits = c("y", "y*", expression(paste("x"[1])), 
     		expression(paste("x"[2])), expression(paste(epsilon)) )
for (i in 1:dim(r.x)[3]) {
	if (i == 2) p <- plot.func.fixed(r.x, var = names(r.x)[i], tit = tits[i], 
		breaks = b.ystar)
	if (i == 5) p <- plot.func.fixed(r.x, var = names(r.x)[i], tit = tits[i], 
		breaks = b.error)
	if( i==1 | i == 3 | i == 4) 
		p <- plot.func(r.x, var = names(r.x)[i], tit = tits[i]) 
	assign(paste0("p",i), p)
}
#png(filename = "F2b.png", w = 9*1.5, h = 2.3*1.2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(1,5)))
print(p1, vp=viewport(layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.col = 4))
print(p5, vp=viewport(layout.pos.col = 5))
#dev.off()

### 3. v issue aggregate the x1 data 
# double the pixel size and run the regression again using NN y and predict 
# get different support y*
rca <- aggregate(rc, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
rca <- resample(rc, rca, "ngb")
rcb <- disaggregate(rca, fact = 2)
dfa <- data.frame(getValues(rcb))
dfa$y <- df$y

lm.v <- lm(rm.standard, data = dfa)
summary(lm.v)
ry <- raster(rcb); ry  <- setValues(ry,  dfa$y)	
ryh<- raster(rcb); ryh <- setValues(ryh, fitted(lm.v))	
rx1<- raster(rcb); rx1 <- setValues(rx1, dfa$x1)
rxe<- raster(rcb); rxe <- setValues(rxe, residuals(lm.v))	
r.v <- brick(ry,ryh,rx1, rxe)

tits = c("y", "y*", expression(paste("x"[1])), expression(paste(epsilon)) )
for (i in 1:dim(r.v)[3]) {
	if (i == 2) p <- plot.func.fixed(r.v, var = names(r.v)[i], tit = tits[i], 
		breaks = b.ystar)
	if (i == 4) p <- plot.func.fixed(r.v, var = names(r.v)[i], tit = tits[i], 
		breaks = b.error)
	if( i==1 | i == 3 | i == 5) 
		p <- plot.func(r.v, var = names(r.v)[i], tit = tits[i]) 
	assign(paste0("p",i), p)
}
#png(filename = "F2c.png", w = 9*1.5, h = 2.3*1.2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(1,5)))
print(p1, vp=viewport(layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.col = 5))
#dev.off()

### 4. m: to show this you could have a slightly altered x1.
# get different support y*
dfa <- df
dfa$x1 <- dfa$x1 + sqrt(dfa$x1)
lm.m <- lm(rm.standard, data = dfa)
summary(lm.m)
ry <- raster(rc); ry  <- setValues(ry,  dfa$y)	
ryh<- raster(rc); ryh <- setValues(ryh, fitted(lm.m))	
rx1<- raster(rc); rx1 <- setValues(rx1, dfa$x1)	
rxe<- raster(rc); rxe <- setValues(rxe, residuals(lm.m))	
r.m <- brick(ry,ryh,rx1, rxe)
tits = c("y", "y*", expression(paste("x"[1])), expression(paste(epsilon)))
for (i in 1:dim(r.m)[3]) {
	if (i == 2) p <- plot.func.fixed(r.m, var = names(r.m)[i], tit = tits[i], 
		breaks = b.ystar)
	if (i == 4) p <- plot.func.fixed(r.m, var = names(r.m)[i], tit = tits[i], 
		breaks = b.error)
	if( i==1 | i == 3 | i == 5) 
		p <- plot.func(r.m, var = names(r.m)[i], tit = tits[i]) 
	assign(paste0("p",i), p)
}

#png(filename = "F2d.png", w = 9*1.5, h = 2.3*1.2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(1,5)))
print(p1, vp=viewport(layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.col = 5))
#dev.off()

### 5. e: would need to add some noise to x1.
dfj <- df
dfj$x1 <- jitter(dfj$x1, 500, 300)
lm.e <- lm(rm.standard, data = dfj)
summary(lm.e)
ry <- raster(rc); ry  <- setValues(ry,  dfj$y)	
ryh<- raster(rc); ryh <- setValues(ryh, fitted(lm.e))	
rx1<- raster(rc); rx1 <- setValues(rx1, dfj$x1)	
rxe<- raster(rc); rxe <- setValues(rxe, residuals(lm.e))	
r.e <- brick(ry,ryh,rx1, rxe)

tits = c("y", "y*", expression(paste("x"[1])), expression(paste(epsilon)))
for (i in 1:dim(r.e)[3]) {
	if (i == 2) p <- plot.func.fixed(r.e, var = names(r.e)[i], tit = tits[i], 
		breaks = b.ystar)
	if (i == 4) p <- plot.func.fixed(r.e, var = names(r.e)[i], tit = tits[i], 
		breaks = b.error)
	if( i==1 | i == 3 | i == 5) 
		p <- plot.func(r.e, var = names(r.e)[i], tit = tits[i]) 
	assign(paste0("p",i), p)
}
#png(filename = "F2e.png", w = 9*1.5, h = 2.3*1.2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(1,5)))
print(p1, vp=viewport(layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.col = 5))
#dev.off()

### 6. f: GWR with boxcar.
rm.standard <- as.formula(y~x1+0)
rc.sp <- as(rc, "SpatialPointsDataFrame")
coords <- coordinates(rc.sp)
dMat <- gw.dist(coords)
#bw = bw.gwr(rm.standard, data = rc.sp, approach="AIC",kernel="bisquare",
#                 adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
gm.f <- gwr.basic(rm.standard, data = rc.sp, bw = 1000, kernel="boxcar", 
				  adaptive=F, dMat = dMat)
ry <- raster(rc); ry  <- setValues(ry,  df$y)	
ryh<- raster(rc); ryh <- setValues(ryh, gm.f$SDF$yhat)	
rx1<- raster(rc); rx1 <- setValues(rx1, df$x1)	
re<- raster(rc); re <- setValues(re, gm.f$SDF$residual)
r.gm <- brick(ry,ryh,rx1,re)

tits = c("y", "y*", expression(paste("x"[1])), expression(paste(epsilon)) )
for (i in 1:dim(r.gm)[3]) {
	if (i == 2) p <- plot.func.fixed(r.gm, var = names(r.gm)[i], tit = tits[i], 
		breaks = b.ystar)
	if (i == 4) p <- plot.func.fixed(r.gm, var = names(r.gm)[i], tit = tits[i], 
		breaks = b.error)
	if( i==1 | i == 3 | i == 5) 
		p <- plot.func(r.gm, var = names(r.gm)[i], tit = tits[i]) 
	assign(paste0("p",i), p)
}
#png(filename = "F2f.png", w = 9*1.5, h = 2.3*1.2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(1,5)))
print(p1, vp=viewport(layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.col = 5))
#dev.off()

### 7. f2: to show this you would need a different fitting model. 
# Spatial linear mixed model
# see Harry email and https://stats.idre.ucla.edu/r/faq/how-do-i-model-a-spatially-autocorrelated-outcome/
df$dummy <- rep(1, nrow(df))
rm.i <- as.formula(y~x1+x2+0)
df$easting <- coords[,1]
df$northing <- coords[,2]
lm.i <- lme(fixed = rm.i, data = df, random = ~ 1 | dummy, method = "REML")
lm.f2i <- update(lm.i, correlation = corExp(form = ~easting+northing), method = "REML")
summary(lm.f2i)
ry <- raster(rc); ry  <- setValues(ry,  df$y)	
ryh<- raster(rc); ryh <- setValues(ryh, as.vector(fitted(lm.f2i)))	
rx1<- raster(rc); rx1 <- setValues(rx1, df$x1)
rx2<- raster(rc); rx2 <- setValues(rx2, df$x2)
re<- raster(rc); re <- setValues(re, as.vector(residuals(lm.f2i)))
r.i <- brick(ry,ryh,rx1,rx2,re)
tits = c("y", "y*", expression(paste("x"[1])), 
         expression(paste("x"[2])), expression(paste(epsilon)) )
for (i in 1:dim(r.m)[3]) {
	if (i == 2) p <- plot.func.fixed(r.i, var = names(r.i)[i], tit = tits[i], 
		breaks = b.ystar)
	if (i == 5) p <- plot.func.fixed(r.i, var = names(r.i)[i], tit = tits[i], 
		breaks = b.error)
	if( i==1 | i == 3 | i == 4) 
		p <- plot.func(r.i, var = names(r.m)[i], tit = tits[i]) 
	assign(paste0("p",i), p)
}
#png(filename = "F2g.png", w = 9*1.5, h = 2.3*1.2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(1,5)))
print(p1, vp=viewport(layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.col = 4))
print(p5, vp=viewport(layout.pos.col = 5))
#dev.off()

### 8. i:  
lm.ii <- update(lm.i, correlation = corExp(form = ~easting+northing), method = "ML")
summary(lm.ii)
ry <- raster(rc); ry  <- setValues(ry,  df$y)	
ryh<- raster(rc); ryh <- setValues(ryh, as.vector(fitted(lm.ii)))	
rx1<- raster(rc); rx1 <- setValues(rx1, df$x1)
rx2<- raster(rc); rx2 <- setValues(rx2, df$x2)
re<- raster(rc); re <- setValues(re, as.vector(residuals(lm.ii)))
r.i <- brick(ry,ryh,rx1,rx2,re)

tits = c("y", "y*", expression(paste("x"[1])), 
         expression(paste("x"[2])), expression(paste(epsilon)) )
for (i in 1:dim(r.m)[3]) {
	if (i == 2) p <- plot.func.fixed(r.i, var = names(r.i)[i], tit = tits[i], 
		breaks = b.ystar)
	if (i == 5) p <- plot.func.fixed(r.i, var = names(r.i)[i], tit = tits[i], 
		breaks = b.error)
	if( i==1 | i == 3 | i == 4) 
		p <- plot.func(r.i, var = names(r.m)[i], tit = tits[i]) 
	assign(paste0("p",i), p)
}
#png(filename = "F2h.png", w = 9*1.5, h = 2.3*1.2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(1,5)))
print(p1, vp=viewport(layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.col = 4))
print(p5, vp=viewport(layout.pos.col = 5))
#dev.off()

# Coef tab
# gwr
pv1<- 2*pt(- abs(gm.f$SDF@data$x1_TV), df=nrow(gm.f$SDF)-1)
gw.df <- apply(gm.f$SDF@data[,c(1,8)], 2, fivenum)
gw.df <- data.frame(gw.df, fivenum(pv1))
arse <- names(summary(lm.stand)$coefficients[, c(1,3,4)])
names(gw.df) <- arse
# reml
f2.df <- coef(summary(lm.f2i))[,c(1,4,5) ]
colnames(f2.df) <- arse
# ml
i.df <- coef(summary(lm.ii))[,c(1,4,5)]
colnames(i.df) <- arse

## Summary table of the models
tab <- rbind(summary(lm.stand)$coefficients[, c(1,3,4)],
             summary(lm.x)$coefficients[, c(1,3,4)],
             summary(lm.v)$coefficients[, c(1,3,4)],
             summary(lm.m)$coefficients[, c(1,3,4)],
             summary(lm.e)$coefficients[, c(1,3,4)],
             gw.df, 
             f2.df,
             i.df)  

tab
#write.csv(tab, "tab1.csv")

### Boxplots of y* / yhat
df3 <- data.frame(Reference = fitted(lm.stand),
                  x = (fitted(lm.x)), 
                  v = (fitted(lm.v)), 
                  m = (fitted(lm.m)), 
                  e = (fitted(lm.e)),
                  f = (gm.f$SDF@data[,c(3)]), 
                  f2 = (fitted(lm.f2i)),
                  i = (fitted(lm.ii))) 
df3 <- melt(df3)
head(df3)
#png(filename = "F3.png", w = 9, h = 2.5, units = "in", res = 600)
tit <- "y*"
ggplot(df3, aes(variable, (value))) +
  geom_boxplot(outlier.colour = "firebrick3",alpha = 0.2, outlier.size = 2) +
  ylab(tit) +
  xlab("Regression Models / Issues") +
  ylim(c(500, 2000)) + 
  theme(
#    axis.title.y =element_text(size=16, face = "bold"),
    axis.text.x=element_text(size=11, face = "bold"))
#dev.off()

### Boxplots of errors
df2 <- data.frame(Reference = residuals(lm.stand),
                  x = (residuals(lm.x)), 
                  v = (residuals(lm.v)), 
                  m = (residuals(lm.m)), 
                  e = (residuals(lm.e)),
                  f = (gm.f$SDF@data[,c(4)]),
                  f2 = (residuals(lm.f2i)),
                  i = (residuals(lm.ii)))  
df2 <- melt(df2)
head(df2)
#png(filename = "F4.png", w = 9, h = 2.5, units = "in", res = 600)
tit <- "Residuals"
ggplot(df2, aes(variable, value)) +
  geom_boxplot(outlier.colour = "firebrick3",alpha = 0.2, outlier.size = 2) +
  ylab(tit) +
  xlab("Regression Models / Issues") +
  ylim(c(750,-750))  + 
  theme(
#    axis.title.y =element_text(size=16, face = "bold"),
    axis.text.x=element_text(size=11, face = "bold"))
#dev.off()

####  do visual comparison of y*
plot.func2 <- function(r = r.stand, var = "x", tit = "y") {
	p <-tm_shape(ryh)+
		tm_raster(col = var, alpha = 1, 
			#style = "kmeans", n = 6, 
			breaks = br,
			palette = "RdBu",
     		title = tit) +
     tm_layout(
     	legend.title.size=1,
     	legend.text.size = 0.6,
     	legend.bg.color = "white",
     	legend.bg.alpha = 1) 
     return(p)
}
### comparing to y
df3 <- data.frame(Reference = fitted(lm.stand),
                  x = (fitted(lm.x)), 
                  v = (fitted(lm.v)), 
                  m = (fitted(lm.m)), 
                  e = (fitted(lm.e)),
                  f = (gm.f$SDF@data[,c(3)]), 
                  f2 = (fitted(lm.f2i)),
                  i = (fitted(lm.ii))) 

for (i in 1:8){
	df3[,i] <- df$y - df3[,i]
}
rr <- raster(rc); rr  <- setValues(rx,  df3$Reference)	
rx <- raster(rc); rx  <- setValues(rx,  df3$x)	
rv <- raster(rc); rv  <- setValues(rv,  df3$v)	
rm <- raster(rc); rm  <- setValues(rm,  df3$m)
re <- raster(rc); re  <- setValues(re,  df3$e)
rf <- raster(rc); rf  <- setValues(rf,  df3$f)	
rf2 <- raster(rc); rf2  <- setValues(rf2,  as.vector(df3$f2))	
ri <- raster(rc); ri  <- setValues(ri,  as.vector(df3$i))	
ryh <- brick(rr,rx,rv,rm,re,rf, rf2,ri)
tits = c("Reference","x", "v", "m", "e", "f", "f2", "i")
names(ryh) <- tits
summary(unlist(df3[,1:8]))
br <- c(-3500, -100, -50, 0, 50, 100, 3000)
for (i in 1:dim(ryh)[3]) {
  p <- plot.func2(ryh, var = names(ryh)[i], tit = tits[i]) 
  assign(paste0("p",i), p)
}
#png(filename = "F5.png", w = 9*1.5*4/5, h = 2.3*1.2*2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(2,4)))
print(p1, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(p5, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p6, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(p7, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
print(p8, vp=viewport(layout.pos.row = 2, layout.pos.col = 4))
#dev.off()

### comapring to x1
### comparing to y
df3 <- data.frame(Reference = fitted(lm.stand),
                  x = (fitted(lm.x)), 
                  v = (fitted(lm.v)), 
                  m = (fitted(lm.m)), 
                  e = (fitted(lm.e)),
                  f = (gm.f$SDF@data[,c(3)]), 
                  f2 = (fitted(lm.f2i)),
                  i = (fitted(lm.ii))) 

for (i in 1:8){
	df3[,i] <- df$x1 - df3[,i]
}
rr <- raster(rc); rr  <- setValues(rx,  df3$Reference)	
rx <- raster(rc); rx  <- setValues(rx,  df3$x)	
rv <- raster(rc); rv  <- setValues(rv,  df3$v)	
rm <- raster(rc); rm  <- setValues(rm,  df3$m)
re <- raster(rc); re  <- setValues(re,  df3$e)
rf <- raster(rc); rf  <- setValues(rf,  df3$f)	
rf2 <- raster(rc); rf2  <- setValues(rf2,  as.vector(df3$f2))	
ri <- raster(rc); ri  <- setValues(ri,  as.vector(df3$i))	
ryh <- brick(rr,rx,rv,rm,re,rf, rf2,ri)
tits = c("Reference","x", "v", "m", "e", "f", "f2", "i")
names(ryh) <- tits
summary(unlist(df3[,1:8]))
br <- c(-3500, -100, -50, 0, 50, 100, 3000)
for (i in 1:dim(ryh)[3]) {
  p <- plot.func2(ryh, var = names(ryh)[i], tit = tits[i]) 
  assign(paste0("p",i), p)
}
#png(filename = "F6.png", w = 9*1.5*4/5, h = 2.3*1.2*2, units = "in", res = 600)
pushViewport(viewport(layout=grid.layout(2,4)))
print(p1, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(p4, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(p5, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p6, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(p7, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
print(p8, vp=viewport(layout.pos.row = 2, layout.pos.col = 4))
#dev.off()

###### END PART 1 - Demonstration #######


###### START PART 2 - Loess Plateau case study #######
## 1. Load Data 
source_data("https://github.com/lexcomber/regression_semantics/blob/master/lp.RData?raw=True")

## convert aspect to Northness and Eastness
# https://groups.google.com/forum/#!msg/Maxent/5VA0ATKJrPI/utSe0PaV2W4J
# aspect
a <- data.30$aspect
data.30$east <- sin(a)
data.30$north <- cos(a)
a <- data.90$aspect
data.90$east <- sin(a)
data.90$north <- cos(a)
# flow direction
a <- ((data.30$flowdir/128) * 360) 
data.30$flowne <- sin(a)
a <- ((data.90$flowdir/128) * 360) 
data.90$flowne <- sin(a)

data.30$Ks <-data.30$Ks.0.10cm.
data.90$Ks <-data.90$Ks.0.10cm.

# see http://spatialreference.org/ref/?search=krasovsky+1940+Albers
tmp.proj <- CRS("+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +x_0=4000000 +y_0=0 +datum=WGS84 +units=m +no_defs")
crs.val <- CRS("+proj=longlat +datum=WGS84 ")
out.ll <- spTransform(out, crs.val)
data.ll <- spTransform(data.30, crs.val)
# Study Area Figure
ul <- as.vector(cbind(bbox(out.ll)[2,2], 
                      bbox(out.ll)[1,1]))
lr <- as.vector(cbind(bbox(out.ll)[2,1], 
                      bbox(out.ll)[1,2]))
MyMap <- openmap(ul,lr, type = "stamen-terrain")
### Figure 7
#png("F7.png", width = 10, height = 6, units = 'in', res = 300)
par(mar = c(0,0,0,0))
plot(MyMap, removeMargin=FALSE)
plot(spTransform(out.ll, osm()), add = TRUE, lwd = 2)
plot(spTransform(data.ll, osm()),bg = "#25252580", 
     add = TRUE, pch = 21, cex = 1.5)
scalebar(200000, type = "bar", xy = c(11329039, 4081403), 
         divs = 4, below = "km",
         label = c(0,100,200) )
#dev.off()

## 2. Data Prep, model chopice and Descriptive statistics and plots
# Model choice
tmp =  data.30@data[, c(4,7,10,16,19,22,31:40)]
summary(stepAIC(lm(Ks~.+0, data = tmp), trace = F))
summary(stepAIC(lm(Ks~., data = tmp), trace = F))
# Clay, Silt, SSWC, east flowne
# Data Prep
df <- data.30@data[, c(40,4:9,22:24,37,39)]
df <- data.frame(df, data.90@data[, c(37,39)])
names(df)[11:14] <- c("east30", "flowne30","east90", "flowne90")
names(df)[2:4] <- c("clay10", "clay20","clay40")
names(df)[5:7] <- c("silt10", "silt20","silt40")
names(df)[8:10] <- c("sswc10", "sswc20","sswc40")
# get rid of 20cm
df <- df[, c(1,2,4,5,7,8,10:14)]
data <- SpatialPointsDataFrame(data.30, data = data.frame(df),proj4string = tmp.proj)
# Descriptive plots: soil physical
tmp <- data@data[,2:7]
depths <- c("0 to 10", "20 to 40")
m<- melt(tmp)
head(m)
depth.i <- rep(NA, nrow(m))
i <- grep("10", m$variable)
depth.i[i] <- depths[1]
i <- grep("40", m$variable)
depth.i[i] <- depths[2]
m <- data.frame(m, depth = depth.i)
head(m)
m$variable <- gsub("10", "",m$variable )
m$variable <- gsub("40", "",m$variable )
m$depth <- factor(m$depth, levels = c("20 to 40","0 to 10" ))

# check 
table(m$variable)
#png(filename = "F9.png", w = 9, h = 2.5, units = "in", res = 300)
ggplot(m, aes(
      y = as.numeric(value),
      x = as.factor(depth),
      fill = as.factor(depth))
  )+ 
  geom_boxplot(show.legend=FALSE, outlier.colour = "firebrick3",
  	outlier.alpha = 0.35, outlier.size = 2) +
  #scale_y_continuous(limits = quantile(m$value, c(0.1, 0.9)))+
  coord_flip() +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  scale_fill_brewer(type="qual", palette = "Reds") +
  ylab("Variable Value") +
  xlab("Measurement Depth") +
  theme_bw() 
¢dev.off()
# Descriptive plots: terrain
tmp <- data@data[,8:11]
res <- c("30m", "90m")
m<- melt(tmp)
head(m)
res.i <- rep(NA, nrow(m))
i <- grep("30", m$variable)
res.i[i] <- res[1]
i <- grep("90", m$variable)
res.i[i] <- res[2]
m <- data.frame(m, res = res.i)
m$variable <- gsub("30", "",m$variable )
m$variable <- gsub("90", "",m$variable )
head(m)
table(m$variable)
#png(filename = "F8a.png", w = 9/2, h = 2.5, units = "in", res = 300)
ggplot(m[m$variable == "east",], aes(x = value)) +
  geom_histogram(aes(y=..density..),
                   binwidth=0.2,colour="white") +
    geom_density(alpha=.2, fill="#FF6666") +
    facet_wrap( ~ res, ncol = 4) +
    theme(axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab("Variable") +
  	xlab("Flow-NE") 
#dev.off()  

#png(filename = "F8b.png", w = 9/2, h = 2.5, units = "in", res = 300)
ggplot(m[m$variable != "east",], aes(x = value)) +
  geom_histogram(aes(y=..density..),
                   binwidth=0.2,colour="white") +
    geom_density(alpha=.2, fill="#FF6666") +
    facet_wrap( ~ res, ncol = 4) +
    theme(axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab("Variable") +
  	xlab("East-ness") 
#dev.off()  
 
## 3. Regressions
# 3.1 Reference: y, y*, x1 where y predicted using a simple linear regression model only.
rm.0 <- as.formula(Ks~sswc10+east30+flowne30+0)
lm.0 <- lm(rm.0, data = data)
summary(lm.0)
# 3.1 The x issue: the number and choice of potential covariates
rm.x <- as.formula(Ks~clay10+silt10+sswc10+east30+flowne30+0)
lm.x <- lm(rm.x, data = data)
summary(lm.x)
# 3.3 The v issue: observation support: just aggregate the x1 data 
rm.v <- as.formula(Ks~sswc10+east90+flowne90+0)
lm.v <- lm(rm.v, data = data)
summary(lm.v)
# 3.4 The m issue: the method of measurement
rm.m <- as.formula(Ks~sswc40+east30+flowne30+0)
#rm.m <- as.formula(Ks~sswc40+east30+flowne30+0)
lm.m <- lm(rm.m, data = data)
summary(lm.m)
# 3.5 The first f  issue: model specification need to select a different model. – GWR coefficients
dMat2 = gw.dist(coordinates(data)) 
bw = bw.gwr(rm.0, adaptive = F, kernel = "gaussian", dMat = dMat2, data = data)
gw.f2 <- gwr.basic(rm.0, bw = bw, adaptive = T, kernel = "gaussian", data = data, dMat = dMat2)
summary(gw.f2$SDF@data[, c(1:3,5)])
# 3.6 The second f issue
# Spatial linear mixed model
data$dummy <- rep(1, nrow(data))
rm.x <- as.formula(Ks~clay10+silt10+sswc10+east30+flowne30+0)
data$easting <- coordinates(data)[,1]
data$northing <- coordinates(data)[,2]
lm.f2 <- lme(fixed = rm.x, data = data@data, random = ~ 1 | dummy, method = "REML")
summary(lm.f2)
lm.f2i <- update(lm.f2, correlation = corExp(form = ~easting+northing), method = "REML")
summary(lm.f2i)
coef(summary(lm.f2i))[,c(1,2,4,5)]
# 3.7 The i issue: model identification -  a mixed-effects regression 
lm.i <- update(lm.f2, correlation = corExp(form = ~easting+northing), method = "ML")
summary(lm.i)
coef(summary(lm.i))[,c(1,2,4,5)]

## P values for gwr coeffs
pv1<- 2*pt(- abs(gw.f2$SDF@data$sswc10_TV), df=nrow(gw.f2$SDF)-1)
pv2<- 2*pt(- abs(gw.f2$SDF@data$east30_TV), df=nrow(gw.f2$SDF)-1)
pv3<- 2*pt(- abs(gw.f2$SDF@data$flowne30_TV), df=nrow(gw.f2$SDF)-1)
pv.df <- data.frame(pv1, pv2, pv3)
pv.iqr <- apply(pv.df, 2, IQR)
tv.iqr <- apply(gw.f2$SDF@data[, c(12:14)], 2, IQR)
coef.iqr <- apply(gw.f2$SDF@data[, c(1:3)], 2, IQR)
gw.df <- data.frame(coef.iqr, tv.iqr, pv.iqr)
arse <- colnames(summary(lm.0)$coefficients[, c(1,3,4)])
names(gw.df) <- arse
# other names
f2.df <- coef(summary(lm.f2i))[,c(1,4,5)]
i.df <- coef(summary(lm.i))[,c(1,4,5)]
colnames(f2.df) <- arse
colnames(i.df) <- arse

## Summary table of the models
tab <- rbind(summary(lm.0)$coefficients[c(1,3,2), c(1,3,4)],
             summary(lm.x)$coefficients[c(3,5,4,1,2), c(1,3,4)],
             summary(lm.v)$coefficients[c(1,3,2), c(1,3,4)],
             summary(lm.m)$coefficients[c(1,3,2), c(1,3,4)],
             gw.df[c(1,3,2),],
             f2.df[c(3,5,4,1,2),],
             i.df[c(3,5,4,1,2),])
#mod.list <- c("lm.stand", "lm.x", "lm.x", "lm.v", "lm.m", "lm.e", "gm.f", "gm.f","gm.f","gm.f",
#              "gm.f","lm.f2", "lm.f2","lm.i", "lm.i")
#rownames(tab) <- mod.list
#write.csv(tab, "tab2.csv")

df <- data.frame(y.0 = fitted(lm.0),
				 y.x = fitted(lm.x),
				 y.v = fitted(lm.v),
				 y.m = fitted(lm.m),
				 y.f = gw.f2$SDF@data$yhat,
				 y.f2 = as.vector(fitted(lm.f2i)),
				 y.i = as.vector(fitted(lm.i)))
df.sp <- SpatialPointsDataFrame(data, df, proj4string = crs.val)
#png(filename = "F10.png", w = 7, h = 10, units = "in", res = 300)
tm_shape(out) +
     tm_polygons(colour="grey") +
     tm_shape(df.sp) +
     	tm_dots(col = c("y.0", "y.x", "y.v","y.m","y.f","y.f2", "y.i"), size = 0.6, 
     		alpha = 1, palette = "-RdBu",
     		#style = "equal", n = 5,
     		breaks = c(-0.15,0.3,0.5,1,1.7),midpoint = NA,
     		#title = paste("Model", c(1,2,3,4,6,7))) +		
     		title = append("Reference", paste0(c("x ","v ", "m ", "f ","f2 ", "i "), "issue"))) +
     tm_facets(ncol=2) +
     tm_layout(legend.bg.color = "grey95")
#dev.off()

###### END PART 2 - Loess Plateau case study #######
