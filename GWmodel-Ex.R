pkgname <- "GWmodel"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('GWmodel')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("DubVoter")
### * DubVoter

flush(stderr()); flush(stdout())

### Name: DubVoter
### Title: Voter turnout data in Greater Dublin(SpatialPolygonsDataFrame)
### Aliases: DubVoter Dub.voter
### Keywords: data

### ** Examples


data(DubVoter)
ls()
## Not run: 
##D spplot(Dub.voter,names(Dub.voter)[4:12])
## End(Not run)



cleanEx()
nameEx("EN_CB")
### * EN_CB

flush(stderr()); flush(stdout())

### Name: EN_CB
### Title: England county boundaries
### Aliases: cty_eng
### Keywords: data

### ** Examples

data(EN_CB)
ls()



cleanEx()
nameEx("EWHP")
### * EWHP

flush(stderr()); flush(stdout())

### Name: EWHP
### Title: House price data set (DataFrame) in England and Wales
### Aliases: EWHP ewhp
### Keywords: data

### ** Examples

###
data(EWHP)
head(ewhp)
houses.spdf <- SpatialPointsDataFrame(ewhp[, 1:2], ewhp)
 ####Get the border of England and Wales
data(EWOutline)
plot(ewoutline)
plot(houses.spdf, add = TRUE, pch = 16)



cleanEx()
nameEx("GE2015")
### * GE2015

flush(stderr()); flush(stdout())

### Name: GE2015
### Title: General Election outcome data at constituency level in England
### Aliases: ge2015
### Keywords: data

### ** Examples

data(GE2015)
ls()



cleanEx()
nameEx("Georgia")
### * Georgia

flush(stderr()); flush(stdout())

### Name: Georgia
### Title: Georgia census data set (csv file)
### Aliases: Georgia Gedu.df
### Keywords: data

### ** Examples

data(Georgia)
ls()
coords <- cbind(Gedu.df$X, Gedu.df$Y)
educ.spdf <- SpatialPointsDataFrame(coords, Gedu.df)
spplot(educ.spdf, names(educ.spdf)[4:10])



cleanEx()
nameEx("GeorgiaCounties")
### * GeorgiaCounties

flush(stderr()); flush(stdout())

### Name: GeorgiaCounties
### Title: Georgia counties data (SpatialPolygonsDataFrame)
### Aliases: Gedu.counties
### Keywords: data

### ** Examples

data(GeorgiaCounties)
plot(Gedu.counties)
data(Georgia)
coords <- cbind(Gedu.df$X, Gedu.df$Y)
educ.spdf <- SpatialPointsDataFrame(coords, Gedu.df)
plot(educ.spdf, add=TRUE)




cleanEx()
nameEx("LondonHP")
### * LondonHP

flush(stderr()); flush(stdout())

### Name: LondonHP
### Title: London house price data set (SpatialPointsDataFrame)
### Aliases: LondonHP londonhp
### Keywords: data

### ** Examples

data(LondonHP)
data(LondonBorough)
ls()
plot(londonborough)
plot(londonhp, add=TRUE)



cleanEx()
nameEx("USelect")
### * USelect

flush(stderr()); flush(stdout())

### Name: USelect
### Title: Results of the 2004 US presidential election at the county level
###   (SpatialPolygonsDataFrame)
### Aliases: USelect2004
### Keywords: data

### ** Examples


data(USelect)
ls()



cleanEx()
nameEx("ggwr.basic")
### * ggwr.basic

flush(stderr()); flush(stdout())

### Name: ggwr.basic
### Title: Generalised GWR models with Poisson and Binomial options
### Aliases: ggwr.basic gwr.generalised gwr.binomial gwr.binomial.wt
###   gwr.poisson gwr.poisson.wt gwr.fitted print.ggwrm
### Keywords: generalised GWR

### ** Examples

data(LondonHP)
## Not run: 
##D DM<-gw.dist(dp.locat=coordinates(londonhp))
##D bw.f1 <- bw.ggwr(BATH2~FLOORSZ,data=londonhp, dMat=DM)
##D res.poisson<-ggwr.basic(BATH2~FLOORSZ, bw=bw.f1,data=londonhp, dMat=DM)
##D bw.f2 <- bw.ggwr(BATH2~FLOORSZ,data=londonhp, dMat=DM,family ="binomial")
##D res.binomial<-ggwr.basic(BATH2~FLOORSZ, bw=bw.f2,data=londonhp, dMat=DM,
##D               family ="binomial")
## End(Not run)



cleanEx()
nameEx("gw.dist")
### * gw.dist

flush(stderr()); flush(stdout())

### Name: gw.dist
### Title: Distance matrix calculation
### Aliases: gw.dist coordinate_rotate eu_dist_mat eu_dist_smat eu_dist_vec
###   mk_dist_mat mk_dist_smat mk_dist_vec cd_dist_mat cd_dist_smat
###   cd_dist_vec md_dist_mat md_dist_smat md_dist_vec
### Keywords: GW tools

### ** Examples

dp<-cbind(sample(100),sample(100))
rp<-cbind(sample(10),sample(10))
#Euclidean distance metric is used.
dist.v1<-gw.dist(dp.locat=dp, focus=5, p=2, theta=0, longlat=FALSE)
#Manhattan distance metric is used.
#The coordinate system is rotated by an angle 0.5 in radian.
dist.v2<-gw.dist(dp.locat=dp, focus=5, p=1, theta=0.5)
#Great Circle distance metric is used.
dist.v3<-gw.dist(dp.locat=dp, focus=5, longlat=TRUE)
#A generalized Minkowski distance metric is used with p= 0.75 .
#The coordinate system is rotated by an angle 0.8 in radian.
dist.v4<-gw.dist(dp.locat=dp,rp.locat=rp, focus=5, p=0.75,theta=0.8)
################################
#matrix is calculated
#Euclidean distance metric is used.
dist.m1<-gw.dist(dp.locat=dp, p=2, theta=0, longlat=FALSE)
#Manhattan distance metric is used.
#The coordinate system is rotated by an angle 0.5 in radian.
dist.m2<-gw.dist(dp.locat=dp, p=1, theta=0.5)
#Great Circle distance metric is used.
#dist.m3<-gw.dist(dp.locat=dp, longlat=TRUE)
#A generalized Minkowski distance metric is used with p= 0.75 .
#The coordinate system is rotated by an angle 0.8 in radian.
dist.m4<-gw.dist(dp.locat=dp,rp.locat=rp, p=0.75,theta=0.8)



cleanEx()
nameEx("gwda")
### * gwda

flush(stderr()); flush(stdout())

### Name: gwda
### Title: GW Discriminant Analysis
### Aliases: gwda print.gwda grouping.xy wqda wlda splitx wmean wvarcov
###   wprior confusion.matrix
### Keywords: GWDA

### ** Examples

## Not run: 
##D  require(tmap)
##D  data(ge2015)
##D  data(cty_eng)
##D  ge2015 <- ge2015[ge2015$WINNER ##D 
##D  dMat <- gw.dist(coordinates(ge2015))
##D  bw <- bw.gwda(WINNER~Age65over+OwnOcc+NoQual+Unemp+NonWhite+LoneParHH,data=ge2015,
##D  adaptive=TRUE,dMat=dMat)
##D  ge.gwda <- gwda(WINNER~Age65over+OwnOcc+NoQual+Unemp+NonWhite+LoneParHH,data=ge2015,
##D  bw=bw,adaptive=TRUE,dMat=dMat)
##D  table(ge2015$WINNER,ge.gwda$SDF$group.predicted)
##D  tm_shape(ge.gwda$SDF)+tm_fill("entropy")+tm_shape(cty_eng)+tm_borders()
##D  
## End(Not run)



cleanEx()
nameEx("gwpca")
### * gwpca

flush(stderr()); flush(stdout())

### Name: gwpca
### Title: GWPCA
### Aliases: gwpca robustSvd rwpca wpca wt.median print.gwpca
### Keywords: GWPCA

### ** Examples

## Not run: 
##D if(require("mvoutlier") && require("RColorBrewer"))
##D {
##D   data(bsstop)
##D   Data.1 <- bsstop[, 1:14]
##D   colnames(Data.1)
##D   Data.1.scaled <- scale(as.matrix(Data.1[5:14]))  # standardised data...
##D   rownames(Data.1.scaled) <- Data.1[, 1]
##D   #compute principal components:
##D   pca <- princomp(Data.1.scaled, cor = FALSE, scores = TRUE)  
##D   # use covariance matrix to match the following...
##D   pca$loadings
##D   data(bss.background)
##D   backdrop <- function() 
##D    plot(bss.background, asp = 1, type = "l", xaxt = "n", yaxt = "n", 
##D    xlab = "", ylab = "", bty = "n", col = "grey")
##D   pc1 <- pca$scores[, 1]
##D   backdrop()
##D   points(Data.1$XCOO[pc1 > 0], Data.1$YCOO[pc1 > 0], pch = 16, col = "blue")
##D   points(Data.1$XCOO[pc1 < 0], Data.1$YCOO[pc1 < 0], pch = 16, col = "red")
##D   
##D   #Geographically Weighted PCA and mapping the local loadings
##D   # Coordinates of the sites
##D   Coords1 <- as.matrix(cbind(Data.1$XCOO,Data.1$YCOO)) 
##D   d1s <- SpatialPointsDataFrame(Coords1,as.data.frame(Data.1.scaled))
##D   pca.gw <- gwpca(d1s,vars=colnames(d1s@data),bw=1000000,k=10)
##D   local.loadings <- pca.gw$loadings[, , 1]  
##D   
##D   # Mapping the winning variable with the highest absolute loading
##D   # note first component only - would need to explore all components..
##D   
##D   lead.item <- colnames(local.loadings)[max.col(abs(local.loadings))]
##D   df1p = SpatialPointsDataFrame(Coords1, data.frame(lead = lead.item))
##D   backdrop()
##D   colour <- brewer.pal(8, "Dark2")[match(df1p$lead, unique(df1p$lead))]
##D   plot(df1p, pch = 18, col = colour, add = TRUE)
##D   legend("topleft", as.character(unique(df1p$lead)), pch = 18, col = 
##D       brewer.pal(8, "Dark2"))
##D   backdrop()
##D   
##D   #Glyph plots give a view of all the local loadings together
##D   glyph.plot(local.loadings, Coords1, add = TRUE)
##D   
##D   #it is not immediately clear how to interpret the glyphs fully, 
##D   #so inter-actively identify the full loading information using:
##D   check.components(local.loadings, Coords1)
##D   
##D   # GWPCA with an optimal bandwidth
##D   bw.choice <- bw.gwpca(d1s,vars=colnames(d1s@data),k=2) 
##D   pca.gw.auto  <- gwpca(d1s,vars=colnames(d1s@data),bw=bw.choice,k=2)
##D   # note first component only - would need to explore all components..
##D   local.loadings <- pca.gw.auto$loadings[, , 1]  
##D   
##D   lead.item <- colnames(local.loadings)[max.col(abs(local.loadings))]
##D   df1p = SpatialPointsDataFrame(Coords1, data.frame(lead = lead.item))
##D   backdrop()
##D   colour <- brewer.pal(8, "Dark2")[match(df1p$lead, unique(df1p$lead))]
##D   plot(df1p, pch = 18, col = colour, add = TRUE)
##D   legend("topleft", as.character(unique(df1p$lead)), pch = 18, 
##D   col = brewer.pal(8, "Dark2"))
##D   
##D   # GWPCPLOT for investigating the raw multivariate data
##D   gw.pcplot(d1s, vars=colnames(d1s@data),focus=359, bw = bw.choice) 
##D }
## End(Not run)



cleanEx()
nameEx("gwpca.montecarlo.1")
### * gwpca.montecarlo.1

flush(stderr()); flush(stdout())

### Name: gwpca.montecarlo.1
### Title: Monte Carlo (randomisation) test for significance of GWPCA
###   eigenvalue variability for the first component only - option 1
### Aliases: gwpca.montecarlo.1 montecarlo.gwpca.1 plot.mcsims
### Keywords: GWPCA

### ** Examples

## Not run: 
##D data(DubVoter)
##D DM<-gw.dist(dp.locat=coordinates(Dub.voter))
##D gmc.res<-gwpca.montecarlo.1(data=Dub.voter, vars=c("DiffAdd", "LARent",
##D "SC1", "Unempl", "LowEduc"), bw=20,dMat=DM,adaptive=TRUE)
##D gmc.res
##D plot(gmc.res)
## End(Not run)



cleanEx()
nameEx("gwpca.montecarlo.2")
### * gwpca.montecarlo.2

flush(stderr()); flush(stdout())

### Name: gwpca.montecarlo.2
### Title: Monte Carlo (randomisation) test for significance of GWPCA
###   eigenvalue variability for the first component only - option 2
### Aliases: gwpca.montecarlo.2 montecarlo.gwpca.2
### Keywords: GWPCA

### ** Examples

## Not run: 
##D data(DubVoter)
##D DM<-gw.dist(dp.locat=coordinates(Dub.voter))
##D gmc.res.autow<-gwpca.montecarlo.2(data=Dub.voter, vars=c("DiffAdd", "LARent",
##D "SC1", "Unempl", "LowEduc"), dMat=DM,adaptive=TRUE)
##D gmc.res.autow
##D plot.mcsims(gmc.res.autow)
## End(Not run)



cleanEx()
nameEx("gwr.basic")
### * gwr.basic

flush(stderr()); flush(stdout())

### Name: gwr.basic
### Title: Basic GWR model
### Aliases: gwr.basic gw_reg gwr_diag Ci_mat F1234.test print.gwrm
### Keywords: GWR

### ** Examples

data(LondonHP)
DM<-gw.dist(dp.locat=coordinates(londonhp))
##Compare the time consumed with and without a specified distance matrix
## Not run: 
##D system.time(gwr.res<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=1000,
##D             kernel = "gaussian"))
##D system.time(DM<-gw.dist(dp.locat=coordinates(londonhp)))
##D system.time(gwr.res<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=1000,
##D             kernel = "gaussian", dMat=DM))
##D 
##D ## specify an optimum bandwidth by cross-validation appraoch
##D bw1<-bw.gwr(PURCHASE~FLOORSZ, data=londonhp, kernel = "gaussian",dMat=DM)
##D gwr.res1<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=bw1,kernel = "gaussian", 
##D                    dMat=DM)
##D gwr.res1 
## End(Not run)
data(LondonBorough)

nsa = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(561900,200900), 
scale = 500, col=1)
## Not run: 
##D if(require("RColorBrewer"))
##D {
##D   mypalette<-brewer.pal(6,"Spectral")
##D   x11()
##D   spplot(gwr.res1$SDF, "FLOORSZ", key.space = "right", cex=1.5, cuts=10,
##D   ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
##D   main="GWR estimated coefficients for FLOORSZ with a fixed bandwidth", 
##D   col.regions=mypalette, sp.layout=list(nsa, londonborough))}
## End(Not run)
## Not run: 
##D bw2<-bw.gwr(PURCHASE~FLOORSZ,approach="aic",adaptive=TRUE, data=londonhp, 
##D             kernel = "gaussian", dMat=DM)
##D gwr.res2<-gwr.basic(PURCHASE~FLOORSZ, data=londonhp, bw=bw2,adaptive=TRUE,
##D                     kernel = "gaussian", dMat=DM)
##D gwr.res2
##D if(require("RColorBrewer"))
##D {
##D   x11()
##D   spplot(gwr.res2$SDF, "FLOORSZ", key.space = "right", cex=1.5, cuts=10,
##D   ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
##D   main="GWR estimated coefficients for FLOORSZ with an adaptive bandwidth", 
##D   col.regions=mypalette, sp.layout=list(nsa,londonborough))}
## End(Not run)



cleanEx()
nameEx("gwr.bootstrap")
### * gwr.bootstrap

flush(stderr()); flush(stdout())

### Name: gwr.bootstrap
### Title: Bootstrap GWR
### Aliases: gwr.bootstrap print.gwrbsm generate.lm.data parametric.bs
###   parametric.bs.local se.bs bias.bs ci.bs pval.bs gwrtvar gwrt.mlr
###   gwrt.lag gwrt.err gwrt.sma bw.gwr3
### Keywords: GWR

### ** Examples

## Not run: 
##D #Example with the Georgia educational attainment data
##D data(Georgia)
##D data(GeorgiaCounties)
##D coords <- cbind(Gedu.df$X, Gedu.df$Y)
##D Gedu.spdf <- SpatialPointsDataFrame(coords, Gedu.df)
##D #Make a SpatialPolygonDataFrame
##D require(RColorBrewer)
##D gSRDF <- SpatialPolygonsDataFrame(polygons(Gedu.counties), over(Gedu.counties, 
##D                                   Gedu.spdf),match.ID=T)  
##D mypalette.1 <- brewer.pal(11,"Spectral")
##D X11(width=9,height=8)                   
##D spplot(gSRDF, names(gSRDF)[c(5,7:9)], col.regions=mypalette.1,
##D cuts=10, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Georgia educational attainment predictor data")))
##D bsm.res <- gwr.bootstrap(PctBach~PctRural+PctEld+PctFB+PctPov, gSRDF, 
##D                          R=999, longlat=T)
##D bsm.res
##D #local bootstrap tests with respect to: MLR, ERR, SMA and LAG models.
##D mypalette.local.test <- brewer.pal(10,"Spectral")
##D X11(width=12,height=16)
##D spplot(bsm.res$SDF, names(bsm.res$SDF)[14:17], col.regions=mypalette.local.test,
##D cuts=9, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Local p-values for each coefficient of the MLR model 
##D                        null hypothesis")))
##D 
##D X11(width=12,height=16)
##D spplot(bsm.res$SDF, names(bsm.res$SDF)[19:22], col.regions=mypalette.local.test,
##D cuts=9, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Local p-values for each coefficient of the ERR model 
##D                        null hypothesis")))
##D X11(width=12,height=16)
##D spplot(bsm.res$SDF, names(bsm.res$SDF)[24:27], col.regions=mypalette.local.test,
##D cuts=9, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Local p-values for each coefficient of the SMA model null
##D                        hypothesis")))
##D 
##D X11(width=12,height=16)
##D spplot(bsm.res$SDF, names(bsm.res$SDF)[29:32], col.regions=mypalette.local.test,
##D cuts=9, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Local p-values for each coefficient of the LAG model null
##D                        hypothesis")))
##D ################################################################################
##D #Example with Dublin voter data
##D data(DubVoter)
##D X11(width=9,height=8)                   
##D spplot(Dub.voter, names(Dub.voter)[c(5,7,9,10)], col.regions=mypalette.1,
##D cuts=10, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Dublin voter turnout predictor data")))
##D bsm.res1 <- gwr.bootstrap(GenEl2004~LARent+Unempl+Age18_24+Age25_44, Dub.voter
##D                          , R=999)
##D bsm.res1
##D 
##D #local bootstrap tests with respect to: MLR, ERR, SMA and LAG models.
##D X11(width=11,height=8)
##D spplot(bsm.res1$SDF, names(bsm.res1$SDF)[14:17], col.regions=mypalette.local.test,
##D cuts=9, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Local p-values for each coefficient of the MLR model null
##D                         hypothesis")))
##D X11(width=11,height=8)
##D spplot(bsm.res1$SDF, names(bsm.res1$SDF)[19:22], col.regions=mypalette.local.test,
##D cuts=9, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Local p-values for each coefficient of the ERR model null
##D                         hypothesis")))
##D X11(width=11,height=8)
##D spplot(bsm.res1$SDF, names(bsm.res1$SDF)[24:27], col.regions=mypalette.local.test,
##D cuts=9, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Local p-values for each coefficient of the SMA model 
##D                             null hypothesis")))
##D X11(width=11,height=8)
##D spplot(bsm.res1$SDF, names(bsm.res1$SDF)[29:32], col.regions=mypalette.local.test,
##D cuts=9, par.settings=list(fontsize=list(text=15)),
##D main=expression(paste("Local p-values for each coefficient of the LAG model 
##D                             null hypothesis")))
## End(Not run)



cleanEx()
nameEx("gwr.lcr")
### * gwr.lcr

flush(stderr()); flush(stdout())

### Name: gwr.lcr
### Title: GWR with a locally-compensated ridge term
### Aliases: gwr.lcr ridge.lm print.gwrlcr
### Keywords: GWR-LCR

### ** Examples

data(DubVoter)
require(RColorBrewer)

# Function to find the global condition number (CN)
BKW_cn <- function (X) {
  p <- dim(X)[2]
  Xscale <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  Xsvd <- svd(Xscale)$d
  cn <- Xsvd[1] / Xsvd[p]
  cn
}
#
X <- cbind(1,Dub.voter@data[,3:10])
head(X)
CN.global <- BKW_cn(X)
CN.global
## Not run: 
##D # gwr.lcr function with a global bandwidth to check that the global CN is found
##D gwr.lcr1 <- gwr.lcr(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
##D +Age25_44+Age45_64, data=Dub.voter, bw=10000000000)
##D summary(gwr.lcr1$SDF$Local_CN)
##D 
##D # Find and map the local CNs from a basic GWR fit using the lcr-gwr function 
##D #(note this is NOT the locally-compensated ridge GWR fit as would need to set 
##D #lambda.adjust=TRUE and cn.thresh=30, say)
##D 
##D bw.lcr2 <- bw.gwr.lcr(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
##D +Age25_44+Age45_64, data=Dub.voter, kernel="bisquare", adaptive=TRUE)
##D gwr.lcr2 <- gwr.lcr(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
##D +Age25_44+Age45_64, data=Dub.voter, bw=bw.lcr2, kernel="bisquare", adaptive=TRUE)
##D if(require("RColorBrewer"))
##D   spplot(gwr.lcr2$SDF,"Local_CN",col.regions=brewer.pal(9,"YlOrRd"),cuts=8,
##D   main="Local CN")
## End(Not run)



cleanEx()
nameEx("gwr.model.view")
### * gwr.model.view

flush(stderr()); flush(stdout())

### Name: gwr.model.view
### Title: Visualise the GWR models from 'gwr.model.selection'
### Aliases: gwr.model.view model.view.gwr
### Keywords: GWR

### ** Examples

## Not run: 
##D data(LondonHP)
##D DM<-gw.dist(dp.locat=coordinates(londonhp))
##D DeVar<-"PURCHASE"
##D InDeVars<-c("FLOORSZ","GARAGE1","BLDPWW1","BLDPOSTW")
##D model.sel<-gwr.model.selection(DeVar,InDeVars, data=londonhp,
##D kernel = "gaussian", dMat=DM,bw=5000)
##D model.list<-model.sel[[1]]
##D gwr.model.view(DeVar, InDeVars, model.list=model.list)
## End(Not run)



cleanEx()
nameEx("gwr.montecarlo")
### * gwr.montecarlo

flush(stderr()); flush(stdout())

### Name: gwr.montecarlo
### Title: Monte Carlo (randomisation) test for significance of GWR
###   parameter variability
### Aliases: gwr.montecarlo montecarlo.gwr
### Keywords: GWR

### ** Examples

## Not run: 
##D data(LondonHP)
##D DM<-gw.dist(dp.locat=coordinates(londonhp))
##D bw<-bw.gwr(PURCHASE~FLOORSZ,data=londonhp,dMat=DM, kernel="gaussian")
##D #See any difference in the next two commands and why?
##D res.mont1<-gwr.montecarlo(PURCHASE~PROF+FLOORSZ, data = londonhp,dMat=DM,
##D nsim=99, kernel="gaussian", adaptive=FALSE, bw=3000)
##D res.mont2<-gwr.montecarlo(PURCHASE~PROF+FLOORSZ, data = londonhp,dMat=DM,
##D nsim=99, kernel="gaussian", adaptive=FALSE, bw=300000000000)
## End(Not run)



cleanEx()
nameEx("gwr.multiscale")
### * gwr.multiscale

flush(stderr()); flush(stdout())

### Name: gwr.multiscale
### Title: Multiscale GWR
### Aliases: gwr.multiscale gwr.q2 print.multiscalegwr gwr.backfit
### Keywords: multiscale GWR

### ** Examples

data(LondonHP)
EUDM <- gw.dist(coordinates(londonhp))
#No bandwidth is selected, and bws0 values are used
## Not run: 
##D ###Similar as the basic GWR
##D res1<-gwr.multiscale(PURCHASE~FLOORSZ+PROF, data=londonhp, criterion="dCVR",kernel="gaussian", 
##D adaptive=T, bws0=c(100, 100, 100),bw.seled=rep(T, 3), dMats=list(EUDM,EUDM,EUDM))
##D #FBGWR
##D res2<-gwr.multiscale(PURCHASE~FLOORSZ+PROF, data=londonhp, criterion="dCVR",kernel="gaussian",
##D adaptive=T, bws0=c(100, 100, 100), dMats=list(EUDM,EUDM,EUDM))
##D #Mixed GWR
##D res3<-gwr.multiscale(PURCHASE~FLOORSZ+PROF, data=londonhp, bws0=c(Inf, 100, 100, Inf),
##D                bw.seled=rep(T, 3),kernel="gaussian", dMats=list(EUDM,EUDM,EUDM))
##D #PSDM GWR
##D res4<- gwr.multiscale(PURCHASE~FLOORSZ+PROF, data=londonhp, kernel="gaussian", p.vals=c(1,2,3))
## End(Not run)



cleanEx()
nameEx("gwr.predict")
### * gwr.predict

flush(stderr()); flush(stdout())

### Name: gwr.predict
### Title: GWR used as a spatial predictor
### Aliases: gwr.predict gw.reg1 print.gwrm.pred
### Keywords: GWR

### ** Examples

## Not run: 
##D data(LondonHP)
##D gwr.pred<-gwr.predict(PURCHASE~FLOORSZ, data=londonhp, bw=2000,kernel = "gaussian")
##D gwr.pred
##D #########Global OLS regression results and comparison with gstat functions
##D if(require("gstat"))
##D {
##D   mlr.g <- gstat(id = "xx1", formula = PURCHASE~FLOORSZ,data=londonhp)
##D   mlr.g1 <- predict(mlr.g, newdata = londonhp, BLUE = TRUE)
##D   mlr.g1
##D }
##D ############
##D ols.pred<-gwr.predict(PURCHASE~FLOORSZ, data=londonhp, bw=100000000000000000000000)
##D ols.pred$SDF
## End(Not run)



cleanEx()
nameEx("gwr.robust")
### * gwr.robust

flush(stderr()); flush(stdout())

### Name: gwr.robust
### Title: Robust GWR model
### Aliases: gwr.robust
### Keywords: robust GWR

### ** Examples

## Not run: 
##D data(DubVoter)
##D bw.a <- bw.gwr(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
##D +Age25_44+Age45_64,
##D data=Dub.voter,approach="AICc",kernel="bisquare",adaptive=TRUE)
##D bw.a
##D gwr.res <- gwr.basic(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
##D +Age25_44+Age45_64,
##D data=Dub.voter,bw=bw.a,kernel="bisquare",adaptive=TRUE,F123.test=TRUE)
##D print(gwr.res)
##D 
##D # Map of the estimated coefficients for LowEduc
##D names(gwr.res$SDF)
##D if(require("RColorBrewer"))
##D {
##D   mypalette<-brewer.pal(6,"Spectral")
##D   X11(width=10,height=12)
##D   spplot(gwr.res$SDF,"LowEduc",key.space = "right",
##D   col.regions=mypalette,at=c(-8,-6,-4,-2,0,2,4),
##D   main="Basic GW regression coefficient estimates for LowEduc")
##D }
##D # Robust GW regression and map of the estimated coefficients for LowEduc
##D rgwr.res <- gwr.robust(GenEl2004~DiffAdd+LARent+SC1+Unempl+LowEduc+Age18_24
##D +Age25_44+Age45_64, data=Dub.voter,bw=bw.a,kernel="bisquare",
##D adaptive=TRUE,F123.test=TRUE)
##D print(rgwr.res)
##D if(require("RColorBrewer"))
##D {
##D   X11(width=10,height=12)
##D   spplot(rgwr.res$SDF, "LowEduc", key.space = "right",
##D   col.regions=mypalette,at=c(-8,-6,-4,-2,0,2,4),
##D   main="Robust GW regression coefficient estimates for LowEduc")
##D }
## End(Not run)



cleanEx()
nameEx("gwr.scalable")
### * gwr.scalable

flush(stderr()); flush(stdout())

### Name: gwr.scalable
### Title: Scalable GWR
### Aliases: gwr.scalable scgwr_pre scgwr_loocv scgwr_reg AICc1
###   gwr.scalable.loocv gwr_diag1 print.scgwrm
### Keywords: Scalable GWR

### ** Examples

## Not run: 
##D require(spData)
##D data(boston)
##D boston <- boston.c
##D coordinates(boston) <- ~ LON + LAT
##D res <- gwr.scalable(formula = MEDV ~ CRIM + ZN + INDUS + CHAS + AGE, data = boston, bw.adapt = 100)
##D res
## End(Not run)



cleanEx()
nameEx("gwss")
### * gwss

flush(stderr()); flush(stdout())

### Name: gwss
### Title: Geographically weighted summary statistics (GWSS)
### Aliases: gwss local.corr print.gwss
### Keywords: GWSS

### ** Examples

## Not run: 
##D data(EWHP)
##D data(EWOutline)
##D head(ewhp)
##D houses.spdf <- SpatialPointsDataFrame(ewhp[, 1:2], ewhp)
##D localstats1 <- gwss(houses.spdf, vars = c("PurPrice", "FlrArea"), bw = 50000)
##D head(data.frame(localstats1$SDF))
##D localstats1
##D ##A function for mapping data
##D if(require("RColorBrewer"))
##D {
##D    quick.map <- function(spdf,var,legend.title,main.title) 
##D    {
##D      x <- spdf@data[,var]
##D      cut.vals <- pretty(x)
##D      x.cut <- cut(x,cut.vals)
##D      cut.levels <- levels(x.cut)
##D      cut.band <- match(x.cut,cut.levels)
##D      colors <- brewer.pal(length(cut.levels), "YlOrRd")
##D      colors <- rev(colors)
##D      par(mar=c(1,1,1,1))
##D      plot(ewoutline,col="olivedrab",bg="lightblue1")
##D      title(main.title)
##D      plot(spdf,add=TRUE,col=colors[cut.band],pch=16)
##D      legend("topleft",cut.levels,col=colors,pch=16,bty="n",title=legend.title)
##D   }
##D   quick.map(localstats1$SDF, "PurPrice_LM", "1000's Uk Pounds", 
##D   "Geographically Weighted Mean")
##D   par(mfrow = c(1, 2))
##D   quick.map(localstats1$SDF, "PurPrice_LSKe", "Skewness Level", "Local Skewness")
##D   quick.map(localstats1$SDF, "PurPrice_LSD", "1000's Pounds", "Local Standard Deviation")
##D   #Exploring Non-Stationarity of Relationships
##D   quick.map(localstats1$SDF, "Corr_PurPrice.FlrArea", expression(rho), 
##D   "Geographically Weighted Pearson Correlation")
##D   #Robust, Quantile Based Local Summary Statistics
##D   localstats2 <- gwss(houses.spdf, vars = c("PurPrice", "FlrArea"), 
##D   bw = 50000, quantile = TRUE)
##D   quick.map(localstats2$SDF, "PurPrice_Median", "1000 UK Pounds", 
##D   "Geographically Weighted Median House Price")
##D }
## End(Not run)



cleanEx()
nameEx("gwss.montecarlo")
### * gwss.montecarlo

flush(stderr()); flush(stdout())

### Name: gwss.montecarlo
### Title: Monte Carlo (randomisation) test for gwss
### Aliases: gwss.montecarlo montecarlo.gwss
### Keywords: GWSS

### ** Examples

## Not run: 
##D data(LondonHP)
##D DM<-gw.dist(dp.locat=coordinates(londonhp))
##D test.lss<-gwss.montecarlo(data=londonhp, vars=c("PURCHASE","FLOORSZ"), bw=5000,
##D           kernel ="gaussian", dMat=DM,nsim=99)
##D test.lss
## End(Not run)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
