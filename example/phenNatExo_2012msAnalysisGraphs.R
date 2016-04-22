## Started 21 August 2012 ##
## By Lizzie ##


options(stringsAsFactors=FALSE)
library(nlme)
library(picante)
library(car)
library(ggplot2)
# library(caper)

setwd("/Users/Lizzie/Documents/R/NCEAS/Phenology/NatExotic")
source("/Users/Lizzie/Documents/R/fromNicolas/pglm.R")
source("phenNatExo_dataprep.R")
source("phenNatExo_helperfxs.R")

# get rid of ugly ggplot default theme
ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw()}


# some f(x)s to clean up graphic code mess
plotmodels <- function(filename, pglmodel, lmmodel, xrangeuse, colz, xlab){
pdf(paste("output/graphs/",filename,".pdf", sep=""), width=4, height=2, pointsize=12)
par(cex=cex,xpd=TRUE,yaxt="n")
tryme <- plot(xrangeuse,c(0, 2.25),type="n",
        xlab=xlab,
        ylab="")
# abline(h=0, col="gray")
text(xtxt,rep(-7.5,2),as.vector(c("exotic", "native")), srt=45, cex=1)
y<-c(1.5,1.5)
x<-as.vector(c(summary(lmmodel)$coef[1], summary(lmmodel)$coef[1]+
    summary(lmmodel)$coef[2]))
xsem<-as.vector(c(rep(summary(lmmodel)$coef[4], 2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(16), col=colz, cex=cexpoints)
# pglm
y<-c(0.75,0.75)
x<-as.vector(c(coef(pglmodel)[[1]], coef(pglmodel)[[1]]+coef(pglmodel)[[2]]))
xsem<-as.vector(c(rep(pglmodel$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
dev.off()
}

plotlmmodel <- function(filename, lmmodhere, xrangeuse, colz, xlab){
pdf(paste("output/graphs/",filename,".pdf", sep=""), width=4, height=2, pointsize=12)
par(cex=cex,xpd=TRUE,yaxt="n")
tryme <- plot(xrangeuse,c(0, 2.25),type="n",
        xlab=xlab,
        ylab="")
# abline(h=0, col="gray")
text(xtxt,rep(-7.5,2),as.vector(c("exotic", "native")), srt=45, cex=1)
y<-c(1.5,1.5)
x<-as.vector(c(summary(lmmodhere)$coef[1], summary(lmmodhere)$coef[1]+summary(lmmodhere)$coef[2]))
xsem<-as.vector(c(rep(summary(lmmodhere)$coef[4], 2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(16), col=colz, cex=cexpoints)
dev.off()
}

natexobyX <- function(data, xvar, yvar, xlim, ylim, col){
    natives <- subset(data, native.exotic=="native")
    exotics <- subset(data, native.exotic=="exotic")
    plot(data[[yvar]]~data[[xvar]], type="n", ylim=ylim, xlim=xlim)
    points(natives[[yvar]]~natives[[xvar]], col=col[2],
        ylim=ylim, xlim=xlim)
    points(exotics[[yvar]]~exotics[[xvar]], col=col[1],
        ylim=ylim,xlim=xlim)
  }



###########################
## Get the data & models ##
###########################

conc <- conctraits("input/sensitivities/concord_ALLSEAS.", "concord", conc.natexo, noxconc)
concdoy <- conctraits.doy(conc.natexo, noxconc, concchange)

conctref.mat <- conctrefum.mat
concgdd.dat <- subset(as.data.frame(conc$gdd[-341,]),is.na(native.exotic)==FALSE)
concgdd.dat <- concgdd.dat[-14,]

conctree <- matchtreephy(conctref.mat, concgdd.dat,"latbi")

conc.dat <- subset(as.data.frame(conc$mat[-341,]),is.na(native.exotic)==FALSE)
conc.dat <- conc.dat[-14,]

concdata <- matchtreedat(conctref.mat, conc.dat,"latbi",
    c("full.code", "nat.exoinv", "inv", "native.exotic", "native.exotic12", "inv.non12", "meanFFD", "daysperC"))
conctree <- matchtreephy(conctref.mat, conc.dat,"latbi") # overwriting above, but seems fine since species don't change

conc.doy <- subset(concdoy, is.na(native.exotic)==FALSE)
conc.doy <-  conc.doy[-13,]

concdoych <- matchtreedat(conctref.mat, conc.doy,"latbi",
    c("full.code", "nat.exoinv", "inv", "native.exotic", "doychange"))
conctreedoy <- matchtreephy(conctref.mat, conc.doy,"latbi")

vcv.conc <- vcv.phylo(conctree, model="Brownian", corr=FALSE)
vcv.doy <- vcv.phylo(conctreedoy, model="Brownian", corr=FALSE)
    
conc.ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=concdata)
conc.ffd.gee <-pglmEstLambda(meanFFD~as.factor(native.exotic), data=concdata, vcv.conc)
conc.doy.lm <- lm(doychange~as.factor(native.exotic), data=concdoych)
conc.doy.gee <- pglmEstLambda(doychange~as.factor(native.exotic), data=concdoych, vcv.doy) 
conc.mat.lm <- lm(daysperC~as.factor(native.exotic), data=concdata)
conc.mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=concdata, vcv.conc)

# bam, bam Konza!
konztref.mat <- konztrefum.mat
konz <- konzatraits("input/sensitivities/konza_kans_usa_ALLSEAS.", "konza", gentraitstra)
konzdoy <- konza.doy(gentraitstra,konzachange,  "konza")
konzmatic <- konzagoo.nosig(konz, konzdoy, konztref.mat)

vcv.konz <- vcv.phylo(konzmatic$tree.mat, model="Brownian", corr=FALSE)
    
konz.ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=konzmatic$matmod)
konz.ffd.gee <-pglmEstLambda(meanFFD~as.factor(native.exotic), 
    data=konzmatic$matmod, vcv.konz)
konz.mat.lm <- lm(daysperC~as.factor(native.exotic), data=konzmatic$matmod)
konz.mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
    data=konzmatic$matmod, vcv.konz)

konz.prec.lm <- lm(daysperMI~as.factor(native.exotic), data=konzmatic$gddpp)


# bam, bam Fargo!
fargotref.mat <- fargotrefum.mat
fargo <- fargotraits.wprecip.nocrazychenopodium("input/sensitivities/travers_nd_usa_ALLSEAS.", "fargo", far.natexo, noxfar)
fargodoy <- fargotraits.doy.nocrazychenopodium(far.natexo, noxfar, fargochange)
far.doy <-  subset(fargodoy ,is.na(native.exotic)==FALSE)
fardoy <-  subset(far.doy, native.exotic!="native/exotic" & native.exotic!="")
fargomatic <- mostsitesgoo.wprecip.nosig(fargo, fardoy, "doychange", fargotref.mat)

vcv.far <- vcv.phylo(fargomatic$tree.mat, model="Brownian", corr=FALSE)
vcv.fardoy <- vcv.phylo(fargomatic$treedoy, model="Brownian", corr=FALSE)
vcv.farpp <- vcv.phylo(fargomatic$tree.gddpp, model="Brownian", corr=FALSE)
    
far.ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=fargomatic$matmod)
far.ffd.gee <-pglmEstLambda(meanFFD~as.factor(native.exotic), data=fargomatic$matmod,
    vcv.far)
far.doy.lm <- lm(doychange~as.factor(native.exotic), data=fargomatic$doymod)
far.doy.gee <- pglmEstLambda(doychange~as.factor(native.exotic),
    data=fargomatic$doymod, vcv.fardoy)

far.mat.lm <- lm(daysperC~as.factor(native.exotic), data=fargomatic$matmod)
far.mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
    data=fargomatic$matmod, vcv.far)

far.prec.lm <- lm(daysperMI~as.factor(native.exotic), data=fargomatic$gddpp)
far.prec.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
    data=fargomatic$gddpp, vcv.farpp)


# bam, bam Chinnor!
fittref.mat <- fittrefum.mat
fitt <- fittertraits.wprecip("input/sensitivities/fitter_chn_grb_ALLSEAS.", "fitter", fittertra)
fittdoy <- fitter.doy(fittertra, fitterchange)
fittermatic <- mostsitesgoo.wprecip.nosig(fitt, fittdoy, "changeovertime",  fittref.mat)

vcv.fit <- vcv.phylo(fittermatic$tree.mat, model="Brownian", corr=FALSE)
vcv.fitdoy <- vcv.phylo(fittermatic$treedoy, model="Brownian", corr=FALSE)
vcv.fitgddpp <- vcv.phylo(fittermatic$tree.gddpp, model="Brownian", corr=FALSE)

fit.ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=fittermatic$matmod)
fit.ffd.gee <-pglmEstLambda(meanFFD~as.factor(native.exotic), data=fittermatic$matmod,
    vcv.fit)
fit.doy.lm <- lm(changeovertime~as.factor(native.exotic), data=fittermatic$doymod)
fit.doy.gee <- pglmEstLambda(changeovertime~as.factor(native.exotic),
    data=fittermatic$doymod, vcv.fitdoy)

fit.mat.lm <- lm(daysperC~as.factor(native.exotic), data=fittermatic$matmod)
fit.mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
    data=fittermatic$matmod, vcv.fit)

fit.pp.lm <- lm(daysperMI~as.factor(native.exotic), data=fittermatic$gddpp)
fit.pp.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
    data=fittermatic$gddpp, vcv.fitgddpp)

# bam, bam dumb DC
dctref.mat <- dctrefum.mat
sppdctraits <- subset(sppdctraits, native.exotic != "<NA>")
sppdctraits.ncul <- subset(sppdctraits, CUL=="N")
washdc <- washdctraits.wprecip("input/sensitivities/abu_wsh_usa_ALLSEAS.", "washdc", sppdctraits.ncul)
washdcdoy <- washdc.doy(sppdctraits.ncul, washdcchange)
washdcdoy <- subset(washdcdoy, is.na(changeovertime)==FALSE)
washdcmatic <- dcgoo.wprecip.nosig(washdc, washdcdoy, "changeovertime",   dctref.mat)

vcv.dc <- vcv.phylo(washdcmatic$tree.mat, model="Brownian", corr=FALSE)
vcv.dcdoy <- vcv.phylo(washdcmatic$treedoy, model="Brownian", corr=FALSE)
vcv.dcgddpp <- vcv.phylo(washdcmatic$tree.gddpp, model="Brownian", corr=FALSE)

dc.ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=washdcmatic$matmod)
dc.ffd.gee <-pglmEstLambda(meanFFD~as.factor(native.exotic), data=washdcmatic$matmod,
    vcv.dc)
dc.doy.lm <- lm(changeovertime~as.factor(native.exotic), data=washdcmatic$doymod)
dc.doy.gee <- pglmEstLambda(changeovertime~as.factor(native.exotic),
    data=washdcmatic$doymod, vcv.dcdoy)

dc.mat.lm <- lm(daysperC~as.factor(native.exotic), data=washdcmatic$matmod)
dc.mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
    data=washdcmatic$matmod, vcv.dc)

dc.pp.lm <- lm(daysperMI~as.factor(native.exotic), data=washdcmatic$gddpp)
dc.pp.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
    data=washdcmatic$gddpp, vcv.dcgddpp)

#########################
## Konza soil moisture ##
#########################

vcv.konzsma <- vcv.phylo(konzmatic$tree.sma, model="Brownian", corr=FALSE)
vcv.konzsmab <- vcv.phylo(konzmatic$tree.smab, model="Brownian", corr=FALSE)
vcv.konzsmanoout <- vcv.phylo(konzmatic$tree.smanoout, model="Brownian", corr=FALSE)


konz.sma.lm <- lm(daysperdev~as.factor(native.exotic), data=konzmatic$smamod)
konz.sma.gee <-pglmEstLambda(daysperdev~as.factor(native.exotic), 
    data=konzmatic$smamod, vcv.konzsma)

konz.smab.lm <- lm(daysperdev~as.factor(native.exotic), data=konzmatic$smabmod)
konz.smab.gee <-pglmEstLambda(daysperdev~as.factor(native.exotic), 
    data=konzmatic$smabmod, vcv.konzsmab)

konz.smanoout.lm <- lm(daysperdev~as.factor(native.exotic), data=konzmatic$smanooutmod)
konz.smanoout.gee <-pglmEstLambda(daysperdev~as.factor(native.exotic), 
    data=konzmatic$smanooutmod, vcv.konzsmanoout)

# merge files
konzsmxsheds <- merge(konzmatic$smamod, konzmatic$smabmod, by="row.names",
    suffixes=c(".1d", ".20b"))
konzsmxsheds$daysperdevxsheds <- (konzsmxsheds$daysperdev.1d+konzsmxsheds$daysperdev.20b)/2
row.names(konzsmxsheds) <- konzsmxsheds$Row.names

konzsmxsheds.o <- merge(konzmatic$smanooutmod, konzmatic$smabmod, by="row.names",
    suffixes=c(".1d", ".20b"))
konzsmxsheds.o$daysperdevxsheds <- (konzsmxsheds.o$daysperdev.1d+konzsmxsheds.o$daysperdev.20b)/2
row.names(konzsmxsheds.o) <- konzsmxsheds.o$Row.names

qqnorm(konzsmxsheds$daysperdev.20b)
qqnorm(konzsmxsheds$daysperdev.1d)
qqnorm(konzsmxsheds$daysperdevxsheds)

plot(daysperdev.1d~daysperdev.20b, data=konzsmxsheds)
summary(lm(daysperdev.1d~daysperdev.20b, data=konzsmxsheds))

konz.smxw.lm <-  lm(daysperdevxsheds~as.factor(native.exotic.1d), data=konzsmxsheds)
# konz.smxw.lmo <-  lm(daysperdevxsheds~as.factor(native.exotic.1d),
#    data=subset(konzsmxsheds,daysperdevxsheds<30)) # superseded by new smanoout file
konz.smxw.gee <-pglmEstLambda(daysperdevxsheds~as.factor(native.exotic.1d), 
    data=konzsmxsheds, vcv.konzsmab)

konz.smxwo.lm <-  lm(daysperdevxsheds~as.factor(native.exotic.1d), data=konzsmxsheds.o)
konz.smxwo.gee <-pglmEstLambda(daysperdevxsheds~as.factor(native.exotic.1d), 
    data=konzsmxsheds.o, vcv.konzsmanoout)

########################
## changeovertime & n ##
########################

conc.doyN <- merge(conc.doy , conc.dat, by=c("latbi", "native.exotic"))
summary(lm(doychange~n, data=conc.doyN)) # positive relationship
summary(lm(doychange~n*native.exotic, data=conc.doyN)) # all sig., intxn is (-0.8 value)

far.doyN <- merge(fardoy, fargo$mat, by=c("latbi", "native.exotic"))
summary(lm(doychange~n, data=far.doyN)) # no relationship (Anova)
summary(lm(doychange~n*native.exotic, data=far.doyN)) #intxn (-0.5 value)

fitt.doyN <- merge(fittdoy, fitt$mat, by=c("latbi", "native.exotic"))
summary(lm(changeovertime~n, data=fitt.doyN)) # slight negative
summary(lm(changeovertime~n*native.exotic, data=fitt.doyN)) # just main effects

wash.doyN <- merge(washdcdoy, washdc$mat, by=c("latbi", "native.exotic"))
summary(lm(changeovertime~n, data=wash.doyN)) # slight positive
summary(lm(changeovertime~n*native.exotic, data=wash.doyN)) # just main effect of n

#########################
## regression analyses ##
#########################

concdatadoy <- merge(conc.dat, conc.doy, by=c("latbi", "native.exotic"))

concdoych <- matchtreedat(conctref.mat, concdatadoy,"latbi",
    c("meanFFD","daysperC", "native.exotic", "doychange"))
conctreedatdoy <- matchtreephy(conctref.mat, concdatadoy,"latbi")

vcv.concReg <- vcv.phylo(conctreedatdoy, model="Brownian", corr=FALSE)

conc.fddmat.lm <- lm(daysperC~meanFFD, data=concdata) # sig
conc.fddmatINTER.lm <- lm(daysperC~meanFFD*native.exotic, data=concdata)
# meandFFD sig, native.exotic at 0.1 and intxn at 0.06
conc.fddmat.gee <- pglmEstLambda(daysperC~meanFFD, data=concdata, vcv.conc)

conc.match.lm <- lm(doychange~daysperC, data=concdoych) # sig
conc.matchINTER.lm <- lm(doychange~daysperC*native.exotic, data=concdoych)
    # sig inter  and sig effect of daysperC, exotics have steeper slope (as expected)
conc.match.gee <- pglmEstLambda(doychange~daysperC, data=concdoych, vcv.concReg)

conc.fddmat.gee$lambda # zero basically
conc.match.gee$lambda # zero basically


# konza
konz.fddmat.lm <- lm(daysperC~meanFFD, data=konzmatic$matmod) # sig
konz.fddmatINTER.lm <- lm(daysperC~meanFFD*native.exotic, data=konzmatic$matmod) # meanFFD sig
konz.fddmat.gee <- pglmEstLambda(daysperC~meanFFD, data=konzmatic$matmod, vcv.konz)

# fargo
fargomatic$mat[["latbi"]] <- row.names(fargomatic$mat)
fardatadoy <- merge(fargomatic$mat, fardoy , by=c("latbi", "native.exotic"))

fardoych <- matchtreedat(fargotref.mat, fardatadoy,"latbi",
    c("meanFFD","daysperC", "native.exotic", "doychange"))
fartreedatdoy <- matchtreephy(fargotref.mat, fardatadoy,"latbi")

vcv.farReg <- vcv.phylo(fartreedatdoy, model="Brownian", corr=FALSE)

far.fddmat.lm <- lm(daysperC~meanFFD, data=fargomatic$matmod) # NS
far.fddmatINTER.lm <- lm(daysperC~meanFFD*native.exotic, data=fargomatic$matmod) # all NS
far.fddmat.gee <- pglmEstLambda(daysperC~meanFFD, data=fargomatic$matmod, vcv.far)

far.match.lm <- lm(doychange~daysperC, data=fardoych) # sig
far.matchINTER.lm <- lm(doychange~daysperC*native.exotic, data=fardoych) # sig intxn and sens
far.match.gee <- pglmEstLambda(doychange~daysperC, data=fardoych, vcv.farReg)

far.fddmat.gee$lambda # zero basically
far.match.gee$lambda # zero basically

# chinnor
fittermatic$mat[["latbi"]] <- row.names(fittermatic$mat)
fitdatadoy <- merge(fittermatic$mat, fittdoy , by=c("latbi", "native.exotic"))

fitdoych <- matchtreedat(fittref.mat, fitdatadoy,"latbi",
    c("meanFFD","daysperC", "native.exotic", "changeovertime"))
fittreedatdoy <- matchtreephy(fittref.mat, fitdatadoy,"latbi")

vcv.fitReg <- vcv.phylo(fittreedatdoy, model="Brownian", corr=FALSE)

fit.fddmat.lm <- lm(daysperC~meanFFD, data=fittermatic$matmod) # sig
fit.fddmatINTER.lm <- lm(daysperC~meanFFD*native.exotic, data=fittermatic$matmod)
   # main effects sig
fit.fddmat.gee <- pglmEstLambda(daysperC~meanFFD, data=fittermatic$matmod, vcv.fit)

fit.match.lm <- lm(changeovertime~daysperC, data=fitdoych) # sig
fit.matchINTER.lm <- lm(changeovertime~daysperC*native.exotic, data=fitdoych) # sig daysperC only
fit.match.gee <- pglmEstLambda(changeovertime~daysperC, data=fitdoych, vcv.fitReg)

fit.fddmat.gee$lambda # 0.075
fit.match.gee$lambda # zero basically

coef(fit.fddmat.lm)[["meanFFD"]]
coef(fit.fddmat.gee)[["meanFFD"]]

# DC
washdcmatic$mat[["latbi"]] <- row.names(washdcmatic$mat)
dcdatadoy <- merge(washdcmatic$mat, washdcdoy , by=c("latbi", "native.exotic"))

dcdoych <- matchtreedat(dctref.mat, dcdatadoy,"latbi",
    c("meanFFD","daysperC", "native.exotic", "changeovertime"))
dctreedatdoy <- matchtreephy(dctref.mat, dcdatadoy,"latbi")

vcv.dcReg <- vcv.phylo(dctreedatdoy, model="Brownian", corr=FALSE)

dc.fddmat.lm <- lm(daysperC~meanFFD, data=washdcmatic$matmod) # sig
dc.fddmatINTER.lm <- lm(daysperC~meanFFD*native.exotic, data=washdcmatic$matmod)
    # sig meanFFD and natexo at 0.09
dc.fddmat.gee <- pglmEstLambda(daysperC~meanFFD, data=washdcmatic$matmod, vcv.dc)

dc.match.lm <- lm(changeovertime~daysperC, data=dcdoych) # sig
dc.matchINTER.lm <- lm(changeovertime~daysperC*native.exotic, data=dcdoych)
    # sig all (0.08 for native.exotic)
dc.match.gee <- pglmEstLambda(changeovertime~daysperC, data=dcdoych, vcv.dcReg)

dc.fddmat.gee$lambda # 0.076
dc.match.gee$lambda # zero basically

coef(dc.fddmat.lm)[["meanFFD"]]
coef(dc.fddmat.gee)[["meanFFD"]]

##########################
## Abundance at Concord ##
##########################

concabundoy <- merge(conc.doy, concAbunsm, by="latbi")
concdoychabun <- matchtreedat(conctref.mat, concabundoy,"latbi",
     c("full.code", "nat.exoinv", "inv", "native.exotic",
     "doychange", "DeltaAbun_cont"))

conctreedoychabun <- matchtreephy(conctref.mat, concabundoy,"latbi")

concdoychabun$abun_in3 <- as.numeric(concdoychabun$DeltaAbun_cont)
concdoychabun$abun_in3[concdoychabun$abun_in3>0 & concdoychabun$abun_in3<7] <- 1
concdoychabun$abun_in3[concdoychabun$abun_in3<0 & concdoychabun$abun_in3>-7] <- -1

vcv.concabun <- vcv.phylo(conctreedoychabun, model="Brownian", corr=FALSE)

conc.doyabun3.lm <- lm(doychange~as.factor(native.exotic)*as.factor(abun_in3), data=concdoychabun)
conc.doyabun3.gee <- pglmEstLambda(doychange~as.factor(native.exotic)*as.factor(abun_in3),
    data=concdoychabun, vcv.concabun) # lambda close to zero!

conc.doyabun.lm <- lm(doychange~as.factor(native.exotic)*as.numeric(DeltaAbun_cont),
    data= concdoychabun)


###################
## Pooling sites ##
###################

concmatmod <- subset(concdata, select=c("native.exotic","native.exotic12","meanFFD",
    "daysperC"))
concmatmod$site <- "concord"

fittmatmod <- subset(as.data.frame(fittermatic$matmod), 
    select=c("native.exotic","native.exotic12","meanFFD","daysperC"))
fittmatmod$site <- "fitter"
farmatmod <- subset(as.data.frame(fargomatic$matmod),
    select=c("native.exotic","native.exotic12","meanFFD","daysperC"))
farmatmod$site <- "fargo"
konzmatmod <- subset(as.data.frame(konzmatic$matmod),
    select=c("native.exotic","native.exotic12","meanFFD","daysperC"))
konzmatmod$site <- "konza"
washdcmatmod <- subset(as.data.frame(washdcmatic$matmod),
    select=c("native.exotic","native.exotic12","meanFFD","daysperC"))
washdcmatmod$site <- "dc"

allsitesgoo <- rbind(fittmatmod, concmatmod, farmatmod, konzmatmod,
     washdcmatmod)
allsitesgoo$latbi <- row.names(allsitesgoo)

concdoych$site <- "concord"
names(concdoych)[names(concdoych)=="doychange"] <- "changeovertime"
fardoych$site <- "fargo"
names(fardoych)[names(fardoych)=="doychange"] <- "changeovertime"
fitdoych$site <- "fitter"
dcdoych$site <- "washdc"

allsitesdoy <- rbind(concdoych, fardoych, fitdoych, dcdoych)
allsitesdoy$latbi <- row.names(allsitesdoy)

##
FFDsitesmod <- lm(meanFFD~native.exotic*site, data=allsitesgoo)
FFDsitememod <- lme(meanFFD~native.exotic, random=~1|site, data=allsitesgoo)

senssitesmod <- lm(daysperC~native.exotic*site, data=allsitesgoo)
senssitesmemod <- lme(daysperC~native.exotic, random=~1|site, data=allsitesgoo)

doysitesmod <- lm(changeovertime~native.exotic*site, data=allsitesdoy)
doysitesmemod <- lme(changeovertime~native.exotic, random=~1|site, data=allsitesdoy)

FFDXsens.sitesmod <- lm(daysperC~meanFFD*site, data=allsitesgoo)
FFDXsens.sitesmemod <- lme(daysperC~meanFFD, random=~1|site, data=allsitesgoo)

sensXdoy.sitesmod <- lm(changeovertime~daysperC*site, data=allsitesdoy)
sensXdoy.sitesmemod <- lme(changeovertime~daysperC, random=~1|site, data=allsitesdoy)



####################################################
## Comparing Chinnor & Concord: 42 shared species ##
####################################################

concgdd <- conc$gdd
fittgdd <- fitt$gdd
fittgddsm <- subset(fittgdd, select=c("latbi", "daysperHI", "meanFFD",
    "native.exotic"))
concgddsm <- subset(concgdd, select=c("latbi", "daysperHI", "meanFFD",
    "native.exotic"))

spoverlap <- merge(fittgddsm, concgddsm, by="latbi", suffixes=c(".fitt", ".conc"))

spoverlap$natexo.sq <- paste(spoverlap$native.exotic.fitt,
    spoverlap$native.exotic.conc, sep="-")

nativesinFitt <- subset(spoverlap, natexo.sq=="native-exotic") # 42 species!

fittgdd$nativesinFitt <- "no"
fittgdd$nativesinFitt[which(fittgdd$latbi %in% nativesinFitt$latbi)] <- "yes"

concgdd$nativesinFitt <- "no"
concgdd$nativesinFitt[which(concgdd$latbi %in% nativesinFitt$latbi)] <- "yes"

# get rid of other exotics at Concord and all exotics at Fitter:
concgdd <- concgdd[which(!(concgdd$nativesinFitt=="no"
    & concgdd$native.exotic=="exotic")),]
fittgdd <- subset(fittgdd, native.exotic=="native")

concdoy <- conctraits.doy(conc.natexo, noxconc, concchange)
fittdoy <- fitter.doy(fittertra, fitterchange)

fittdoy$nativesinFitt <- "no"
fittdoy$nativesinFitt[which(fittdoy$latbi %in% nativesinFitt$latbi)] <- "yes"

concdoy$nativesinFitt <- "no"
concdoy$nativesinFitt[which(concdoy$latbi %in% nativesinFitt$latbi)] <- "yes"

# get rid of other exotics at Concord and all exotics at Fitter:
concdoy <- concdoy[which(!(concdoy$nativesinFitt=="no"
    & concdoy$native.exotic=="exotic")),]
fittdoy <- subset(fittdoy, native.exotic=="native")

spo.fittmod <- lm(daysperHI~nativesinFitt, data=fittgdd)
spo.fittmod.ffd <- lm(meanFFD~nativesinFitt, data=fittgdd)
spo.fittmod.doy <- lm(changeovertime~nativesinFitt, data=fittdoy)
spo.concmod <- lm(daysperHI~nativesinFitt, data=concgdd)
spo.concmod.ffd <- lm(meanFFD~nativesinFitt, data=concgdd)
spo.concmod.doy <- lm(doychange~nativesinFitt, data=concdoy)

# get the 42 species and see if temperature sensitivities differ
nativesinFitt$diff <- nativesinFitt$daysperHI.fitt -  nativesinFitt$daysperHI.conc
binom.test(c(nrow(subset(nativesinFitt, diff<0)), nrow(subset(nativesinFitt, diff>0))))



########################
### Graphs are here! ###
########################

colourz <- c("firebrick", "dodgerblue4")
colz<- c("firebrick", "dodgerblue2")
cex <- 0.8
cexpoints <- 1.5
xtxt <- c(1,2)

#####
## X~Y regression plots
#####

## Concord
pdf(paste("output/graphs/concXsensFFD.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(concdata, "meanFFD", "daysperC", c(50,280), c(-45,35), colz)
abline(conc.fddmat.lm)
dev.off()

pdf(paste("output/graphs/concXsensFFDitxn.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(concdata, "meanFFD", "daysperC", c(50,280), c(-45,35), colz)
abline(coef(conc.fddmat.lm)[1], coef(conc.fddmat.lm)[2], col=colz[1])
abline(coef(conc.fddmat.lm)[1], coef(conc.fddmat.lm)[2]+coef(conc.fddmat.lm)[4], col=colz[2])
dev.off()

pdf(paste("output/graphs/concXchangesens.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(concdoych, "daysperC","doychange",  c(-40,30), c(-110, 60), colz)
abline(conc.match.lm)
dev.off()

pdf(paste("output/graphs/concXchangesensitxn.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(concdoych, "daysperC","doychange",  c(-40,30), c(-110, 60), colz)
abline(coef(conc.matchINTER.lm)[1], coef(conc.matchINTER.lm)[2], col=colz[1])
abline(coef(conc.matchINTER.lm)[1], coef(conc.matchINTER.lm)[2]+coef(conc.matchINTER.lm)[4], col=colz[2])
dev.off()


## Konza
pdf(paste("output/graphs/konzXsensFFD.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(as.data.frame(konzmatic$matmod), "meanFFD", "daysperC", c(50,280),
    c(-25,15), colz)
abline(konz.fddmat.lm)
dev.off()

## Fargo
pdf(paste("output/graphs/farXsensFFD.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(as.data.frame(fargomatic$matmod), "meanFFD", "daysperC", c(80,280),
    c(-30,20), colz)
# abline(far.fddmat.lm)
dev.off()

pdf(paste("output/graphs/farXchangesens.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(fardoych, "daysperC","doychange", c(-15, 20), c(-35, 30), colz)
abline(far.match.lm)
dev.off()

pdf(paste("output/graphs/farXchangesensitxn.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(fardoych, "daysperC","doychange", c(-15, 20), c(-35, 30), colz)
abline(coef(far.matchINTER.lm)[1], coef(far.matchINTER.lm)[2], col=colz[1])
abline(coef(far.matchINTER.lm)[1], coef(far.matchINTER.lm)[2]+coef(far.matchINTER.lm)[4], col=colz[2])
dev.off()

## Fitter (Chinnor)
pdf(paste("output/graphs/chinnXsensFFD.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(as.data.frame(fittermatic$matmod), "meanFFD", "daysperC", c(10,280),
    c(-40,20), colz)
abline(fit.fddmat.lm)
dev.off()

pdf(paste("output/graphs/chinnXchangesens.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(fitdoych, "daysperC","changeovertime", c(-40, 20), c(-2, 1), colz)
abline(fit.match.lm)
dev.off()


pdf(paste("output/graphs/chinnXchangesensitxn.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(fitdoych, "daysperC","changeovertime", c(-40, 20), c(-2, 1), colz)
abline(coef(fit.matchINTER.lm)[1], coef(fit.matchINTER.lm)[2], col=colz[1])
abline(coef(fit.matchINTER.lm)[1], coef(fit.matchINTER.lm)[2]+coef(fit.matchINTER.lm)[4], col=colz[2])
dev.off()



## DC
pdf(paste("output/graphs/dcXsensFFD.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(as.data.frame(washdcmatic$matmod), "meanFFD", "daysperC", c(30,160),
    c(-60,40), colz)
abline(dc.ffdmat.lm)
dev.off()

pdf(paste("output/graphs/dcXchangesens.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(dcdoych, "daysperC","changeovertime", c(-60, 30), c(-17,10), colz)
abline(dc.match.lm)
dev.off()

pdf(paste("output/graphs/dcXchangesensitxn.pdf", sep=""), width=6, height=4, pointsize=12)
natexobyX(dcdoych, "daysperC","changeovertime", c(-60, 30), c(-17,10), colz)
abline(coef(dc.matchINTER.lm)[1], coef(dc.matchINTER.lm)[2], col=colz[1])
abline(coef(dc.matchINTER.lm)[1], coef(dc.matchINTER.lm)[2]+coef(dc.matchINTER.lm)[4], col=colz[2])
dev.off()


#####
## Basic LM AND PGLM model plots and histograms
#####

plotmodels("concMAT", conc.mat.gee, conc.mat.lm, c(-6,0), colz,
    "sensitivity (changes in days per C")
plotmodels("concFFD", conc.ffd.gee, conc.ffd.lm, c(145,185), colz,
    "mean flowering day of year")
plotmodels("concchange", conc.doy.gee, conc.doy.lm, c(-8,3), colz,
    "change in flowering over time") 

plotmodels("konzMAT", konz.mat.gee, konz.mat.lm, c(-6,0), colz,
    "sensitivity (changes in days per C")
plotmodels("konzFFD", konz.ffd.gee, konz.ffd.lm, c(110,185), colz,
    "mean flowering day of year")
plotmodels("konzSM", konz.smxw.gee, konz.smxw.lm, c(-5,5), colz,
    "sensitivity (change in days per standardized soil moisture unit)")

plotmodels("fargoMAT", far.mat.gee, far.mat.lm, c(-2.5,0), colz,
    "sensitivity (changes in days per C")
plotmodels("fargoFFD", far.ffd.gee, far.ffd.lm, c(150,165), colz,
    "mean flowering day of year")
plotmodels("fargochange", far.doy.gee, far.doy.lm, c(-2,6), colz,
    "change in flowering over time")
plotmodels("fargoPREC", far.prec.gee, far.prec.lm, c(-3.5,3.5), colz,
    "sensitivity (change in days per standardized precipitation unit)")

plotmodels("chinnMAT", fit.mat.gee, fit.mat.lm, c(-12,-5), colz,
    "sensitivity (changes in days per C")
plotmodels("chinnFFD", fit.ffd.gee, fit.ffd.lm, c(120, 155), colz,
    "mean flowering day of year")
plotmodels("chinnchange", fit.doy.gee, fit.doy.lm, c(-0.3,0.1), colz,
    "change in flowering over time")
plotmodels("chinPREC", fit.pp.gee, fit.pp.lm, c(-0.3,1.5), colz,
    "sensitivity (precipitation)")

plotmodels("dcMAT", dc.mat.gee, dc.mat.lm, c(-5,-2), colz, 
    "sensitivity (changes in days per C")
plotmodels("dcFFD", dc.ffd.gee, dc.ffd.lm, c(110, 125), colz,
    "mean flowering day of year")
plotmodels("dcchange", dc.doy.gee, dc.doy.lm, c(-0.3,0.1), colz,
    "change in flowering over time")
plotmodels("dcPREC", dc.pp.gee, dc.pp.lm, c(-0.3,3), colz,
    "sensitivity (precipitation)")


##
## histograms 
##
daterhere <- concdata
pdf(paste("output/graphs/concMAThistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=daysperC), xlab="sensitivity (change in days per deg C)") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

pdf(paste("output/graphs/concFFDhistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=meanFFD), xlab="mean first flowering day of year") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

daterhere <- concdoych
pdf(paste("output/graphs/concchangehistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=doychange), xlab="change in flowering") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

# konza
daterhere <- konzmatic$matmod
pdf(paste("output/graphs/konzMAThistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=daysperC), xlab="sensitivity (change in days per deg C)") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

pdf(paste("output/graphs/konzFFDhistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=meanFFD), xlab="mean first flowering day of year") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

daterhere <- konzsmxsheds
pdf(paste("output/graphs/konzSMhistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=daysperdevxsheds), xlab="sensitivity (change in days per standardized soil moisture unit)") + 
    geom_histogram(data=subset(daterhere, native.exotic.1d=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic.1d=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()


# fargo
daterhere <- fargomatic$matmod
pdf(paste("output/graphs/fargoMAThistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=daysperC), xlab="sensitivity (change in days per deg C)") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

pdf(paste("output/graphs/fargoFFDhistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=meanFFD), xlab="mean first flowering day of year") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

daterhere <- fargomatic$doymod
pdf(paste("output/graphs/fargochangehistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=doychange), xlab="change in flowering") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

daterhere <- fargomatic$gddpp
pdf(paste("output/graphs/fargoPREChistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=daysperMI), xlab="sensitivity (change in days per standardized precipitation unit)")+ 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) +
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

# chinnor
daterhere <- fittermatic$matmod
pdf(paste("output/graphs/chinnorMAThistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=daysperC), xlab="sensitivity (change in days per deg C)") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

pdf(paste("output/graphs/chinnorFFDhistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=meanFFD), xlab="mean first flowering day of year") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

daterhere <- fittermatic$doymod
pdf(paste("output/graphs/chinnorchangehistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=changeovertime), xlab="change in flowering") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

daterhere <- fittermatic$gddpp
pdf(paste("output/graphs/chinnorPREChistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=daysperMI), xlab="sensitivity (change in days per precipitation unit)") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()



# DC
daterhere <- washdcmatic$matmod
pdf(paste("output/graphs/dcMAThistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=daysperC), xlab="sensitivity (change in days per deg C)") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

pdf(paste("output/graphs/dcFFDhistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=meanFFD), xlab="mean first flowering day of year") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

daterhere <- washdcmatic$doymod
pdf(paste("output/graphs/dcchangehistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=changeovertime), xlab="change in flowering") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()

daterhere <- washdcmatic$gddpp
pdf(paste("output/graphs/dcPREChistograms.pdf", sep=""), width=4, height=3.5, pointsize=12)
ggplot(daterhere, aes(x=daysperMI), xlab="change in flowering") + 
    geom_histogram(data=subset(daterhere, native.exotic=="exotic"), 
        fill=colourz[1], alpha = 0.8) + 
    geom_histogram(data=subset(daterhere, native.exotic=="native"), 
        fill=colourz[2], alpha = 0.5)
dev.off()



##############################
## lambda stuff 23 Oct 2012 ##
##############################

fitContinuous(conctree, concdata$meanFFD, model="lambda")
 # for phylomatic: 0.9455048, for phlawd: 0.1218433
lambdaTree(conctree, 0) -> concLambda0
fitContinuous(concLambda0, concdata$meanFFD)

##############################
## all pglms on one graph  ##
##############################

## sensitivities
##
xrangeuse <- c(-12,1)
pdf("output/graphs/pglmsSens.pdf", width=4, height=6, pointsize=12)
par(cex=cex,xpd=TRUE,yaxt="n")
tryme <- plot(xrangeuse,c(0.1, 0.7),type="n",
        xlab=xlab,
        ylab="")
# abline(h=0, col="gray")
text(xtxt,rep(-13,2),as.vector(c("exotic", "native")), srt=45, cex=1)
# pglm chinnor
y<-c(0.2,0.2)
x<-as.vector(c(coef(fit.mat.gee)[[1]], coef(fit.mat.gee)[[1]]+coef(fit.mat.gee)[[2]]))
xsem<-as.vector(c(rep(fit.mat.gee$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
# pglm concord
y<-c(0.3,0.3)
x<-as.vector(c(coef(conc.mat.gee)[[1]], coef(conc.mat.gee)[[1]]+coef(conc.mat.gee)[[2]]))
xsem<-as.vector(c(rep(conc.mat.gee$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
# pglm Washdc
y<-c(0.4, 0.4)
x<-as.vector(c(coef(dc.mat.gee)[[1]], coef(dc.mat.gee)[[1]]+coef(dc.mat.gee)[[2]]))
xsem<-as.vector(c(rep(dc.mat.gee$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
# pglm Fargo
y<-c(0.5,0.5)
x<-as.vector(c(coef(far.mat.gee)[[1]], coef(far.mat.gee)[[1]]+coef(far.mat.gee)[[2]]))
xsem<-as.vector(c(rep(far.mat.gee$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
# pglm Konza
y<-c(0.6,0.6)
x<-as.vector(c(coef(konz.mat.gee)[[1]], coef(konz.mat.gee)[[1]]+coef(konz.mat.gee)[[2]]))
xsem<-as.vector(c(rep(konz.mat.gee$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
dev.off()

## changes
##
xrangeuse <- c(-8,6)
pdf("output/graphs/pglmsChangeDays.pdf", width=4, height=3, pointsize=12)
par(cex=cex,xpd=TRUE,yaxt="n")
tryme <- plot(xrangeuse,c(0.25, 0.45),type="n",
        xlab=xlab,
        ylab="")
# abline(h=0, col="gray")
text(xtxt,rep(-13,2),as.vector(c("exotic", "native")), srt=45, cex=1)
# pglm concord
y<-c(0.3,0.3)
x<-as.vector(c(coef(conc.doy.gee)[[1]], coef(conc.doy.gee)[[1]]+coef(conc.doy.gee)[[2]]))
xsem<-as.vector(c(rep(conc.doy.gee$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
# pglm Fargo
y<-c(0.4,0.4)
x<-as.vector(c(coef(far.doy.gee)[[1]], coef(far.doy.gee)[[1]]+coef(far.doy.gee)[[2]]))
xsem<-as.vector(c(rep(far.doy.gee$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
dev.off()

xrangeuse <- c(-0.3,0.21)
pdf("output/graphs/pglmsChangeYrs.pdf", width=4, height=3, pointsize=12)
par(cex=cex,xpd=TRUE,yaxt="n")
tryme <- plot(xrangeuse,c(0.1, 0.4),type="n",
        xlab=xlab,
        ylab="")
# abline(h=0, col="gray")
text(xtxt,rep(-13,2),as.vector(c("exotic", "native")), srt=45, cex=1)
# pglm chinnor
y<-c(0.2,0.2)
x<-as.vector(c(coef(fit.doy.gee)[[1]], coef(fit.doy.gee)[[1]]+coef(fit.doy.gee)[[2]]))
xsem<-as.vector(c(rep(fit.doy.gee$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
# pglm Washdc
y<-c(0.3, 0.3)
x<-as.vector(c(coef(dc.doy.gee)[[1]], coef(dc.doy.gee)[[1]]+coef(dc.doy.gee)[[2]]))
xsem<-as.vector(c(rep(dc.doy.gee$sterr[2],2)))
arrows(x-xsem,y,x+xsem,y,code=3,angle=90,length=0.05, col=colz)
points(x,y, pch=c(17), col=colz, cex=cexpoints)
dev.off()

## note to self on imposable histograms:
# library(lattice)
# histogram(~dcme$changeovertime|dcme$native.exotic)
