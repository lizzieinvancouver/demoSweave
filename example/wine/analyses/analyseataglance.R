### Started 9.December.2014 ###
### By Lizzie ###

## This file does a quick analysis of the 'at a glance' vintage data
# from Broadbent's book ##

## Updated 17 March 2016 ##
## To correct broken paths and get mean timing of harvest ##

## Updated 16 June 2015 ##
## To reference the correct climate data files ##

## It uses the at a glance data with data added in for others years by Jehane #
# but I still need to finish cleaning the at a glance a little more ##
# So I should update this file when I have done that (see cleanataglance.R). ##

##################

# housekeeping #
rm(list=ls()) 
options(stringsAsFactors=FALSE)

library(ggplot2)
library(lattice)
library(reshape)
library(ordinal)
library(xtable)

# stop("stopping here Lizzie ...")

setwd("~/Documents/git/teaching/demoSweave/example/wine/analyses")
source("source/analysesfunctions.R")

# data #
# vintages updated 10 April 2015 #
# climate data updated 16 June 2015 #
# vintage data
allstars4reg <- read.csv("data/output/allstars4reg.csv", header=TRUE)
# climate data from Ben
temp <- read.csv("data/climate/seas_temp_MJJ.onedeg.csv",
    header=TRUE)
colnames(temp)[1] <- "year"
prec <- read.csv("data/climate/seas_prec_MJJ.onedeg.csv",
    header=TRUE)
colnames(prec)[1] <- "year"
pdsi <- read.csv("data/climate/seas_pdsi_MJJ.onedeg.csv",
    header=TRUE)
colnames(pdsi)[1] <- "year"
daux <- read.csv("data/input/dauxdata.csv", header=TRUE, skip=2)
daux <- daux[,1:28]
names(daux)[names(daux)=="Abb."] <- "year"
ghd <- as.data.frame(daux)
ghd <- subset(ghd, select=c("year", "Bor", "Bur"))

envdata <- merge(ghd, temp, by="year", suffixes=c(".ghd", ".temp"))
envdata <- merge(envdata, prec, by="year", suffixes=c("", ".prec"))
envdata <- merge(envdata, pdsi, by="year", suffixes=c("", ".pdsi"))

# 1600-1980 mean timing for Bordeaux and Burgundy #
dauxearlytimes <- subset(daux, year>1599 & year<1981)
dauxlatertimes <- subset(daux, year>1980)

mean(dauxearlytimes$Bur, na.rm=TRUE) # 27 so September 28
mean(dauxlatertimes$Bur, na.rm=TRUE) # 19 so September 20

mean(dauxearlytimes$Bor, na.rm=TRUE) # 25 so September 26
mean(dauxlatertimes$Bor, na.rm=TRUE) # 16.5 so September 17 or 18

subset(dauxlatertimes, year==2003) # Bur was -10 (so Aug 21); Bor was -4.3 or August 28

# break out data by region for analyses, more plots # 
burgbord <- subset(allstars4reg, region=="burgundy" | region=="bordeaux")
bbshape <- recast(burgbord, id.var=c("year", "colour.name", "region"), 
    measure.var="stars", formula=year~colour.name+region)

# How related are things across regions? #
summary(lm(red_burgundy~white_burgundy, bbshape)) # R2 of 0.5
summary(lm(red_bordeaux~white_bordeaux, bbshape)) # R2 of 0.4
summary(lm(red_burgundy~red_bordeaux, bbshape)) # R2 of 0.4
summary(lm(white_burgundy~white_bordeaux, bbshape)) # R2 of 0.4

# merge in environmental data, and away we go #
bbshapetemp <- merge(temp, bbshape, by="year")
bbshapeprec <- merge(prec, bbshape, by="year")
bbshapepdsi <- merge(pdsi, bbshape, by="year")
bbshapeghd <- merge(ghd, bbshape, by="year")

bbshapeenv <- merge(bbshape, envdata, by="year")

# summaries #
# set year
yearhere <- 1980

# this steps through the CLM by temp, precop and pdsi
# it's what I did before I built before the f(x) in the source file (above)
bordtbf80 <- clm(as.factor(red_bordeaux)~Bor, data=
    subset(bbshapetemp, year<=yearhere))
bordtaft80 <- clm(as.factor(red_bordeaux)~Bor, data=
    subset(bbshapetemp, year>yearhere))
bordprbf80 <- clm(as.factor(red_bordeaux)~Bor, data=
    subset(bbshapeprec, year<=yearhere))
bordpraft80 <- clm(as.factor(red_bordeaux)~Bor, data=
    subset(bbshapeprec, year>yearhere))
bordpdbf80 <- clm(as.factor(red_bordeaux)~Bor, data=
    subset(bbshapepdsi, year<=yearhere))
bordpdaft80 <- clm(as.factor(red_bordeaux)~Bor, data=
    subset(bbshapepdsi, year>yearhere))

bordredrow <-cbind(bordtbf80$coefficients["Bor"],
    bordtaft80$coefficients["Bor"], bordprbf80$coefficients["Bor"],
    bordpraft80$coefficients["Bor"], bordpdbf80$coefficients["Bor"],
    bordpdaft80$coefficients["Bor"])

##################
## CLM models ####
##################

## alert! The CLM f(x)s: GHD comes first, then the climate vars

rbor <- getclmrows("red_bordeaux", "Bor", 1980)
wbor <- getclmrows("white_bordeaux", "Bor", 1980)
rbur <- getclmrows("red_burgundy", "Bur", 1980)
wbur <- getclmrows("white_burgundy", "Bur", 1980)

clm.mods <- rbind(rbor, wbor, rbur, wbur)

clm.models <- data.frame(clm.mods , row.names = c("Red Bordeaux",
    "White Bordeaux", "Red Burgundy", "White Burgundy"))
names(clm.models) <- c("GHD: pre", "GHD: post", "Temp: pre-80",
    "Temp: post-80", "Prec: pre-80", "Prec: post-80", "PDSI: pre-80",
    "PDSI: post-80")

## CLM models p values 
rborp <- getclmrows.wpval("red_bordeaux", "Bor", 1980)
wborp <- getclmrows.wpval("white_bordeaux", "Bor", 1980)
rburp <- getclmrows.wpval("red_burgundy", "Bur", 1980)
wburp <- getclmrows.wpval("white_burgundy", "Bur", 1980)

clm.mods.wpval <- rbind(rborp, wborp, rburp, wburp)

clm.models.wpval <- data.frame(clm.mods.wpval , row.names = c("Red Bordeaux",
    "White Bordeaux", "Red Burgundy", "White Burgundy"))
names(clm.models.wpval) <- c("GHD: pre", "GHD: post", "Temp: pre-80", 
    "Temp: post-80", "Prec: pre-80", "Prec: post-80", "PDSI: pre-80",
    "PDSI: post-80")

clm.pvals.simple <- clm.models.wpval
for (i in c(1:8)) {
    clm.pvals.simple[,i] <- ifelse(clm.pvals.simple[,i] >0.01,
        round(clm.pvals.simple[,i], 2), "<0.01")
}

clm.models.full <- clm.models
for (j in c(1:4, 7:8)) {
clm.models.full[,j] <- paste(round(clm.models[,j],3), " (",
    clm.pvals.simple[,j],")", sep=(""))
}
for (j in c(5:6)) {
clm.models.full[,j] <- paste(round(clm.models[,j],3), " (",
    clm.pvals.simple[,j],")", sep=(""))
}

clm.models.part1 <- clm.models.full[,1:4]
names(clm.models.part1) <- c("GHD: 1900-1980", "GHD: 1981-2001",
    "Temp: 1900-1980", "Temp: 1981-2001")
clm.models.part2 <- clm.models.full[,5:8]
names(clm.models.part2) <- c("Prec: 1900-1980", "Prec: 1981-2001",
    "PDSI: 1900-1980", "PDSI: 1981-2001")

####################
## lm regressions ##

rborp <- getlmrows.wpval("red_bordeaux", "Bor", 1979)
wborp <- getlmrows.wpval("white_bordeaux", "Bor", 1979)
rburp <- getlmrows.wpval("red_burgundy", "Bur", 1979)
wburp <- getlmrows.wpval("white_burgundy", "Bur", 1979)
lm.mods.wpval <- rbind(rborp, wborp, rburp, wburp)

##
## plots!
colz3 <- c("dodgerblue", "darkturquoise", "firebrick3")

threebreaksplot <- function(dater, whichwine, whichx, colpalette, maintitle){
    plot(dater[[whichwine]]~dater[[whichx]], type="n", ylab="Broadbent quality",
        xlab=whichx, main=maintitle)
    data1950=subset(dater, year<1949)
    points(data1950[[whichwine]]~data1950[[whichx]], pch=16, col=colpalette[1])
    abline(lm(data1950[[whichwine]]~data1950[[whichx]]), col=colpalette[1])
    data195080=subset(dater, year>1949 & year<1980)
    points(data195080[[whichwine]]~data195080[[whichx]], pch=16, col=colpalette[2])
    abline(lm(data195080[[whichwine]]~data195080[[whichx]]), col=colpalette[2])
    data1980=subset(dater, year>1979)
    points(data1980[[whichwine]]~data1980[[whichx]], pch=16, col=colpalette[3])
    abline(lm(data1980[[whichwine]]~data1980[[whichx]]), col=colpalette[3])
  }


##
## plots!
## simple time-series looksees #
par(mfrow=c(2,1))
ylim=c(-1,6)

plot(stars~year, data=subset(allstars4reg, region=="bordeaux" &
    colour.name=="red"), type="l", col="darkred", main="Bordeaux Red & White",
     ylim=ylim)
lines(stars~year, data=subset(allstars4reg, region=="bordeaux" &
    colour.name=="white"), type="l", col="skyblue", ylim=ylim)
plot(stars~year, data=subset(allstars4reg, region=="burgundy" &
    colour.name=="red"), type="l", col="darkred", main="Burgundy Red & White",
     ylim=ylim)
lines(stars~year, data=subset(allstars4reg, region=="burgundy" &
    colour.name=="white"), type="l", col="skyblue", ylim=ylim)
plot(stars~year, data=subset(allstars4reg, region=="alsace"), type="l",  
     col="darkblue", main="Alsace",ylim=ylim)
plot(stars~year, data=subset(allstars4reg, region=="champagne"), type="l",
     col="pink", main="Champagne!", ylim=ylim)

##
## plots!

par(mar = c(5,4,2,1))
par(mfrow=c(2,2)) # GHD
threebreaksplot(bbshapeenv, "red_bordeaux", "Bor.ghd", colz3, "Red Bordeaux")
threebreaksplot(bbshapeenv, "red_burgundy", "Bur.ghd", colz3, "Red Burgundy")
threebreaksplot(bbshapeenv, "white_bordeaux", "Bor.ghd", colz3, "White Bordeaux")
threebreaksplot(bbshapeenv, "white_burgundy", "Bur.ghd", colz3, "White Burgundy")

par(mfrow=c(2,2)) # temp
threebreaksplot(bbshapeenv, "red_bordeaux", "Bor.temp", colz3, "Red Bordeaux")
threebreaksplot(bbshapeenv, "red_burgundy", "Bur.temp", colz3, "Red Burgundy")
threebreaksplot(bbshapeenv, "white_bordeaux", "Bor.temp", colz3, "White Bordeaux")
threebreaksplot(bbshapeenv, "white_burgundy", "Bur.temp", colz3, "White Burgundy")

par(mfrow=c(2,2)) # pdsi
threebreaksplot(bbshapeenv, "red_bordeaux", "Bor.pdsi", colz3, "Red Bordeaux")
threebreaksplot(bbshapeenv, "red_burgundy", "Bur.pdsi", colz3, "Red Burgundy")
threebreaksplot(bbshapeenv, "white_bordeaux", "Bor.pdsi", colz3, "White Bordeaux")
threebreaksplot(bbshapeenv, "white_burgundy", "Bur.pdsi", colz3, "White Burgundy")

par(mfrow=c(1,3))
plot(red_bordeaux~Bor.ghd, data=bbshapeenv)
abline(lm(red_bordeaux~Bor.ghd, data=bbshapeenv))
summary(lm(red_bordeaux~Bor.ghd, data=bbshapeenv))

plot(red_bordeaux~Bor.temp, data=bbshapeenv)
abline(lm(red_bordeaux~Bor.temp, data=bbshapeenv))
summary(lm(red_bordeaux~Bor.temp, data=bbshapeenv))

plot(red_bordeaux~Bor.pdsi, data=bbshapeenv)
abline(lm(red_bordeaux~Bor.pdsi, data=bbshapeenv))
summary(lm(red_bordeaux~Bor.pdsi, data=bbshapeenv))

par(mfrow=c(1,2))
plot(red_bordeaux~Bor.temp, data=subset(bbshapeenv, year<1980 & year>1950))
abline(lm(red_bordeaux~Bor.temp, data=subset(bbshapeenv, year<1980 & year>1950)))
summary(lm(red_bordeaux~Bor.temp, data=subset(bbshapeenv, year<1980 & year>1950)))

plot(red_bordeaux~Bor.temp, data=subset(bbshapeenv, year>1980))
abline(lm(red_bordeaux~Bor.temp, data=subset(bbshapeenv, year>1980)))
summary(lm(red_bordeaux~Bor.temp, data=subset(bbshapeenv, year>1980)))

par(mfrow=c(1,2))
plot(red_bordeaux~Bor.pdsi, data=subset(bbshapeenv, year<1980 & year>1950))
abline(lm(red_bordeaux~Bor.pdsi, data=subset(bbshapeenv, year<1980 & year>1950)))
summary(lm(red_bordeaux~Bor.pdsi, data=subset(bbshapeenv, year<1980 & year>1950)))

plot(red_bordeaux~Bor.pdsi, data=subset(bbshapeenv, year>1980))
abline(lm(red_bordeaux~Bor.pdsi, data=subset(bbshapeenv, year>1980)))
summary(lm(red_bordeaux~Bor.pdsi, data=subset(bbshapeenv, year>1980)))
