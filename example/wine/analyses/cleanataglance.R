### Started 8.December.2014 ###
### By Lizzie ###

## This file takes some from Veg05_Veg06_CompareFile.R (my PhD) #
## to do a QA/QC on the 'at a glance' vintage data from Broadbent's book ##
## It also add in the stars <3 (since at a glance is 3-5 stars) for ...
# several regions for 1900 onward (but not before) ##

## Updated as of 8 April 2015: Nicely cleaned for now ##
## It also includes some work to pull out which years we needed more data for ##
## For that: see also needstars_README.txt ##

# housekeeping #
rm(list=ls()) 
options(stringsAsFactors=FALSE)

library(plyr)
library(reshape)

# setwd("~/Documents/git/teaching/demoSweave/example/wine/analyses")

# get data #
aagsg.raw <- read.csv("data/input/vintages_sg.csv", header=TRUE)
aaghe.raw <- read.csv("data/input/vintages_he.csv", header=TRUE)

# get rid of some random NA data
aaghe <- subset(aaghe.raw, region!="")
aagsg <- subset(aagsg.raw, region!="")

# make things match
aaghe$region <- tolower(aaghe$region) 
aaghe$colour.name <- tolower(aaghe$colour.name)

sort(unique(aaghe$region))
sort(unique(aagsg$region))

# small issues: space in nothern rhone for Sally and...
# For Italy, Harold entered "piedmont italy" and "tuscany italy"
# while Sally wrote each in the notes (as well as some (v) for variable)
aagsg$region[aagsg$region=="northern rhone "] <- "northern rhone"
aagsg[grep("pied", aagsg$notes),]$region <- "piedmont italy"
aagsg[grep("tusc", aagsg$notes),]$region <- "tuscany italy"

# For 'colour.name,' they refer to certain things the same, but using different words
aaghe$colour.name[aaghe$colour.name=="not sure" &
    aaghe$region=="alsace"] <- "unclear"
aaghe$colour.name[aaghe$colour.name=="multiple" &
    aaghe$region=="california"] <- "varies"
aaghe$colour.name[aaghe$colour.name=="white" &
    aaghe$region=="champagne"] <- "champagne"
aaghe$colour.name[aaghe$colour.name=="mulltiple" &
    aaghe$region=="loire"] <- "varies"
aaghe$colour.name[aaghe$colour.name=="multiple" &
    aaghe$region=="madeira"] <- "varies"
aaghe$colour.name[aaghe$colour.name=="unclear" &
    aaghe$region=="port"] <- "red"
aaghe$colour.name[aaghe$colour.name=="white" &
    aaghe$region=="tokaji"] <- "sweet"

## merge! to find errors ##
aag <- merge(aagsg, aaghe, by=c("region", "colour.name", "year", "stars"),
    all.x=TRUE, all.y=TRUE, suffixes=c(".sg", ".he"))


##################################
## manually fixing some errors ##
##################################
## Now fix the errors (see ataglancemerged_cleaningnotes.csv for Jehane's full notes) ##
## Non-matches that do not need fixing:
# 1875 Burgundy: 1875 red burgundy is 5 stars (correct by Sally, sg). Harold (he) forgot to enter.
# 1889 Tokaji sweet: he forgot to enter, sg correct
# 29 rows where Sally skipped a page, no change needed!
# Non-matches that require some deleting (at least one always got it right, so no work other than deleting!)

# 1933 Burgundy white: sg correct, delete 1833 mistyped entry by Harold
aag.cleaned <- aag[!aag$year==1833,] 
# 1876 Burgundy: No information for 1876 red burgundy in book. he probably meant to enter 1875
aag.cleaned <- aag.cleaned[!(aag.cleaned$year==1876 & aag.cleaned$region=="burgundy"),]
    # there's no white in 1876 so this is okay
aag.cleaned <- aag.cleaned[!(aag.cleaned$region=="germany" & aag.cleaned$year==1942
    & aag.cleaned$stars==4),] # errors in book, he entry appears correct
aag.cleaned <- aag.cleaned[!(aag.cleaned$region=="germany" & aag.cleaned$year==1959
    & aag.cleaned$stars==4),] # errors in book, he entry appears correct
aag.cleaned <- aag.cleaned[!(aag.cleaned$region=="northern rhone" & aag.cleaned$year==1900),]
    # he typed 1900 instead of 1990. sg's entry (N. Rhone, red, 1990, 5 stars) is correct.
aag.cleaned <- aag.cleaned[!(aag.cleaned$region=="tuscany italy" & aag.cleaned$year==1952
    & aag.cleaned$stars==5),] # sg typed 1952 instead of 1962. he entry for 1962 (5 stars) correct. 
aag.cleaned <- aag.cleaned[!(aag.cleaned$region=="germany" & aag.cleaned$year==1963
    & aag.cleaned$stars==4),] # he accidentally made two entries (and none for 1964) so remove this one 


#########################
## final little tweaks ##
#########################
# select down to fewer columns and decide what to do about few small data issues
aag.sm <- subset(aag.cleaned, select=c("region", "year", "stars", "colour.name"))
aag.checkagg <- aggregate(aag.sm["stars"], aag.sm[c("region", "colour.name", "year")], FUN=length)
# okay, the 5 duplicates are outside of Burgundy and Bordeaux ...
# ... so ignore for now
unique(aag.checkagg$stars)
subset(aag.checkagg, stars==2)

# aggregate over colour (though why I would want this I am not sure?!)
aag.agg <- aggregate(aag["stars"], aag[c("region", "year")], FUN=mean)

write.csv(aag.agg, "data/output/ataglancemixedcolors.csv", row.names=FALSE)
write.csv(aag, "data/output/ataglancemerged.csv", row.names=FALSE)


###################################
### Figuring out needstarscolor ###
###################################
# The 'at a glance' data are only for stars 3-5
# so any vintages that are NA or 1-2 are missing

# Back in Dec 2014 I thought I would add in 1.5 as quality for other years maybe?
# But need to worry about how short some time series are ...
ddply(aag,~region,summarise, start=min(year), end=max(year), mean=mean(year))

addinyrs <- data.frame(year=rep(c(1900:2000), 14), stars=rep(rep(1.5,
    length(1900:2000)), 14), region=rep(unique(aag$region),
    each=length(1900:2000)))

aag.agg1900on <- merge(aag.agg, addinyrs, by=c("year", "region"), all.y=TRUE)
aag.agg1900on$stars.x[which(is.na(aag.agg1900on$stars.x)==TRUE)] <- 1.5
aag.agg1900on$stars.y <- NULL

aag.agg1900on <- subset(aag.agg1900on, region!="california" & region !="loire" &
    region != "piedmont italy" & region != "southern rhone" & region != "tokaji" &
    region !="tuscany italy")

ddply(aag.agg1900on,~region,summarise, start=min(year), end=max(year), mean=mean(year))

# write.csv(aag.agg1900on, "data/output/ataglanceaddindata.csv", row.names=FALSE)

# Do it for dataframe with color for Bordeaux & Burgundy #
# The below is me figuring out which data to ask Jehane Samaha to get #
# As opposed to adding in data #
burborwcolor <- subset(aag.sm, region=="bordeaux"|
    region=="burgundy")

addinyrscol <- data.frame(year=rep(c(1900:2000), 4), stars=rep(rep(1.5,
    length(1900:2000)), 4), region=rep(rep(unique(burborwcolor$region),
    each=length(1900:2000)), 2), colour.name=rep(c("white", "white", "red", "red"),
    each=length(1900:2000)))

aag.aggwcolor.1900 <- merge(aag.sm, addinyrscol,
    by=c("year", "region", "colour.name"), all.y=TRUE)
aag.aggwcolor.1900$stars.x[which(is.na(aag.aggwcolor.1900$stars.x)==TRUE)] <- 1.5
aag.aggwcolor.1900$stars.y <- NULL
names(aag.aggwcolor.1900)[names(aag.aggwcolor.1900)=="stars.x"] <- "stars"

write.csv(aag.aggwcolor.1900, "data/output/needstarscolor.csv", row.names=FALSE)

burborwcolor <- subset(aag.sm, region=="bordeaux"|region=="burgundy")

addinyrscol <- data.frame(year=rep(c(1900:2000), 4), stars=rep(rep(1.5,
    length(1900:2000)), 4), region=rep(rep(unique(burborwcolor$region),
    each=length(1900:2000)), 2), colour.name=rep(c("white", "white", "red", "red"),
    each=length(1900:2000)))

aag.aggwcolor.1900 <- merge(aag.sm, addinyrscol,
    by=c("year", "region", "colour.name"), all.y=TRUE)
aag.aggwcolor.1900$stars.x[which(is.na(aag.aggwcolor.1900$stars.x)==TRUE)] <- 1.5
aag.aggwcolor.1900$stars.y <- NULL
names(aag.aggwcolor.1900)[names(aag.aggwcolor.1900)=="stars.x"] <- "stars"

write.csv(aag.aggwcolor.1900, "data/output/needstarscolor.csv", row.names=FALSE)


################################################
## Make up complete at a glance for 4 regions ##
################################################

# Data for Bordeuax & Burgundy, complete below:
# The below are data from Jehane after ...

# Below has all data but some of it was pre-full cleaning ...
# so, remove 3-5 stars and merge them back in
aag.bb.all <- read.csv("data/input/needstars_full.csv", header=TRUE)
aag.bb.poor <- subset(aag.bb.all, stars<3|is.na(stars)==TRUE)
burborwcolor.tobind <- data.frame(year=burborwcolor$year, region=burborwcolor$region, 
    colour.name=burborwcolor$colour.name, stars=burborwcolor$stars)
burborwcolor.tobind <- subset(burborwcolor.tobind, year>1899)
    # remember, Jehane only pulled data 1900 onward for stars 0-2
aag.bb <- rbind(aag.bb.poor[,1:4], burborwcolor.tobind)

# data for Alsace/Champange, only the 0-3 stars data (that is data not in at a glance)
aag.acinput <- read.csv("data/input/needstars.csv", header=TRUE)
aag.ac <- subset(aag.acinput, region=="champagne" | region=="alsace")
aag.agg.ac <- subset(aag.agg, region=="champagne" | region=="alsace")

aag.acsm <- subset(aag.ac, select=c("region", "year", "stars"))
aagfull.ac <- rbind(aag.agg.ac, aag.acsm)

aag.acbind <- data.frame(year=aagfull.ac$year, region=aagfull.ac$region, 
    colour.name=NA, stars=aagfull.ac$stars)
aag.bbsm <- aag.bb[,1:4]

allstars4reg <- rbind(aag.bbsm, aag.acbind)

write.csv(allstars4reg, "data/output/allstars4reg.csv", row.names=FALSE)

###########
###########

