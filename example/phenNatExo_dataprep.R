### Started 13 July 2012 (Friday the 13th) ###
### By Lizzie ###

# setwd("/Users/Lizzie/Documents/R/NCEAS/Phenology/NatExotic")
# library(picante)
# options(stringsAsFactors=FALSE)

source("phen_changeovertimefxs.R")
source("phenNatExo_helperfxs_cleanFittertraits.R")

scrubtotrees <- function(data, genuscolumn){
    data[[genuscolumn]][data[[genuscolumn]]=="Triodanus"] <- "Triodanes"
    data[[genuscolumn]][data[[genuscolumn]]=="Amphicarpa"] <- "Amphicarpaea"
    data[[genuscolumn]][data[[genuscolumn]]=="Roripa"] <- "Rorippa"
    data[[genuscolumn]][data[[genuscolumn]]=="Raphiolepis"] <- "Rhaphiolepis"
    data[[genuscolumn]][data[[genuscolumn]]=="Tamus"] <- "Dioscorea"
    # seriously, it's the only Tamus in the ffd file, but double-check the DC full data sometime
    return(data)
  }

##
## get trees
maticpath <- "input/trees/Phylomatic_trees_for_phenology_June2012/"
conctre.mat <- read.tree(paste(maticpath,"Concord_tree_Phylomatic_June2012.txt", sep=""))
dctre.mat <- read.tree(paste(maticpath,"DC_tree_Phylomatic_June2012.txt", sep=""))
fargotre.mat <- read.tree(paste(maticpath,"Fargo_tree_Phylomatic_June2012.txt", sep=""))
fittre.mat <- read.tree(paste(maticpath,"Fitter_tree_Phylomatic_June2012.txt", sep=""))
konztre.mat <- read.tree(paste(maticpath,"Konza_tree_Phylomatic_June2012.txt", sep=""))

maticfpath <- "input/trees/PhylomaticFULL_trees_for_phenology_June2012/"
conctref.mat <- read.tree(paste(maticfpath,"Phylomatic_Concord.txt", sep=""))
dctref.mat <- read.tree(paste(maticfpath,"Phylomatic_WashDC_tree.txt", sep=""))
fargotref.mat <- read.tree(paste(maticfpath,"Phylomatic_Fargo_tree.txt", sep=""))
fittref.mat <- read.tree(paste(maticfpath,"Phylomatic_Fitter_tree.txt", sep=""))
konztref.mat <- read.tree(paste(maticfpath,"Phylomatic_Konza_tree.txt", sep=""))

maticfpathUM <- "input/trees/PhylomaticFULL_trees_for_phenology_June2012/ultrametric/"
conctrefum.mat <- read.tree(paste(maticfpathUM,"Concord_ultra_tree.txt", sep=""))
dctrefum.mat <- read.tree(paste(maticfpathUM,"WashDC_ultra_tree.txt", sep=""))
fargotrefum.mat <- read.tree(paste(maticfpathUM,"Fargo_ultra_tree.txt", sep=""))
fittrefum.mat <- read.tree(paste(maticfpathUM,"Fitter_ultra_tree.txt", sep=""))
konztrefum.mat <- read.tree(paste(maticfpathUM,"Konza_ultra_tree.txt", sep=""))

plawpath <- "input/trees/Phlawd_trees_for_phenology_June2012/"
conctre.pl <- read.tree(paste(plawpath,"Concord_tree_Phlawd_June2012.txt", sep=""))
dctre.pl <- read.tree(paste(plawpath,"DC_tree_Phlawd_June2012.txt", sep=""))
fargotre.pl <- read.tree(paste(plawpath,"fargo_Phlawd_tree_June2012.txt", sep=""))
fittre.pl <- read.tree(paste(plawpath,"Fitter_tree_Phlawd_June2012.txt", sep=""))
konztre.pl <- read.tree(paste(plawpath,"Konza_Phlawd_tree_June2012.txt", sep=""))

plawbladjpath <- "input/trees/Phlawd_wbladj/"
conctre.pl.bl <- read.tree(paste(plawbladjpath,"Concord_ultrametric_tree.txt", sep=""))
dctre.pl.bl <- read.tree(paste(plawbladjpath,"DC_ultrametric_tree.txt", sep=""))
fargotre.pl.bl <- read.tree(paste(plawbladjpath,"Fargo_ultrametric_tree.txt", sep=""))
fittre.pl.bl <- read.tree(paste(plawbladjpath,"Fitter_ultrametric_tree.txt", sep=""))
konztre.pl.bl <- read.tree(paste(plawbladjpath,"Konza_ultrametric_tree.txt", sep=""))

conctre.pl$tip.label <- tolower(conctre.pl$tip.label)
dctre.pl$tip.label <- tolower(dctre.pl$tip.label)
fargotre.pl$tip.label <- tolower(fargotre.pl$tip.label)
fittre.pl$tip.label <- tolower(fittre.pl$tip.label)
konztre.pl$tip.label <- tolower(konztre.pl$tip.label)

conctre.pl.bl$tip.label <- tolower(conctre.pl.bl$tip.label)
dctre.pl.bl$tip.label <- tolower(dctre.pl.bl$tip.label)
fargotre.pl.bl$tip.label <- tolower(fargotre.pl.bl$tip.label)
fittre.pl.bl$tip.label <- tolower(fittre.pl.bl$tip.label)
konztre.pl.bl$tip.label <- tolower(konztre.pl.bl$tip.label)

# throw gymnosperms out of DC trees
dcgymnos <- read.table("input/trees/DCgymnosperms.txt", header=TRUE)
dcgymnos$latbi <- tolower(dcgymnos$latbi)

drop.tip(dctre.mat, which(dctre.mat[["tip.label"]] %in%
    dcgymnos[["latbi"]]))-> dctre.mat 
drop.tip(dctref.mat, which(dctref.mat[["tip.label"]] %in%
    dcgymnos[["latbi"]]))-> dctref.mat 
drop.tip(dctrefum.mat, which(dctrefum.mat[["tip.label"]] %in%
    dcgymnos[["latbi"]]))-> dctrefum.mat 
drop.tip(dctre.pl, which(dctre.pl[["tip.label"]] %in%
    dcgymnos[["latbi"]]))-> dctre.pl 
drop.tip(dctre.pl.bl, which(dctre.pl.bl[["tip.label"]] %in%
    dcgymnos[["latbi"]]))-> dctre.pl.bl 



##
## get the nativity and noxious weeds files
usda.washed.out <- read.csv("input/traits/usda/usdanatexo.csv", header=TRUE)
usda.washed.out <- scrubtotrees(usda.washed.out, "genus")
usda.washed.out$latbi <- tolower(paste(usda.washed.out$genus,
    usda.washed.out$species, sep="_"))

noxstate <- read.csv("input/traits/usda/noxious5sites.csv", header=TRUE)
noxstate <- scrubtotrees(noxstate, "genus")
noxstate$latbi <- tolower(paste(noxstate$genus, noxstate$species, sep="_"))

gentraitstra <- read.csv("input/traits/PhenTraits2012.csv", header=TRUE)
gentraitstra <- scrubtotrees(gentraitstra, "genus")
gentraitstra$latbi <- tolower(paste(gentraitstra$genus, gentraitstra$species, sep="_"))

##
## get the species list to help make exotic lists
spp <- read.csv("input/specieslists/other4_allspp_5yrssplist.csv", header=TRUE)
sppdc <- read.csv("input/specieslists/washdc_allspp_5yrssplist.csv", header=TRUE)
sppdc$latbi <- tolower(paste(sppdc$genus, sppdc$species, sep="_"))

spp$latbi <- tolower(paste(spp$genus, spp$species, sep="_"))
sppconc <- subset(spp, site=="concord")
sppfar <- subset(spp, site=="fargo")
sppfit <- subset(spp, site=="fitter")
sppkonz <- subset(spp, site=="konza")

##
## Make up the exotic list for Concord 
noxconc1 <- merge(sppconc, subset(noxstate, site=="concord"), by="latbi", all.x=FALSE)
noxconc2 <- merge(sppconc, subset(noxstate, site=="concord" & (species=="L." |
    species=="Lour." |species=="Thunb.")), by="genus", all.x=FALSE,
    suffixes=c("","s"))
noxconc <- c(noxconc1$latbi, noxconc2$latbi)

## add in abundance data for Concord
concAbun <- read.delim("input/traits/data_6cat.txt", header=TRUE)
concAbun$name <- tolower(concAbun$name)
names(concAbun)[names(concAbun)=="name"] <- "latbi"
concAbunsm <- subset(concAbun, select=c("latbi", "DeltaAbun_cont"))
# note to self: DeltaAbun_cont is Primack abundance MINUS Hosmer abundance
# so declines in abundance are negative #s and increases in abundance are positive #s

## Get the list of add-ons for native.exotic status for Concord ##
## Most of these showed up as 
concadd <- read.csv("input/traits/natexo_addconc.csv", header=TRUE, na.strings="NA")
concadd <- scrubtotrees(concadd, "genus")
concadd$latbi <- tolower(paste(concadd$genus, concadd$species, sep="_"))

conc.natexo <- rbind(usda.washed.out, concadd)

##
## ND does not have any invasives at the genus only level
noxfar.prep <- merge(sppfar, subset(noxstate, site=="fargo"), by="latbi", all.x=FALSE)
noxfar <- noxfar.prep$latbi
faradd <- read.csv("input/traits/natexo_addfargo.csv", header=TRUE, na.strings="NA")
faradd <- scrubtotrees(faradd, "genus")
faradd$latbi <- tolower(paste(faradd$genus, faradd$species, sep="_"))
far.natexo <- rbind(usda.washed.out, faradd)

##
## DC and invaders
## oh dear, one damn species! (the list is only 8 species long, even though it's 2 states!)
noxdc <- merge(sppdc, subset(noxstate, site=="washdc"), by="latbi", all.x=FALSE)

##
## deal with fitter data
fittertra <- read.csv("input/traits/fromCDavis/Nativity_status_Fitter.csv", header=TRUE)
fittertra_fromHulme <- read.csv("input/traits/New_Phytologist_Variables_fitter.csv", header=TRUE)
## clean fitter from Chuck #
# note that this file has all clean names and does not need to be scrubbed to match trees
fittertra.add <- data.frame(phylomatic_names=c("hordeum_vulgare", "larix_decidua", "taxus_baccata"), Status_Key=c("Non-invasive", "Native", "Native"))
fittertra <- rbind(fittertra, fittertra.add)
names(fittertra)[names(fittertra)=="phylomatic_names"] <- "latbi.unclean"
fittertra$native.exotic <- fittertra$Status_Key
fittertra$native.exotic[fittertra$native.exotic=="Native"] <- "native"
fittertra$native.exotic[fittertra$native.exotic=="Non-invasive"] <- "exotic"
fittertra$native.exotic[fittertra$native.exotic=="Invasive"] <- "exotic"
fittertra$inv <- fittertra$Status_Key
fittertra$inv[fittertra$inv=="Native"] <- "noninv"
fittertra$inv[fittertra$inv=="Non-invasive"] <- "noninv"
fittertra$inv[fittertra$inv=="Invasive"] <- "inv"
fittertra <- cleanfittertra(fittertra)
fittertra$latbi <- tolower(fittertra$latbi.unclean)
## clean up some species spelling ##


## clean fitter from Hulme # change stuff to actual values to prevent confusion #
fittertra_fromHulme$Status <- as.character(fittertra_fromHulme$Status)
fittertra_fromHulme$Status[fittertra_fromHulme$Status=="1"] <- "native"
fittertra_fromHulme$Status[fittertra_fromHulme$Status=="2"] <- "neophyte"
fittertra_fromHulme$Status[fittertra_fromHulme$Status=="3"] <- "archaeohyte"
# compare two methods
fittercompare <- merge(fittertra_fromHulme, fittertra, by.y="latbi.unclean", by.x="Name")
uhoh <- cbind(fittercompare$Status, fittercompare$nat.status) # not so good
# and we only have info for 347 (of 384) species from Hulme

##
## Washington, DC -- everyone hang on for the ride
sppdcdat <- read.csv("input/specieslists/washdc_allspp_5yrsdata.csv", header=TRUE)
washdcadd <- read.csv("input/traits/washdccheckeroo.csv", header=TRUE)
washdcadd$native.exotic[washdcadd$native.exotic=="E"] <- "exotic"
washdcadd$native.exotic[washdcadd$native.exotic=="N"] <- "native"
washdcadd$native.exotic[washdcadd$native.exotic=="NE"] <- "<NA>"
washdcaddsm <- subset(washdcadd, select=c("genus", "species", "native.exotic", "latbi"))
usda.washed.out.dc <- rbind(usda.washed.out, washdcaddsm)
usda.washed.out.dc <- usda.washed.out.dc[-(duplicated(usda.washed.out.dc)),]
sppdc.cul <- aggregate(sppdcdat["doy"], sppdcdat[c("genus", "species","CUL")], FUN=length)
sppdc.cul$doy <- NULL
sppdc.cul <- scrubtotrees(sppdc.cul, "genus")
sppdc.cul$latbi <- tolower(paste(sppdc.cul$genus, sppdc.cul$species, sep="_"))
sppdctraits <- merge(sppdc.cul, usda.washed.out.dc, by=c("genus", "species", "latbi"),
    all.x=TRUE)
# sppdctraits.fornow <- subset(sppdccheck, CUL=="N")
# sppdctraits.fornow$site <- "washdc"
sppdctraits$site <- "washdc"
# rm some duplicates
sppdctraits[duplicated(sppdctraits),]
sppdctraits <- sppdctraits[-1253,]
sppdctraits <- sppdctraits[-993,]
# for now, throw out species with noth N and C culativated, and one spp
# (Osmorhiza) with multiple native/exotic status
sppNC <- c("iris_ensata", "magnolia_virginiana", "muscari_armeniacum",
    "osmorhiza_aristata", "salix_caprea")
sppdctraits <- sppdctraits[which(!sppdctraits$latbi %in% sppNC),]

## get the change over time or split info ##
ffdfld <- read.csv("/Users/Lizzie/Documents/Subversion/phenology/Data/pheno_raw.csv", header=TRUE)
ffd <- subset(ffdfld, event=="ffd")
ffd$date <- as.Date(ffd$date)
ffd$doy <- as.numeric(ffd$doy)
ffd <- scrubtotrees(ffd, "genus")
ffd$latbi <- tolower(paste(ffd$genus, ffd$species, sep="_"))

concchange <- changeoversplit(subset(ffd, site=="concord"), 1904)
fargochange <- changeoversplit(subset(ffd, site=="fargo"), 1962)

washdcchange <- changeovertime(subset(ffd, site=="washdc"))
fitterchange <- changeovertime(subset(ffd, site=="fitter"))
konzachange <- changeovertime(subset(ffd, site=="konza"))


## konza soil moisture sensitivities (added 29 August 2012)
konzsma <- read.csv("input/sensitivities/soilmoisture/konza_kans_usa_ALLSEAS.konza1.norm.mod.csv", header=TRUE) # watershed 1D which is burned annually, see http://www.konza.ksu.edu/knz/pages/research/burnhistory.aspx
                   
konzsmab <- read.csv("input/sensitivities/soilmoisture/konza_kans_usa_ALLSEAS.konzab.norm.mod.csv", header=TRUE) # watershed 20B which is only burned every so often accidentally
##
## f(x)s
## to make the world a tidier place
##

konzatraits <- function(path, sitename, traitfile){
    ## note: no KS state noxious weeds at Konza ##
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath, "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep="")) 
    site.gddpp <- read.csv(paste(sitepath, "gddprec.mod.csv", sep=""))
    site.sma <- konzsma
    site.smab <- konzsmab
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.sma <- scrubtotrees(site.sma, "Genus")
    site.smab <- scrubtotrees(site.smab, "Genus")
    # site.pp <- read.csv(paste(sitepath, "prec.mod.csv", sep=""))
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus, site.gddpp$Species, sep="_"))
    site.sma$latbi <- tolower(paste(site.sma$Genus, site.sma$Species, sep="_"))
    site.smab$latbi <- tolower(paste(site.smab$Genus, site.smab$Species, sep="_"))
    ## native/exotic/invasive codings
    ## note: no KS state noxious weeds at Konza
    sitetra <- subset(traitfile, site==sitename)
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi")
    outputter$mat <- merge(sitetra, site.mat, by="latbi")
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi")
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi")
    outputter$sma <- merge(sitetra, site.sma, by="latbi")
    outputter$smab <- merge(sitetra, site.smab, by="latbi")
    return(outputter)
  }


washdctraits <- function(path, sitename, traitfile){
    ## note: only 1 state (Maryland and Virginia) noxious weeds ##
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.csv", sep=""))
    site.mat <- read.csv(paste(sitepath, "annual.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.csv", sep="")) 
    site.gddpp <- read.csv(paste(sitepath, "gddprec.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    # site.pp <- read.csv(paste(sitepath, "prec.mod.csv", sep=""))
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus, site.gddpp$Species, sep="_"))
    ## native/exotic/invasive codings
    sitetra <- subset(traitfile, site==sitename)
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi")
    outputter$mat <- merge(sitetra, site.mat, by="latbi")
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi")
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi")
    return(outputter)
  }


conctraits <- function(path, sitename, traitfile, noxlist){
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.T-ONLY.csv", sep=""))
    site.mat <- read.csv(paste(sitepath,  "annual.temp.T-ONLY.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.T-ONLY.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi", all.y=TRUE)
    outputter$mat <- merge(sitetra, site.mat, by="latbi", all.y=TRUE)
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi", all.y=TRUE)
    return(outputter)
  }

fargotraits <- function(path, sitename, traitfile, noxlist){
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath,  "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep=""))
    site.gddpp <- read.csv(paste(sitepath, "gddprec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus,
        site.gddpp$Species, sep="_"))
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi", all.y=TRUE)
    outputter$mat <- merge(sitetra, site.mat, by="latbi", all.y=TRUE)
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi", all.y=TRUE)
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi", all.y=TRUE)
    return(outputter)
  }


fittertraits <- function(path, sitename, traitfile){
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath,  "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep=""))
    site.gddpp <- read.csv(paste(sitepath, "gddprec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    # Ugh, additional scrubbing specific to Fitter
    site.mo3$Genus[site.mo3$Genus=="Aster"] <- "symphyotrichum"
    site.mat$Genus[site.mat$Genus=="Aster"] <- "symphyotrichum"
    site.gdd$Genus[site.gdd$Genus=="Aster"] <- "symphyotrichum"
    site.gddpp$Genus[site.gddpp$Genus=="Aster"] <- "symphyotrichum"
    site.mo3$Species[site.mo3$Species=="flos-cuculi"] <- "flos_cuculi"
    site.mat$Species[site.mat$Species=="flos-cuculi"] <- "flos_cuculi"
    site.gdd$Species[site.gdd$Species=="flos-cuculi"] <- "flos_cuculi"
    site.gddpp$Species[site.gddpp$Species=="flos-cuculi"] <- "flos_cuculi"
    site.mo3$Species[site.mo3$Species=="nidus-avis"] <- "nidus_avis"
    site.mat$Species[site.mat$Species=="nidus-avis"] <- "nidus_avis"
    site.gdd$Species[site.gdd$Species=="nidus-avis"] <- "nidus_avis"
    site.gddpp$Species[site.gddpp$Species=="nidus-avis"] <- "nidus_avis"
    site.mo3$Species[site.mo3$Species=="uva-crispa"] <- "uva_crispa"
    site.mat$Species[site.mat$Species=="uva-crispa"] <- "uva_crispa"
    site.gdd$Species[site.gdd$Species=="uva-crispa"] <- "uva_crispa"
    site.gddpp$Species[site.gddpp$Species=="uva-crispa"] <- "uva_crispa"
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus,
        site.gddpp$Species, sep="_"))
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi", all.y=TRUE)
    outputter$mat <- merge(sitetra, site.mat, by="latbi", all.y=TRUE)
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi", all.y=TRUE)
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi", all.y=TRUE)
    return(outputter)
  }

conctraits.doy <- function(traitfile, noxlist, changefile){
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    sitetra <- merge(sitetra, changefile, by="latbi")
    return(sitetra)
  }


fargotraits.doy <- function(traitfile, noxlist, changefile){
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    sitetra <- merge(sitetra, changefile, by="latbi")
    return(sitetra)
  }


fitter.doy <- function(traitfile, changefile){
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    sitetra <- merge(sitetra, changefile, by="latbi")
    changefile$latbi[changefile$latbi=="aster_novi-belgii"] <-
        "symphyotrichum_novi-belgii"
    changefile$latbi[changefile$latbi=="lychnis_flos-cuculi"] <-
        "lychnis_flos_cuculi"
    changefile$latbi[changefile$latbi=="ribes_uva-crispa"] <-
        "ribes_uva_crispa"
    return(sitetra)
  }


washdc.doy <- function(traitfile, changefile){
    sitetra <- traitfile
    sitetra <- merge(sitetra, changefile, by="latbi")
    return(sitetra)
  }


konza.doy <- function(traitfile, changefile, sitename){
    ## native/exotic/invasive codings
    sitetra <- subset(traitfile, site==sitename)
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra <- merge(sitetra, changefile, by="latbi", all.y=TRUE)
    return(sitetra)
  }

## For Fargo, possible outlier problem:
## chenopodium_glaucum has a doy change of 67 and daysperHI of almost -40 ##

fargotraits.nocrazychenopodium <- function(path, sitename, traitfile, noxlist){
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath,  "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep=""))
    site.gddpp <- read.csv(paste(sitepath, "gddprec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus,
       site.gddpp$Species, sep="_"))
    site.mo3 <- subset(site.mo3, latbi != "chenopodium_glaucum")
    site.mat <- subset(site.mat, latbi != "chenopodium_glaucum")
    site.gdd <- subset(site.gdd, latbi != "chenopodium_glaucum")
    site.gddpp <- subset(site.gddpp, latbi != "chenopodium_glaucum")
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi", all.y=TRUE)
    outputter$mat <- merge(sitetra, site.mat, by="latbi", all.y=TRUE)
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi", all.y=TRUE)
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi", all.y=TRUE)
    return(outputter)
  }

fargotraits.doy.nocrazychenopodium <- function(traitfile, noxlist, changefile){
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    sitetra <- merge(sitetra, changefile, by="latbi")
    sitetra <- subset(sitetra, latbi != "chenopodium_glaucum")
    return(sitetra)
  }



##
## Add in precip-only models
##

konzatraits.wprecip <- function(path, sitename, traitfile){
    ## note: no KS state noxious weeds at Konza ##
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath, "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep="")) 
    site.gddpp <- read.csv(paste(sitepath, "gddprec.mod.csv", sep=""))
    site.pp <- read.csv(paste(sitepath, "prec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.pp <- scrubtotrees(site.pp, "Genus")
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus, site.gddpp$Species, sep="_"))
    site.pp$latbi <- tolower(paste(site.pp$Genus, site.pp$Species, sep="_"))
    ## native/exotic/invasive codings
    ## note: no KS state noxious weeds at Konza
    sitetra <- subset(traitfile, site==sitename)
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi")
    outputter$mat <- merge(sitetra, site.mat, by="latbi")
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi")
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi")
    outputter$pp <- merge(sitetra, site.pp, by="latbi")
    return(outputter)
  }


washdctraits.wprecip <- function(path, sitename, traitfile){
    ## note: only 1 state (Maryland and Virginia) noxious weeds ##
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.csv", sep=""))
    site.mat <- read.csv(paste(sitepath, "annual.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.csv", sep="")) 
    site.gddpp <- read.csv(paste(sitepath, "gddprec.csv", sep=""))
    site.pp <- read.csv(paste(sitepath, "prec.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.pp <- scrubtotrees(site.pp, "Genus")
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus, site.gddpp$Species, sep="_"))
    site.pp$latbi <- tolower(paste(site.pp$Genus, site.pp$Species, sep="_"))
    ## native/exotic/invasive codings
    sitetra <- subset(traitfile, site==sitename)
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi")
    outputter$mat <- merge(sitetra, site.mat, by="latbi")
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi")
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi")
    outputter$pp <- merge(sitetra, site.pp, by="latbi")
    return(outputter)
  }


# not doing Concord, we did Temp only due to climate data issues
# could redo with post-1887 data, but ugh

fargotraits.wprecip <- function(path, sitename, traitfile, noxlist){
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath,  "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep=""))
    site.gddpp <- read.csv(paste(sitepath, "gddprec.mod.csv", sep=""))
    site.pp <- read.csv(paste(sitepath, "prec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.pp <- scrubtotrees(site.pp, "Genus")
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus,
        site.gddpp$Species, sep="_"))
    site.pp$latbi <- tolower(paste(site.pp$Genus,
        site.pp$Species, sep="_"))
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi", all.y=TRUE)
    outputter$mat <- merge(sitetra, site.mat, by="latbi", all.y=TRUE)
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi", all.y=TRUE)
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi", all.y=TRUE)
    outputter$pp <- merge(sitetra, site.pp, by="latbi", all.y=TRUE)
    return(outputter)
  }


fittertraits.wprecip <- function(path, sitename, traitfile){
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath,  "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep=""))
    site.gddpp <- read.csv(paste(sitepath, "gddprec.mod.csv", sep=""))
    site.pp <- read.csv(paste(sitepath, "prec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.pp <- scrubtotrees(site.pp, "Genus")
    # Ugh, additional scrubbing specific to Fitter
    site.mo3$Genus[site.mo3$Genus=="Aster"] <- "symphyotrichum"
    site.mat$Genus[site.mat$Genus=="Aster"] <- "symphyotrichum"
    site.gdd$Genus[site.gdd$Genus=="Aster"] <- "symphyotrichum"
    site.gddpp$Genus[site.gddpp$Genus=="Aster"] <- "symphyotrichum"
    site.pp$Genus[site.pp$Genus=="Aster"] <- "symphyotrichum"
    site.mo3$Species[site.mo3$Species=="flos-cuculi"] <- "flos_cuculi"
    site.mat$Species[site.mat$Species=="flos-cuculi"] <- "flos_cuculi"
    site.gdd$Species[site.gdd$Species=="flos-cuculi"] <- "flos_cuculi"
    site.gddpp$Species[site.gddpp$Species=="flos-cuculi"] <- "flos_cuculi"
    site.pp$Species[site.pp$Species=="flos-cuculi"] <- "flos_cuculi"
    site.mo3$Species[site.mo3$Species=="nidus-avis"] <- "nidus_avis"
    site.mat$Species[site.mat$Species=="nidus-avis"] <- "nidus_avis"
    site.gdd$Species[site.gdd$Species=="nidus-avis"] <- "nidus_avis"
    site.gddpp$Species[site.gddpp$Species=="nidus-avis"] <- "nidus_avis"
    site.pp$Species[site.pp$Species=="nidus-avis"] <- "nidus_avis"
    site.mo3$Species[site.mo3$Species=="uva-crispa"] <- "uva_crispa"
    site.mat$Species[site.mat$Species=="uva-crispa"] <- "uva_crispa"
    site.gdd$Species[site.gdd$Species=="uva-crispa"] <- "uva_crispa"
    site.gddpp$Species[site.gddpp$Species=="uva-crispa"] <- "uva_crispa"
    site.pp$Species[site.pp$Species=="uva-crispa"] <- "uva_crispa"
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus,
        site.gddpp$Species, sep="_"))
    site.pp$latbi <- tolower(paste(site.pp$Genus, site.pp$Species, sep="_"))
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi", all.y=TRUE)
    outputter$mat <- merge(sitetra, site.mat, by="latbi", all.y=TRUE)
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi", all.y=TRUE)
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi", all.y=TRUE)
    outputter$pp <- merge(sitetra, site.pp, by="latbi", all.y=TRUE)
    return(outputter)
  }


## For Fargo, possible outlier problem:
## chenopodium_glaucum has a doy change of 67 and daysperHI of almost -40 ##

fargotraits.wprecip.nocrazychenopodium <- function(path, sitename, traitfile, noxlist){
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath,  "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep=""))
    site.gddpp <- read.csv(paste(sitepath, "gddprec.mod.csv", sep=""))
    site.pp <- read.csv(paste(sitepath, "prec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.pp <- scrubtotrees(site.pp, "Genus")
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus,
       site.gddpp$Species, sep="_"))
    site.pp$latbi <- tolower(paste(site.pp$Genus, site.pp$Species, sep="_"))      
    site.mo3 <- subset(site.mo3, latbi != "chenopodium_glaucum")
    site.mat <- subset(site.mat, latbi != "chenopodium_glaucum")
    site.gdd <- subset(site.gdd, latbi != "chenopodium_glaucum")
    site.gddpp <- subset(site.gddpp, latbi != "chenopodium_glaucum")
    site.pp <- subset(site.pp, latbi != "chenopodium_glaucum")
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi", all.y=TRUE)
    outputter$mat <- merge(sitetra, site.mat, by="latbi", all.y=TRUE)
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi", all.y=TRUE)
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi", all.y=TRUE)
    outputter$pp <- merge(sitetra, site.pp, by="latbi", all.y=TRUE)
    return(outputter)
  }


##
## new precip models (as of mid-August 2012, these ended up rather useless)
##

konzatraits.fixprec <- function(path, sitename, traitfile){
    ## note: no KS state noxious weeds at Konza ##
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath, "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep="")) 
    site.gddpp <- read.csv(paste(sitepath, "fixgddprec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    # site.pp <- read.csv(paste(sitepath, "prec.mod.csv", sep=""))
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus, site.gddpp$Species, sep="_"))
    ## native/exotic/invasive codings
    ## note: no KS state noxious weeds at Konza
    sitetra <- subset(traitfile, site==sitename)
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi")
    outputter$mat <- merge(sitetra, site.mat, by="latbi")
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi")
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi")
    return(outputter)
  }


fargotraits.nocrazychenopodium.fixprec <- function(path, sitename, traitfile, noxlist){
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath,  "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep=""))
    site.gddpp <- read.csv(paste(sitepath, "fixgddprec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus,
       site.gddpp$Species, sep="_"))
    site.mo3 <- subset(site.mo3, latbi != "chenopodium_glaucum")
    site.mat <- subset(site.mat, latbi != "chenopodium_glaucum")
    site.gdd <- subset(site.gdd, latbi != "chenopodium_glaucum")
    site.gddpp <- subset(site.gddpp, latbi != "chenopodium_glaucum")
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi", all.y=TRUE)
    outputter$mat <- merge(sitetra, site.mat, by="latbi", all.y=TRUE)
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi", all.y=TRUE)
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi", all.y=TRUE)
    return(outputter)
  }


fargotraits.fixprec <- function(path, sitename, traitfile, noxlist){
    outputter <- list()
    sitepath <- path
    site.mo3 <- read.csv(paste(sitepath, "3month.temp.csv", sep=""))
    site.mat <- read.csv(paste(sitepath,  "annual.temp.csv", sep="")) 
    site.gdd <- read.csv(paste(sitepath, "gdd.mod.csv", sep=""))
    site.gddpp <- read.csv(paste(sitepath, "fixgddprec.mod.csv", sep=""))
    site.mo3 <- scrubtotrees(site.mo3, "Genus")
    site.mat <- scrubtotrees(site.mat, "Genus")
    site.gdd <- scrubtotrees(site.gdd, "Genus")
    site.gddpp <- scrubtotrees(site.gddpp, "Genus")
    site.mo3$latbi <- tolower(paste(site.mo3$Genus, site.mo3$Species, sep="_"))
    site.mat$latbi <- tolower(paste(site.mat$Genus, site.mat$Species, sep="_"))
    site.gdd$latbi <- tolower(paste(site.gdd$Genus, site.gdd$Species, sep="_"))
    site.gddpp$latbi <- tolower(paste(site.gddpp$Genus,
        site.gddpp$Species, sep="_"))
    ## native/exotic/invasive codings
    sitetra <- traitfile
    sitetra$inv <- "noninv"
    sitetra$inv[which(sitetra$latbi %in% noxlist)] <- "inv"
    sitetra$full.code <- paste(sitetra$native.exotic, sitetra$inv, sep=".")
    sitetra$full.code[sitetra$full.code=="NA.noninv"] <- NA
    sitetra$full.code[sitetra$full.code=="NA.inv"] <- NA
    sitetra$nat.exoinv <- sitetra$full.code
    sitetra$nat.exoinv[sitetra$native.exotic=="native"] <- "native"
    sitetra$native.exotic12 <- sitetra$native.exotic
    sitetra$native.exotic12[sitetra$native.exotic12=="native"] <- 1
    sitetra$native.exotic12[sitetra$native.exotic12=="exotic"] <- 2
    sitetra$native.exotic12 <- as.numeric(sitetra$native.exotic12)
    sitetra$inv.non12 <- sitetra$inv
    sitetra$inv.non12[sitetra$inv.non12=="inv"] <- 1
    sitetra$inv.non12[sitetra$inv.non12=="noninv"] <- 2
    sitetra$inv.non12 <- as.numeric(sitetra$inv.non12)
    sitetra$nat.exo.inv123 <- sitetra$full.code
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="native.noninv"] <- 1
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.noninv"] <- 2
    sitetra$nat.exo.inv123[sitetra$nat.exo.inv123=="exotic.inv"] <- 3
    sitetra$nat.exo.inv123 <- as.numeric(sitetra$nat.exo.inv123)
    outputter$mo3 <- merge(sitetra, site.mo3, by="latbi", all.y=TRUE)
    outputter$mat <- merge(sitetra, site.mat, by="latbi", all.y=TRUE)
    outputter$gdd <- merge(sitetra, site.gdd, by="latbi", all.y=TRUE)
    outputter$gddpp <- merge(sitetra, site.gddpp, by="latbi", all.y=TRUE)
    return(outputter)
  }
