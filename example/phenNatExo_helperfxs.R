## Some brand new shiny f(x)s showing off my new understanding of list ##
## Some extracted from phenNatExo_2012.R main analysis file in July 2012 ##
## By Lizzie ##
## Summah 2012 ##
## In Vancouver, which is never really summer ##

## f(x)s for data slimming and organizing ##
# latbiname in form "latbi.columnname"
# selectcols in form c("column1", "column2")
pglm.pval <- function(pglm.mod){
    pvalue <- pf((((pglm.mod$NSSQ - pglm.mod$RSSQ) / pglm.mod$RMS) / (pglm.mod$k - 1)),
    pglm.mod$k - 1, pglm.mod$n - pglm.mod$k,  ncp=0, lower.tail = FALSE, log.p = FALSE)
    return(pvalue)
  }

matchtreedat <- function(treedat, traitdat, latbiname, selectcols){
drop.tip(treedat, which(!(treedat[["tip.label"]] %in% traitdat[[latbiname]])) )->weetree
tosswhat <- which(!(traitdat[[latbiname]] %in% treedat[["tip.label"]]))
if (length(which(!(treedat[["tip.label"]] %in% traitdat[[latbiname]])))>0)
  print ("some species don't match tree to species, the tip labels of these are:")
print(treedat[["tip.label"]][which(!(treedat[["tip.label"]] %in% traitdat[[latbiname]]))])
if (length(tosswhat)>0)
datawee <- traitdat[-tosswhat,]
else
datawee <- traitdat
row.names(datawee) <- datawee[[latbiname]]
dataweesorted <- datawee[weetree$tip.label,]
dataweesortedout <- subset(dataweesorted, select=selectcols)
return(dataweesortedout)
}


matchtreephy <- function(treedat, traitdat, latbiname){
drop.tip(treedat, which(!(treedat[["tip.label"]] %in% traitdat[[latbiname]])) )->weetree
weetree <- multi2di(weetree)
return(weetree)
}

treephy.phlawdsize <- function(matictree, phlawdtree){
drop.tip(matictree, which(!(matictree[["tip.label"]] %in% phlawdtree[["tip.label"]])) )->weetree
weetree <- multi2di(weetree)
return(weetree)
}

# spp in data, not in tree
grabdatnomatches <- function(treedat, traitdat, latbi){
tosswhat <- which(!(traitdat[[latbi]] %in% treedat[["tip.label"]]))
tosses <- traitdat[[latbi]][which(!(traitdat[[latbi]] %in% treedat[["tip.label"]]))]
return(tosses)
}

# spp in tree, not in data
grabphynomatches <- function(treedat, traitdat, latbi){
tosswhat <- which(!(traitdat[[latbi]] %in% treedat[["tip.label"]]))
tosses <- treedat[["tip.label"]][which(!(treedat[["tip.label"]] %in% traitdat[[latbi]]))]
return(tosses)
}


# define helper function to replace all occurrences of
# consecutive spaces with a single space, and remove trailing spaces
fix.spaces <- function(x) {
    x <- gsub(" +", " ", x)
    sub(" +$", "", x)
}

##
## zero out non-significant sensitivies
##

senszero <- function(goofile, pvalcutoff){
    goofile$mo3["daysperC"][goofile$mo3["p.val"]>pvalcutoff] <- 0
    goofile$mat["daysperC"][goofile$mat["p.val"]>pvalcutoff] <- 0
    goofile$gdd["daysperHI"][goofile$gdd["p.val"]>pvalcutoff] <- 0
    goofile$gddpp["daysperHI"][goofile$gddpp["p.valGDDterm"]>pvalcutoff] <- 0
    goofile$gddpp["daysperMI"][goofile$gddpp["p.valPRECSUMterm"]>pvalcutoff] <- 0
    goofile$gddpp["daysperINTER"][goofile$gddpp["p.valINTER"]>pvalcutoff] <- 0
    return(goofile)
  }

senszero.conc <- function(goofile, pvalcutoff){
    goofile$mo3["daysperC"][goofile$mo3["p.val"]>pvalcutoff] <- 0
    goofile$mat["daysperC"][goofile$mat["p.val"]>pvalcutoff] <- 0
    goofile$gdd["daysperHI"][goofile$gdd["p.val"]>pvalcutoff] <- 0
    return(goofile)
  }

senszero.wprecip <- function(goofile, pvalcutoff){
    goofile$mo3["daysperC"][goofile$mo3["p.val"]>pvalcutoff] <- 0
    goofile$mat["daysperC"][goofile$mat["p.val"]>pvalcutoff] <- 0
    goofile$gdd["daysperHI"][goofile$gdd["p.val"]>pvalcutoff] <- 0
    goofile$gddpp["daysperHI"][goofile$gddpp["p.valGDDterm"]>pvalcutoff] <- 0
    goofile$gddpp["daysperMI"][goofile$gddpp["p.valPRECSUMterm"]>pvalcutoff] <- 0
    goofile$gddpp["daysperINTER"][goofile$gddpp["p.valINTER"]>pvalcutoff] <- 0
    goofile$pp["daysperMI"][goofile$pp["p.valPRECSUMterm"]>pvalcutoff] <- 0
    return(goofile)
  }

##
##

konzagoo <- function(sitefile, doyfile, treehere){
    goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    sma.dat <- subset(as.data.frame(sitefile[["sma"]]),is.na(native.exotic)==FALSE)
    sma.dat <- subset(sma.dat, native.exotic!="native/exotic" & native.exotic!="")
    smab.dat <- subset(as.data.frame(sitefile[["smab"]]),is.na(native.exotic)==FALSE)
    smab.dat <- subset(smab.dat, native.exotic!="native/exotic" & native.exotic!="")
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperC"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",  "daysperC"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    goo$smamod <- matchtreedat(treehere, sma.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",  "daysperdev", "Month"))
    goo$tree.sma <- matchtreephy(treehere, sma.dat,"latbi")
    goo$smabmod <- matchtreedat(treehere, smab.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",  "daysperdev", "Month"))
    goo$tree.smab <- matchtreephy(treehere, smab.dat,"latbi")
    # goo$doych <- matchtreedat(treehere, doyzer,"latbi",
    #    c("native.exotic", "changeovertime"))
    # goo$treedoy <- matchtreephy(treehere, goo$doych,"latbi")
    # run the phylo signal stuff
    mPgdd <- subset(goo$gdd , select=c("meanFFD", "daysperHI"))
    mPgddpp <- subset(goo$gddpp , select=c("daysperMI", "daysperHI","daysperINTER"))
    mPmo3 <- subset(goo$mo3, select=c("meanFFD", "daysperC"))
    mPmat <- subset(goo$matmod, select=c("daysperC", "native.exotic12"))
    # mPdoy<- subset(goo$doych, select=c("changeovertime"))
    MPS.gdd <- multiPhylosignal(mPgdd, goo$tree.gdd, reps=mPreps)
    MPS.gddpp <- multiPhylosignal(mPgddpp, goo$tree.gddpp, reps=mPreps)
    MPS.mo3 <- multiPhylosignal(mPmo3, goo$tree.mo3, reps=mPreps)
    MPS.mat<- multiPhylosignal(mPmat, goo$tree.mat, reps=mPreps)
    # MPS.doy<- multiPhylosignal(mPdoy, goo$treedoy, reps=mPreps)
    goo$MPS.gdd <- MPS.gdd
    goo$MPS.gddpp <- MPS.gddpp
    goo$MPS.mo3 <- MPS.mo3
    goo$MPS.mat <- MPS.mat
    return(goo)
}

mostsitesgoo <- function(sitefile, doyfile, doychangetype, treehere){
    goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
         c("native.exotic", "full.code", "meanFFD", "daysperHI"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "full.code", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
         c("native.exotic", "full.code", "meanFFD", "daysperC"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "full.code", "meanFFD", "native.exotic12",
        "inv.non12", "daysperC"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    goo$doymod <- matchtreedat(treehere, doyzer,"latbi",
        c("native.exotic", "full.code", doychangetype))
    goo$treedoy <- matchtreephy(treehere, doyzer,"latbi")
    # run the phylo signal stuff
    mPgdd <- subset(goo$gdd , select=c("meanFFD", "daysperHI"))
    mPgddpp <- subset(goo$gddpp , select=c("daysperMI", "daysperHI","daysperINTER"))
    mPmo3 <- subset(goo$mo3, select=c("meanFFD", "daysperC"))
    mPmat <- subset(goo$matmod, select=c("daysperC", "native.exotic12",  "inv.non12"))
    mPdoy<- subset(goo$doymod, select=c(doychangetype))
    MPS.gdd <- multiPhylosignal(mPgdd, goo$tree.gdd, reps=mPreps)
    MPS.gddpp <- multiPhylosignal(mPgddpp, goo$tree.gddpp, reps=mPreps)
    MPS.mo3 <- multiPhylosignal(mPmo3, goo$tree.mo3, reps=mPreps)
    MPS.mat<- multiPhylosignal(mPmat, goo$tree.mat, reps=mPreps)
    MPS.doy<- multiPhylosignal(mPdoy, goo$treedoy, reps=mPreps)
    goo$MPS.gdd <- MPS.gdd
    goo$MPS.gddpp <- MPS.gddpp
    goo$MPS.mo3 <- MPS.mo3
    goo$MPS.mat <- MPS.mat
    goo$MPS.doy <- MPS.doy
    return(goo)
}


dcgoo <- function(sitefile, doyfile, doychangetype, treehere){
    goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
         c("native.exotic", "meanFFD", "daysperHI"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
         c("native.exotic", "meanFFD", "daysperC"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",
        "daysperC"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    goo$doymod <- matchtreedat(treehere, doyzer,"latbi",
        c("native.exotic", doychangetype))
    goo$treedoy <- matchtreephy(treehere, doyzer,"latbi")
    # run the phylo signal stuff
    mPgdd <- subset(goo$gdd , select=c("meanFFD", "daysperHI"))
    mPgddpp <- subset(goo$gddpp , select=c("daysperMI", "daysperHI","daysperINTER"))
    mPmo3 <- subset(goo$mo3, select=c("meanFFD", "daysperC"))
    mPmat <- subset(goo$matmod, select=c("daysperC", "native.exotic12"))
    mPdoy<- subset(goo$doymod, select=c(doychangetype))
    MPS.gdd <- multiPhylosignal(mPgdd, goo$tree.gdd, reps=mPreps)
    MPS.gddpp <- multiPhylosignal(mPgddpp, goo$tree.gddpp, reps=mPreps)
    MPS.mo3 <- multiPhylosignal(mPmo3, goo$tree.mo3, reps=mPreps)
    MPS.mat<- multiPhylosignal(mPmat, goo$tree.mat, reps=mPreps)
    MPS.doy<- multiPhylosignal(mPdoy, goo$treedoy, reps=mPreps)
    goo$MPS.gdd <- MPS.gdd
    goo$MPS.gddpp <- MPS.gddpp
    goo$MPS.mo3 <- MPS.mo3
    goo$MPS.mat <- MPS.mat
    goo$MPS.doy <- MPS.doy
    return(goo)
}



massivetable <- function(goofile){
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- compar.gee(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        family="gaussian", phy= goofile$tree.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- compar.gee(daysperHI~as.factor(native.exotic), data=goofile$gdd, 
        family="gaussian", phy=goofile$tree.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- compar.gee(daysperHI~as.factor(native.exotic), data=goofile$gddpp, 
        family="gaussian", phy= goofile$tree.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- compar.gee(daysperMI~as.factor(native.exotic), data=goofile$gddpp, 
        family="gaussian", phy=goofile$tree.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- compar.gee(daysperINTER~as.factor(native.exotic), 
        data=goofile$gddpp, family="gaussian", phy=goofile$tree.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- compar.gee(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        family="gaussian", phy=goofile$tree.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- compar.gee(daysperC~as.factor(native.exotic), data=goofile$matmod, 
        family="gaussian", phy=goofile$tree.mat)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$coef[1], ffd.gee$coef[1]+ffd.gee$coef[2],
        foo(ffd.gee)[8])
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
         summary(gdd.lm)$coef[8], gdd.gee$coef[1], gdd.gee$coef[1]+gdd.gee$coef[2],
         foo(gdd.gee)[8])
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
         summary(gddpp.lm)$coef[8], gddpp.gee$coef[1],
         gddpp.gee$coef[1]+gddpp.gee$coef[2], foo(gddpp.gee)[8])
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
         summary(gddpp.p.lm)$coef[8], gddpp.p.gee$coef[1],
         gddpp.p.gee$coef[1]+gddpp.p.gee$coef[2], foo(gddpp.p.gee)[8])
     gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
         summary(gddpp.i.lm)$coef[8], gddpp.i.gee$coef[1],
         gddpp.i.gee$coef[1]+gddpp.i.gee$coef[2], foo(gddpp.i.gee)[8])
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         summary(mo3.lm)$coef[8], mo3.gee$coef[1], mo3.gee$coef[1]+mo3.gee$coef[2],
         foo(mo3.gee)[8])
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
          summary(mat.lm)$coef[8], mat.gee$coef[1], mat.gee$coef[1]+mat.gee$coef[2],
          foo(mat.gee)[8])
    daterrows <- rbind(meanffd, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c( "GLM: exo", "GLM: native", "GLM: p",
        "GEE: exotic", "GEE: native", "GEE: p")
return(memodels)
}

moremassivetable <- function(goofile){
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- compar.gee(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        family="gaussian", phy= goofile$tree.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- compar.gee(daysperHI~as.factor(native.exotic), data=goofile$gdd, 
        family="gaussian", phy=goofile$tree.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- compar.gee(daysperHI~as.factor(native.exotic), data=goofile$gddpp, 
        family="gaussian", phy= goofile$tree.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- compar.gee(daysperMI~as.factor(native.exotic), data=goofile$gddpp, 
        family="gaussian", phy=goofile$tree.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- compar.gee(daysperINTER~as.factor(native.exotic), 
        data=goofile$gddpp, family="gaussian", phy=goofile$tree.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- compar.gee(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        family="gaussian", phy=goofile$tree.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- compar.gee(daysperC~as.factor(native.exotic), data=goofile$matmod, 
        family="gaussian", phy=goofile$tree.mat)
    doy.lm <- lm(doychange~as.factor(native.exotic), data=goofile$doymod)
    doy.gee <- compar.gee(doychange~as.factor(native.exotic), data=goofile$doymod,
          family="gaussian", phy=goofile$treedoy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$coef[1], ffd.gee$coef[1]+ffd.gee$coef[2],
        foo(ffd.gee)[8])
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8], gdd.gee$coef[1], gdd.gee$coef[1]+gdd.gee$coef[2],
        foo(gdd.gee)[8])
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8], gddpp.gee$coef[1],
        gddpp.gee$coef[1]+gddpp.gee$coef[2], foo(gddpp.gee)[8])
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8], gddpp.p.gee$coef[1],
        gddpp.p.gee$coef[1]+gddpp.p.gee$coef[2], foo(gddpp.p.gee)[8])
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8], gddpp.i.gee$coef[1],
        gddpp.i.gee$coef[1]+gddpp.i.gee$coef[2], foo(gddpp.i.gee)[8])
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8], mo3.gee$coef[1], mo3.gee$coef[1]+mo3.gee$coef[2],
        foo(mo3.gee)[8])
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8], mat.gee$coef[1], mat.gee$coef[1]+mat.gee$coef[2],
        foo(mat.gee)[8])
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        summary(doy.lm)$coef[8], doy.gee$coef[1], doy.gee$coef[1]+doy.gee$coef[2],
        foo(doy.gee)[8])
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c( "GLM: exotic", "GLM: native", "GLM: p",
        "GEE: exotic", "GEE: native", "GEE: p")
return(memodels)
}

moremassivetable.fc <- function(goofile){
    ffd.lm <- lm(meanFFD~as.factor(full.code), data=goofile$gdd)
    ffd.gee <- compar.gee(meanFFD~as.factor(full.code), data=goofile$gdd,
        family="gaussian", phy= goofile$tree.gdd)
    gdd.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gdd)
    gdd.gee <- compar.gee(daysperHI~as.factor(full.code), data=goofile$gdd, 
        family="gaussian", phy=goofile$tree.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gddpp)
    gddpp.gee <- compar.gee(daysperHI~as.factor(full.code), data=goofile$gddpp, 
        family="gaussian", phy= goofile$tree.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(full.code), data=goofile$gddpp)
    gddpp.p.gee <- compar.gee(daysperMI~as.factor(full.code), data=goofile$gddpp, 
        family="gaussian", phy=goofile$tree.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(full.code), data=goofile$gddpp)
    gddpp.i.gee <- compar.gee(daysperINTER~as.factor(full.code), 
        data=goofile$gddpp, family="gaussian", phy=goofile$tree.gddpp)
    mo3.lm <- lm(daysperC~as.factor(full.code), data=goofile$mo3)
    mo3.gee <- compar.gee(daysperC~as.factor(full.code), data=goofile$mo3, 
        family="gaussian", phy=goofile$tree.mo3)
    mat.lm <- lm(daysperC~as.factor(full.code), data=goofile$matmod)
    mat.gee <- compar.gee(daysperC~as.factor(full.code), data=goofile$matmod, 
        family="gaussian", phy=goofile$tree.mat)
    doy.lm <- lm(doychange~as.factor(full.code), data=goofile$doymod)
    doy.gee <- compar.gee(doychange~as.factor(full.code), data=goofile$doymod,
        family="gaussian", phy=goofile$treedoy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        ffd.lm$coef[1]+ffd.lm$coef[3], anova(ffd.lm)[[5]][1], ffd.gee$coef[1],
        ffd.gee$coef[1]+ffd.gee$coef[2], ffd.gee$coef[1]+ffd.gee$coef[3],
        foo(ffd.gee)[11])
    gddmod <-  c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        gdd.lm$coef[1]+gdd.lm$coef[3], anova(gdd.lm)[[5]][1], gdd.gee$coef[1],
        gdd.gee$coef[1]+gdd.gee$coef[2], gdd.gee$coef[1]+gdd.gee$coef[3],
        foo(gdd.gee)[11])
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        gddpp.lm$coef[1]+gddpp.lm$coef[3], anova(gddpp.lm)[[5]][1], gddpp.gee$coef[1],
        gddpp.gee$coef[1]+gddpp.gee$coef[2], gddpp.gee$coef[1]+gddpp.gee$coef[3],
        foo(gddpp.gee)[11])
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        gddpp.p.lm$coef[1]+gddpp.p.lm$coef[3], anova(gddpp.p.lm)[[5]][1],
        gddpp.p.gee$coef[1],gddpp.p.gee$coef[1]+gddpp.p.gee$coef[2],
        gddpp.p.gee$coef[1]+gddpp.p.gee$coef[3],
        foo(gddpp.p.gee)[11])
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        gddpp.i.lm$coef[1]+gddpp.i.lm$coef[3], anova(gddpp.i.lm)[[5]][1],
        gddpp.i.gee$coef[1],gddpp.i.gee$coef[1]+gddpp.i.gee$coef[2],
        gddpp.i.gee$coef[1]+gddpp.i.gee$coef[3],
        foo(gddpp.i.gee)[11])
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         mo3.lm$coef[1]+mo3.lm$coef[3], anova(mo3.lm)[[5]][1], mo3.gee$coef[1],
         mo3.gee$coef[1]+mo3.gee$coef[2], mo3.gee$coef[1]+mo3.gee$coef[3],
         foo(mo3.gee)[11])
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        mat.lm$coef[1]+mat.lm$coef[3],  anova(mat.lm)[[5]][1], mat.gee$coef[1],
        mat.gee$coef[1]+mat.gee$coef[2], mat.gee$coef[1]+mat.gee$coef[3],
        foo(mat.gee)[11])
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        doy.lm$coef[1]+doy.lm$coef[3],  anova(doy.lm)[[5]][1], doy.gee$coef[1],
        doy.gee$coef[1]+doy.gee$coef[2], doy.gee$coef[1]+doy.gee$coef[3],
        foo(doy.gee)[11])
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("lm:ex-inv", "ex-non","nat-non", "p", "ee:ex-inv",
         "ex-non", "nat-non", "p")
return(memodels)
}



## I don't want to talk about how ugly this is ##
## I really don't ##
## I don't care ##

moremassivetable.overtime <- function(goofile, doychange){
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- compar.gee(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        family="gaussian", phy= goofile$tree.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- compar.gee(daysperHI~as.factor(native.exotic), data=goofile$gdd, 
        family="gaussian", phy=goofile$tree.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- compar.gee(daysperHI~as.factor(native.exotic), data=goofile$gddpp, 
        family="gaussian", phy= goofile$tree.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- compar.gee(daysperMI~as.factor(native.exotic), data=goofile$gddpp, 
        family="gaussian", phy=goofile$tree.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- compar.gee(daysperINTER~as.factor(native.exotic), 
        data=goofile$gddpp, family="gaussian", phy=goofile$tree.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- compar.gee(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        family="gaussian", phy=goofile$tree.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- compar.gee(daysperC~as.factor(native.exotic), data=goofile$matmod, 
        family="gaussian", phy=goofile$tree.mat)
    doy.lm <- lm(changeovertime~as.factor(native.exotic), data=goofile$doymod)
    doy.gee <- compar.gee(changeovertime~as.factor(native.exotic), data=goofile$doymod,
          family="gaussian", phy=goofile$treedoy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$coef[1], ffd.gee$coef[1]+ffd.gee$coef[2],
        foo(ffd.gee)[8])
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8], gdd.gee$coef[1], gdd.gee$coef[1]+gdd.gee$coef[2],
        foo(gdd.gee)[8])
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8], gddpp.gee$coef[1],
        gddpp.gee$coef[1]+gddpp.gee$coef[2], foo(gddpp.gee)[8])
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8], gddpp.p.gee$coef[1],
        gddpp.p.gee$coef[1]+gddpp.p.gee$coef[2], foo(gddpp.p.gee)[8])
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8], gddpp.i.gee$coef[1],
        gddpp.i.gee$coef[1]+gddpp.i.gee$coef[2], foo(gddpp.i.gee)[8])
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8], mo3.gee$coef[1], mo3.gee$coef[1]+mo3.gee$coef[2],
        foo(mo3.gee)[8])
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8], mat.gee$coef[1], mat.gee$coef[1]+mat.gee$coef[2],
        foo(mat.gee)[8])
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        summary(doy.lm)$coef[8], doy.gee$coef[1], doy.gee$coef[1]+doy.gee$coef[2],
        foo(doy.gee)[8])
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c( "GLM: exotic", "GLM: native", "GLM: p",
        "GEE: exotic", "GEE: native", "GEE: p")
return(memodels)
}


moremassivetable.fc.overtime <- function(goofile){
    ffd.lm <- lm(meanFFD~as.factor(full.code), data=goofile$gdd)
    ffd.gee <- compar.gee(meanFFD~as.factor(full.code), data=goofile$gdd,
        family="gaussian", phy=goofile$tree.gdd)
    gdd.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gdd)
    gdd.gee <- compar.gee(daysperHI~as.factor(full.code), data=goofile$gdd, 
        family="gaussian", phy=goofile$tree.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gddpp)
    gddpp.gee <- compar.gee(daysperHI~as.factor(full.code), data=goofile$gddpp, 
        family="gaussian", phy= goofile$tree.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(full.code), data=goofile$gddpp)
    gddpp.p.gee <- compar.gee(daysperMI~as.factor(full.code), data=goofile$gddpp, 
        family="gaussian", phy=goofile$tree.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(full.code), data=goofile$gddpp)
    gddpp.i.gee <- compar.gee(daysperINTER~as.factor(full.code), 
        data=goofile$gddpp, family="gaussian", phy=goofile$tree.gddpp)
    mo3.lm <- lm(daysperC~as.factor(full.code), data=goofile$mo3)
    mo3.gee <- compar.gee(daysperC~as.factor(full.code), data=goofile$mo3, 
        family="gaussian", phy=goofile$tree.mo3)
    mat.lm <- lm(daysperC~as.factor(full.code), data=goofile$matmod)
    mat.gee <- compar.gee(daysperC~as.factor(full.code), data=goofile$matmod, 
        family="gaussian", phy=goofile$tree.mat)
    doy.lm <- lm(changeovertime~as.factor(full.code), data=goofile$doymod)
    doy.gee <- compar.gee(changeovertime~as.factor(full.code), data=goofile$doymod,
        family="gaussian", phy=goofile$treedoy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        ffd.lm$coef[1]+ffd.lm$coef[3], anova(ffd.lm)[[5]][1], ffd.gee$coef[1],
        ffd.gee$coef[1]+ffd.gee$coef[2], ffd.gee$coef[1]+ffd.gee$coef[3],
        foo(ffd.gee)[11])
    gddmod <-  c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        gdd.lm$coef[1]+gdd.lm$coef[3], anova(gdd.lm)[[5]][1], gdd.gee$coef[1],
        gdd.gee$coef[1]+gdd.gee$coef[2], gdd.gee$coef[1]+gdd.gee$coef[3],
        foo(gdd.gee)[11])
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        gddpp.lm$coef[1]+gddpp.lm$coef[3], anova(gddpp.lm)[[5]][1], gddpp.gee$coef[1],
        gddpp.gee$coef[1]+gddpp.gee$coef[2], gddpp.gee$coef[1]+gddpp.gee$coef[3],
        foo(gddpp.gee)[11])
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        gddpp.p.lm$coef[1]+gddpp.p.lm$coef[3], anova(gddpp.p.lm)[[5]][1],
        gddpp.p.gee$coef[1],gddpp.p.gee$coef[1]+gddpp.p.gee$coef[2],
        gddpp.p.gee$coef[1]+gddpp.p.gee$coef[3],
        foo(gddpp.p.gee)[11])
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        gddpp.i.lm$coef[1]+gddpp.i.lm$coef[3], anova(gddpp.i.lm)[[5]][1],
        gddpp.i.gee$coef[1],gddpp.i.gee$coef[1]+gddpp.i.gee$coef[2],
        gddpp.i.gee$coef[1]+gddpp.i.gee$coef[3],
        foo(gddpp.i.gee)[11])
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         mo3.lm$coef[1]+mo3.lm$coef[3], anova(mo3.lm)[[5]][1], mo3.gee$coef[1],
         mo3.gee$coef[1]+mo3.gee$coef[2], mo3.gee$coef[1]+mo3.gee$coef[3],
         foo(mo3.gee)[11])
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        mat.lm$coef[1]+mat.lm$coef[3],  anova(mat.lm)[[5]][1], mat.gee$coef[1],
        mat.gee$coef[1]+mat.gee$coef[2], mat.gee$coef[1]+mat.gee$coef[3],
        foo(mat.gee)[11])
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        doy.lm$coef[1]+doy.lm$coef[3],  anova(doy.lm)[[5]][1], doy.gee$coef[1],
        doy.gee$coef[1]+doy.gee$coef[2], doy.gee$coef[1]+doy.gee$coef[3],
        foo(doy.gee)[11])
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("lm:ex-inv", "ex-non","nat-non", "p", "ee:ex-inv",
         "ex-non", "nat-non", "p")
return(memodels)
}




##
## as above, but using PGLM instead of GEE

massivetable.pglm <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
         data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp,  vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(meanFFD~native.exotic,
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1], coef(ffd.gee)[1]+coef(ffd.gee)[2],
        pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
         summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1], coef(gdd.gee)[1]+coef(gdd.gee)[2],
         pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
         summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
         coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
         summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
         coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
         summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
         coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1], coef(mo3.gee)[1]+coef(mo3.gee)[2],
         pglm.pval(mo3.gee))
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
          summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1], coef(mat.gee)[1]+coef(mat.gee)[2],
          pglm.pval(mat.gee))
    daterrows <- rbind(meanffd, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c( "LM: exotic", "LM: native", "LM: p", "PGLM:lam", 
        "PGLM: exo", "PGLM: nat", "PGLM: p")
return(memodels)
}

moremassivetable.pglm <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(native.exotic), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(doychange~as.factor(native.exotic), data=goofile$doymod)
    doy.gee <- pglmEstLambda(doychange~as.factor(native.exotic),
        data=goofile$doymod, vcv.doy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2], pglm.pval(mo3.gee))
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], pglm.pval(mat.gee))
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        summary(doy.lm)$coef[8], doy.gee$lambda,  coef(doy.gee)[1], coef(doy.gee)[1]+coef(doy.gee)[2],
        pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod,
        gddpp.imod, mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("LM: exotic", "LM: native", "LM: p", "PGLM:lam", 
        "PGLM: exo", "PGLM: nat", "PGLM: p")
return(memodels)
}


moremassivetable.overtime.pglm <- function(goofile, doychange){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(native.exotic), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(changeovertime~as.factor(native.exotic), data=goofile$doymod)
    doy.gee <- pglmEstLambda(changeovertime~as.factor(native.exotic),
          data=goofile$doymod, vcv.doy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2],
        pglm.pval(mo3.gee))
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], pglm.pval(mat.gee))
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        summary(doy.lm)$coef[8], doy.gee$lambda, coef(doy.gee)[1],
        coef(doy.gee)[1]+coef(doy.gee)[2], pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod,
          gddpp.imod, mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c( "LM: exotic", "LM: native", "LM: p", "PGLM:lam",
        "PGLM: exo", "PGLM: nat", "PGLM: p")
return(memodels)
}

moremassivetable.overtime.pglm.forDC <- function(goofile, doychange){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    # ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
    #    vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    doy.lm <- lm(changeovertime~as.factor(native.exotic), data=goofile$doymod)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8],99, 99, 99, 99)
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8],99, 99, 99, 99)
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8],99, 99, 99, 99)
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8],99, 99, 99, 99)
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8],99, 99, 99, 99)
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8],99, 99, 99, 99)
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8],99, 99, 99, 99)
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        summary(doy.lm)$coef[8],99, 99, 99, 99)
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod,
          gddpp.imod, mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c( "LM: exotic", "LM: native", "LM: p", "PGLM:lam",
        "PGLM: exo", "PGLM: nat", "PGLM: p")
return(memodels)
}


moremassivetable.fc.pglm <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(full.code), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(full.code), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(full.code),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(full.code),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(full.code), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(full.code),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(full.code), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(full.code), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(full.code), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(full.code), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(full.code), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(full.code),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(doychange~as.factor(full.code), data=goofile$doymod)
    doy.gee <- pglmEstLambda(doychange~as.factor(full.code),
        data=goofile$doymod, vcv.doy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        ffd.lm$coef[1]+ffd.lm$coef[3], anova(ffd.lm)[[5]][1],coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], coef(ffd.gee)[1]+coef(ffd.gee)[3],
        pglm.pval(ffd.gee))          
    gddmod <-  c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        gdd.lm$coef[1]+gdd.lm$coef[3], anova(gdd.lm)[[5]][1],coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], coef(gdd.gee)[1]+coef(gdd.gee)[3],
        pglm.pval(gdd.gee))
   gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        gddpp.lm$coef[1]+gddpp.lm$coef[3], anova(gddpp.lm)[[5]][1], coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], coef(gddpp.gee)[1]+coef(gddpp.gee)[3],
        pglm.pval(gddpp.gee))  
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        gddpp.p.lm$coef[1]+gddpp.p.lm$coef[3], anova(gddpp.p.lm)[[5]][1], coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[3],
        pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        gddpp.i.lm$coef[1]+gddpp.i.lm$coef[3], anova(gddpp.i.lm)[[5]][1], coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[3],
        pglm.pval(gddpp.i.gee)) 
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         mo3.lm$coef[1]+mo3.lm$coef[3], anova(mo3.lm)[[5]][1], coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2], coef(mo3.gee)[1]+coef(mo3.gee)[3],
        pglm.pval(mo3.gee)) 
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        mat.lm$coef[1]+mat.lm$coef[3],  anova(mat.lm)[[5]][1], coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], coef(mat.gee)[1]+coef(mat.gee)[3],
        pglm.pval(mat.gee))   
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        doy.lm$coef[1]+doy.lm$coef[3],  anova(doy.lm)[[5]][1], coef(doy.gee)[1],
        coef(doy.gee)[1]+coef(doy.gee)[2], coef(doy.gee)[1]+coef(doy.gee)[3],
        pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("lm:ex-inv", "ex-non","nat-non", "p", "pgl:ex-inv",
         "ex-non", "nat-non", "p")
return(memodels)
}

moremassivetable.fc.overtime.pglm <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(full.code), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(full.code), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(full.code),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(full.code),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(full.code), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(full.code),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(full.code), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(full.code), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(full.code), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(full.code), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(full.code), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(full.code),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(changeovertime~as.factor(full.code), data=goofile$doymod)
    doy.gee <- pglmEstLambda(changeovertime~as.factor(full.code),
        data=goofile$doymod, vcv.doy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        ffd.lm$coef[1]+ffd.lm$coef[3], anova(ffd.lm)[[5]][1],coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], coef(ffd.gee)[1]+coef(ffd.gee)[3],
        pglm.pval(ffd.gee))          
    gddmod <-  c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        gdd.lm$coef[1]+gdd.lm$coef[3], anova(gdd.lm)[[5]][1],coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], coef(gdd.gee)[1]+coef(gdd.gee)[3],
        pglm.pval(gdd.gee))
   gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        gddpp.lm$coef[1]+gddpp.lm$coef[3], anova(gddpp.lm)[[5]][1], coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], coef(gddpp.gee)[1]+coef(gddpp.gee)[3],
        pglm.pval(gddpp.gee))  
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        gddpp.p.lm$coef[1]+gddpp.p.lm$coef[3], anova(gddpp.p.lm)[[5]][1], coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[3],
        pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        gddpp.i.lm$coef[1]+gddpp.i.lm$coef[3], anova(gddpp.i.lm)[[5]][1], coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[3],
        pglm.pval(gddpp.i.gee)) 
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         mo3.lm$coef[1]+mo3.lm$coef[3], anova(mo3.lm)[[5]][1], coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2], coef(mo3.gee)[1]+coef(mo3.gee)[3],
        pglm.pval(mo3.gee)) 
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        mat.lm$coef[1]+mat.lm$coef[3],  anova(mat.lm)[[5]][1], coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], coef(mat.gee)[1]+coef(mat.gee)[3],
        pglm.pval(mat.gee))   
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        doy.lm$coef[1]+doy.lm$coef[3],  anova(doy.lm)[[5]][1], coef(doy.gee)[1],
        coef(doy.gee)[1]+coef(doy.gee)[2], coef(doy.gee)[1]+coef(doy.gee)[3],
        pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("lm:ex-inv", "ex-non","nat-non", "p", "pgl:ex-inv",
         "ex-non", "nat-non", "p")
return(memodels)
}

##
## build all the functions without MPS ##
##


konzagoo.nosig <- function(sitefile, doyfile, treehere){
goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    sma.dat <- subset(as.data.frame(sitefile[["sma"]]),is.na(native.exotic)==FALSE)
    sma.dat <- subset(sma.dat, native.exotic!="native/exotic" & native.exotic!="")
    smab.dat <- subset(as.data.frame(sitefile[["smab"]]),is.na(native.exotic)==FALSE)
    smab.dat <- subset(smab.dat, native.exotic!="native/exotic" & native.exotic!="")
    smanoout.dat <- subset(sma.dat, daysperdev<40)
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperC"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",  "daysperC"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    goo$smamod <- matchtreedat(treehere, sma.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",  "daysperdev", "Month"))
    goo$tree.sma <- matchtreephy(treehere, sma.dat,"latbi")
    goo$smabmod <- matchtreedat(treehere, smab.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",  "daysperdev", "Month"))
    goo$tree.smab <- matchtreephy(treehere, smab.dat,"latbi")
    goo$smanooutmod <- matchtreedat(treehere, smanoout.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",  "daysperdev", "Month"))
    goo$tree.smanoout <- matchtreephy(treehere, smanoout.dat,"latbi")
    # goo$doych <- matchtreedat(treehere, doyzer,"latbi",
    #    c("native.exotic", "changeovertime"))
    # goo$treedoy <- matchtreephy(treehere, goo$doych,"latbi")
    return(goo)
}

mostsitesgoo.nosig <- function(sitefile, doyfile, doychangetype, treehere){
    goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
         c("native.exotic", "full.code", "meanFFD", "daysperHI"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "full.code", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
         c("native.exotic", "full.code", "meanFFD", "daysperC"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "full.code", "meanFFD", "native.exotic12",
        "inv.non12", "daysperC"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    goo$doymod <- matchtreedat(treehere, doyzer,"latbi",
        c("native.exotic", "full.code", doychangetype))
    goo$treedoy <- matchtreephy(treehere, doyzer,"latbi")
    return(goo)
}

dcgoo.nosig <- function(sitefile, doyfile, doychangetype, treehere){
    goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
         c("native.exotic", "meanFFD", "daysperHI"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
         c("native.exotic", "meanFFD", "daysperC"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",
        "daysperC"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    goo$doymod <- matchtreedat(treehere, doyzer,"latbi",
        c("native.exotic", doychangetype))
    goo$treedoy <- matchtreephy(treehere, doyzer,"latbi")
    return(goo)
}



##
##
## Trying something out -- can we adapt this to do DC cultivates vs. not cultivated

dcgoo.nosig.cul <- function(sitefile, doyfile, doychangetype, treehere){
    goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
         c("native.exotic", "meanFFD", "daysperHI", "CUL"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER", "CUL"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
         c("native.exotic", "meanFFD", "daysperC", "CUL"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",
        "daysperC", "CUL"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    goo$doymod <- matchtreedat(treehere, doyzer,"latbi",
        c("native.exotic", doychangetype, "CUL"))
    goo$treedoy <- matchtreephy(treehere, doyzer,"latbi")
    return(goo)
}

moremassivetable.pglm.DCcul <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(CUL), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(CUL), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(CUL), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(CUL),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(CUL), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(CUL),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(CUL), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(CUL),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(CUL), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(CUL), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(CUL), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(CUL), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(CUL), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(CUL),
        data=goofile$matmod, vcv.mat)
    ## all the species that made it to doychange are N (non-cultivated)
    # doy.lm <- lm(changeovertime~as.factor(CUL), data=goofile$doymod)
    # doy.gee <- pglmEstLambda(changeovertime~as.factor(CUL),
    #    data=goofile$doymod, vcv.doy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2], pglm.pval(mo3.gee))
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], pglm.pval(mat.gee))
    # doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
    #    summary(doy.lm)$coef[8], doy.gee$lambda,  coef(doy.gee)[1], coef(doy.gee)[1]+coef(doy.gee)[2],
    #    pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, gddmod, gddppmod, gddpp.pmod,
        gddpp.imod, mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("LM: cul", "LM: noncul", "LM: p", "PGLM:lam", 
        "PGLM: cul", "PGLM: noncul", "PGLM: p")
return(memodels)
}



##
## PGLM table f(x)s adding in PRECIP-only model
##

massivetable.wprecip.pglm <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.pp <- vcv.phylo(goofile$tree.pp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
         data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp,  vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    pp.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$pp)
    pp.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$pp, vcv.pp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(meanFFD~native.exotic,
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1], coef(ffd.gee)[1]+coef(ffd.gee)[2],
        pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
         summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1], coef(gdd.gee)[1]+coef(gdd.gee)[2],
         pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
         summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
         coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    ppmod <- c(pp.lm$coef[1], pp.lm$coef[1]+pp.lm$coef[2],
         summary(pp.lm)$coef[8], pp.gee$lambda, coef(pp.gee)[1],
         coef(pp.gee)[1]+coef(pp.gee)[2], pglm.pval(pp.gee))
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
         summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
         coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
         summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
         coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1], coef(mo3.gee)[1]+coef(mo3.gee)[2],
         pglm.pval(mo3.gee))
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
          summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1], coef(mat.gee)[1]+coef(mat.gee)[2],
          pglm.pval(mat.gee))
    daterrows <- rbind(meanffd, gddmod, gddppmod, ppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec)", "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c( "LM: exotic", "LM: native", "LM: p", "PGLM:lam", 
        "PGLM: exo", "PGLM: nat", "PGLM: p")
return(memodels)
}

moremassivetable.wprecip.pglm <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    vcv.pp <- vcv.phylo(goofile$tree.pp, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    pp.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$pp)
    pp.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$pp, vcv.pp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(native.exotic), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(doychange~as.factor(native.exotic), data=goofile$doymod)
    doy.gee <- pglmEstLambda(doychange~as.factor(native.exotic),
        data=goofile$doymod, vcv.doy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    ppmod <- c(pp.lm$coef[1], pp.lm$coef[1]+pp.lm$coef[2],
        summary(pp.lm)$coef[8], pp.gee$lambda, coef(pp.gee)[1],
        coef(pp.gee)[1]+coef(pp.gee)[2], pglm.pval(pp.gee))
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2], pglm.pval(mo3.gee))
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], pglm.pval(mat.gee))
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        summary(doy.lm)$coef[8], doy.gee$lambda,  coef(doy.gee)[1], coef(doy.gee)[1]+coef(doy.gee)[2],
        pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, ppmod, gddpp.pmod,
        gddpp.imod, mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , "sensitivity (prec)",
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("LM: exotic", "LM: native", "LM: p", "PGLM:lam", 
        "PGLM: exo", "PGLM: nat", "PGLM: p")
return(memodels)
}


moremassivetable.wprecip.overtime.pglm <- function(goofile, doychange){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.pp <- vcv.phylo(goofile$tree.pp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    pp.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$pp)
    pp.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$pp, vcv.pp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(native.exotic), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(changeovertime~as.factor(native.exotic), data=goofile$doymod)
    doy.gee <- pglmEstLambda(changeovertime~as.factor(native.exotic),
          data=goofile$doymod, vcv.doy)
    meanffd <- c(ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    ppmod <- c(pp.lm$coef[1], pp.lm$coef[1]+pp.lm$coef[2],
        summary(pp.lm)$coef[8], pp.gee$lambda, coef(pp.gee)[1],
        coef(pp.gee)[1]+coef(pp.gee)[2], pglm.pval(pp.gee))
    gddpp.imod <- c(gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2],
        pglm.pval(mo3.gee))
    matmod <- c(mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], pglm.pval(mat.gee))
    doymodel <- c(doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        summary(doy.lm)$coef[8], doy.gee$lambda, coef(doy.gee)[1],
        coef(doy.gee)[1]+coef(doy.gee)[2], pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, ppmod, gddpp.pmod,
          gddpp.imod, mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , "sensitivity (prec)",
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c( "LM: exotic", "LM: native", "LM: p", "PGLM:lam",
        "PGLM: exo", "PGLM: nat", "PGLM: p")
return(memodels)
}

konzagoo.wprecip.nosig <- function(sitefile, doyfile, treehere){
    goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    pp.dat <- subset(as.data.frame(sitefile[["pp"]]), is.na(native.exotic)==FALSE)
    pp.dat <- subset(pp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$pp <- matchtreedat(treehere, pp.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperMI"))
    goo$tree.pp <- matchtreephy(treehere, pp.dat,"latbi")   
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperC"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",  "daysperC"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    # goo$doych <- matchtreedat(treehere, doyzer,"latbi",
    #    c("native.exotic", "changeovertime"))
    # goo$treedoy <- matchtreephy(treehere, goo$doych,"latbi")
    return(goo)
}

mostsitesgoo.wprecip.nosig <- function(sitefile, doyfile, doychangetype, treehere){
    goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    pp.dat <- subset(as.data.frame(sitefile[["pp"]]),is.na(native.exotic)==FALSE)
    pp.dat <- subset(pp.dat, native.exotic!="native/exotic" & native.exotic!="")    
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
         c("native.exotic", "full.code", "meanFFD", "daysperHI"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$pp <- matchtreedat(treehere, pp.dat,"latbi",
         c("native.exotic", "full.code", "meanFFD", "daysperMI"))
    goo$tree.pp <- matchtreephy(treehere, pp.dat,"latbi")   
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "full.code", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
         c("native.exotic", "full.code", "meanFFD", "daysperC"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "full.code", "meanFFD", "native.exotic12",
        "inv.non12", "daysperC"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    goo$doymod <- matchtreedat(treehere, doyzer,"latbi",
        c("native.exotic", "full.code", doychangetype))
    goo$treedoy <- matchtreephy(treehere, doyzer,"latbi")
    return(goo)
}

dcgoo.wprecip.nosig <- function(sitefile, doyfile, doychangetype, treehere){
    goo <- list()
    # get all the files together
    gdd.dat <- subset(as.data.frame(sitefile[["gdd"]]),is.na(native.exotic)==FALSE)
    gdd.dat <- subset(gdd.dat, native.exotic!="native/exotic" & native.exotic!="")
    pp.dat <- subset(as.data.frame(sitefile[["pp"]]),is.na(native.exotic)==FALSE)
    pp.dat <- subset(pp.dat, native.exotic!="native/exotic" & native.exotic!="")
    gddpp.dat <- subset(as.data.frame(sitefile[["gddpp"]]),
        is.na(native.exotic)==FALSE)
    gddpp.dat <- subset(gddpp.dat, native.exotic!="native/exotic" & native.exotic!="")
    mo3.dat <- subset(as.data.frame(sitefile[["mo3"]]),is.na(native.exotic)==FALSE)
    mo3.dat <- subset(mo3.dat, native.exotic!="native/exotic" & native.exotic!="")
    mat.dat <- subset(as.data.frame(sitefile[["mat"]]),is.na(native.exotic)==FALSE)
    mat.dat <- subset(mat.dat, native.exotic!="native/exotic" & native.exotic!="")
    doyzer <-subset(doyfile, is.na(native.exotic)==FALSE)
    # get the tree
    treehere <- treehere
    # now match everything up
    goo$gdd <- matchtreedat(treehere, gdd.dat,"latbi",
         c("native.exotic", "meanFFD", "daysperHI"))
    goo$tree.gdd <- matchtreephy(treehere, gdd.dat,"latbi")
    goo$pp <- matchtreedat(treehere, pp.dat,"latbi",
         c("native.exotic", "meanFFD", "daysperMI"))
    goo$tree.pp <- matchtreephy(treehere, pp.dat,"latbi")  
    goo$gddpp <- matchtreedat(treehere, gddpp.dat,"latbi",
        c("native.exotic", "meanFFD", "daysperHI",
        "daysperMI", "daysperINTER"))
    goo$tree.gddpp <- matchtreephy(treehere, gddpp.dat,"latbi")
    goo$mo3 <- matchtreedat(treehere, mo3.dat,"latbi",
         c("native.exotic", "meanFFD", "daysperC"))
    goo$tree.mo3 <- matchtreephy(treehere, mo3.dat,"latbi")
    goo$matmod <- matchtreedat(treehere, mat.dat,"latbi",
        c("native.exotic", "meanFFD", "native.exotic12",
        "daysperC"))
    goo$tree.mat <- matchtreephy(treehere, mat.dat,"latbi")
    goo$doymod <- matchtreedat(treehere, doyzer,"latbi",
        c("native.exotic", doychangetype))
    goo$treedoy <- matchtreephy(treehere, doyzer,"latbi")
    return(goo)
}


##
## New tables that add in DF for LM
## sure this is ugly but at least all the old Sweave code will keep running
##


massivetable.pglmDF <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
         data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp,  vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(meanFFD~native.exotic,
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    ##
    meanffd <- c(ffd.lm$df, ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1], coef(ffd.gee)[1]+coef(ffd.gee)[2],
        pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$df, gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
         summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1], coef(gdd.gee)[1]+coef(gdd.gee)[2],
         pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$df, gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
         summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
         coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$df, gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
         summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
         coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$df, gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
         summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
         coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$df, mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1], coef(mo3.gee)[1]+coef(mo3.gee)[2],
         pglm.pval(mo3.gee))
    matmod <- c(mat.lm$df, mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
          summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1], coef(mat.gee)[1]+coef(mat.gee)[2],
          pglm.pval(mat.gee))
    daterrows <- rbind(meanffd, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("DF",  "LM: exotic", "LM: native", "LM: p", "PGLM: lambda", 
        "PGLM: exotic", "PGLM: native", "PGLM: p")
return(memodels)
}

moremassivetable.pglmDF <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(native.exotic), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(doychange~as.factor(native.exotic), data=goofile$doymod)
    doy.gee <- pglmEstLambda(doychange~as.factor(native.exotic),
        data=goofile$doymod, vcv.doy)
    ##
    meanffd <- c(ffd.lm$df, ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$df, gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$df, gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$df, gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$df, gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$df, mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2], pglm.pval(mo3.gee))
    matmod <- c(mat.lm$df, mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], pglm.pval(mat.gee))
    doymodel <- c(doy.lm$df, doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        summary(doy.lm)$coef[8], doy.gee$lambda,  coef(doy.gee)[1], coef(doy.gee)[1]+coef(doy.gee)[2],
        pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod,
        gddpp.imod, mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("DF",  "LM: exotic", "LM: native", "LM: p", "PGLM: lambda", 
        "PGLM: exotic", "PGLM: native", "PGLM: p")
return(memodels)
}


moremassivetable.overtime.pglmDF <- function(goofile, doychange){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(native.exotic), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(changeovertime~as.factor(native.exotic), data=goofile$doymod)
    doy.gee <- pglmEstLambda(changeovertime~as.factor(native.exotic),
          data=goofile$doymod, vcv.doy)
    ##
    meanffd <- c(ffd.lm$df, ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$df, gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$df, gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$df, gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$df, gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$df, mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
        summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2],
        pglm.pval(mo3.gee))
    matmod <- c(mat.lm$df, mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], pglm.pval(mat.gee))
    doymodel <- c(doy.lm$df, doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        summary(doy.lm)$coef[8], doy.gee$lambda, coef(doy.gee)[1],
        coef(doy.gee)[1]+coef(doy.gee)[2], pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod,
          gddpp.imod, mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("DF",  "LM: exotic", "LM: native", "LM: p", "PGLM: lambda", 
        "PGLM: exotic", "PGLM: native", "PGLM: p")
return(memodels)
}


moremassivetable.fc.pglmDF <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(full.code), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(full.code), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(full.code),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(full.code),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(full.code), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(full.code),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(full.code), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(full.code), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(full.code), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(full.code), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(full.code), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(full.code),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(doychange~as.factor(full.code), data=goofile$doymod)
    doy.gee <- pglmEstLambda(doychange~as.factor(full.code),
        data=goofile$doymod, vcv.doy)
    ##
    meanffd <- c(ffd.lm$df, ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        ffd.lm$coef[1]+ffd.lm$coef[3], anova(ffd.lm)[[5]][1],coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], coef(ffd.gee)[1]+coef(ffd.gee)[3],
        pglm.pval(ffd.gee))          
    gddmod <-  c(gdd.lm$df, gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        gdd.lm$coef[1]+gdd.lm$coef[3], anova(gdd.lm)[[5]][1],coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], coef(gdd.gee)[1]+coef(gdd.gee)[3],
        pglm.pval(gdd.gee))
   gddppmod <- c(gddpp.lm$df, gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        gddpp.lm$coef[1]+gddpp.lm$coef[3], anova(gddpp.lm)[[5]][1], coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], coef(gddpp.gee)[1]+coef(gddpp.gee)[3],
        pglm.pval(gddpp.gee))  
    gddpp.pmod <- c(gddpp.p.lm$df, gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        gddpp.p.lm$coef[1]+gddpp.p.lm$coef[3], anova(gddpp.p.lm)[[5]][1], coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[3],
        pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$df, gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        gddpp.i.lm$coef[1]+gddpp.i.lm$coef[3], anova(gddpp.i.lm)[[5]][1], coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[3],
        pglm.pval(gddpp.i.gee)) 
    mo3mod <- c(mo3.lm$df, mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         mo3.lm$coef[1]+mo3.lm$coef[3], anova(mo3.lm)[[5]][1], coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2], coef(mo3.gee)[1]+coef(mo3.gee)[3],
        pglm.pval(mo3.gee)) 
    matmod <- c(mat.lm$df, mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        mat.lm$coef[1]+mat.lm$coef[3],  anova(mat.lm)[[5]][1], coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], coef(mat.gee)[1]+coef(mat.gee)[3],
        pglm.pval(mat.gee))   
    doymodel <- c(doy.lm$df, doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        doy.lm$coef[1]+doy.lm$coef[3],  anova(doy.lm)[[5]][1], coef(doy.gee)[1],
        coef(doy.gee)[1]+coef(doy.gee)[2], coef(doy.gee)[1]+coef(doy.gee)[3],
        pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("DF", "lm:ex-inv", "ex-non","nat-non", "p", "pgl:ex-inv",
         "ex-non", "nat-non", "p")
return(memodels)
}

moremassivetable.fc.overtime.pglmDF <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.doy <- vcv.phylo(goofile$treedoy, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(full.code), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(full.code), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(full.code),
        data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(full.code), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(full.code),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(full.code), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(full.code),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(full.code), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(full.code), 
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(full.code), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(full.code), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(full.code), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(full.code),
        data=goofile$matmod, vcv.mat)
    doy.lm <- lm(changeovertime~as.factor(full.code), data=goofile$doymod)
    doy.gee <- pglmEstLambda(changeovertime~as.factor(full.code),
        data=goofile$doymod, vcv.doy)
    ##
    meanffd <- c(ffd.lm$df, ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        ffd.lm$coef[1]+ffd.lm$coef[3], anova(ffd.lm)[[5]][1],coef(ffd.gee)[1],
        coef(ffd.gee)[1]+coef(ffd.gee)[2], coef(ffd.gee)[1]+coef(ffd.gee)[3],
        pglm.pval(ffd.gee))          
    gddmod <-  c(gdd.lm$df, gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
        gdd.lm$coef[1]+gdd.lm$coef[3], anova(gdd.lm)[[5]][1],coef(gdd.gee)[1],
        coef(gdd.gee)[1]+coef(gdd.gee)[2], coef(gdd.gee)[1]+coef(gdd.gee)[3],
        pglm.pval(gdd.gee))
   gddppmod <- c(gddpp.lm$df, gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
        gddpp.lm$coef[1]+gddpp.lm$coef[3], anova(gddpp.lm)[[5]][1], coef(gddpp.gee)[1],
        coef(gddpp.gee)[1]+coef(gddpp.gee)[2], coef(gddpp.gee)[1]+coef(gddpp.gee)[3],
        pglm.pval(gddpp.gee))  
    gddpp.pmod <- c(gddpp.p.lm$df, gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
        gddpp.p.lm$coef[1]+gddpp.p.lm$coef[3], anova(gddpp.p.lm)[[5]][1], coef(gddpp.p.gee)[1],
        coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[3],
        pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$df, gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
        gddpp.i.lm$coef[1]+gddpp.i.lm$coef[3], anova(gddpp.i.lm)[[5]][1], coef(gddpp.i.gee)[1],
        coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[3],
        pglm.pval(gddpp.i.gee)) 
    mo3mod <- c(mo3.lm$df, mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         mo3.lm$coef[1]+mo3.lm$coef[3], anova(mo3.lm)[[5]][1], coef(mo3.gee)[1],
        coef(mo3.gee)[1]+coef(mo3.gee)[2], coef(mo3.gee)[1]+coef(mo3.gee)[3],
        pglm.pval(mo3.gee)) 
    matmod <- c(mat.lm$df, mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
        mat.lm$coef[1]+mat.lm$coef[3],  anova(mat.lm)[[5]][1], coef(mat.gee)[1],
        coef(mat.gee)[1]+coef(mat.gee)[2], coef(mat.gee)[1]+coef(mat.gee)[3],
        pglm.pval(mat.gee))   
    doymodel <- c(doy.lm$df, doy.lm$coef[1], doy.lm$coef[1]+doy.lm$coef[2],
        doy.lm$coef[1]+doy.lm$coef[3],  anova(doy.lm)[[5]][1], coef(doy.gee)[1],
        coef(doy.gee)[1]+coef(doy.gee)[2], coef(doy.gee)[1]+coef(doy.gee)[3],
        pglm.pval(doy.gee))
    daterrows <- rbind(meanffd, doymodel, gddmod, gddppmod, gddpp.pmod, gddpp.imod,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD", "flowering time shift",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)",
        "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("DF", "lm:ex-inv", "ex-non","nat-non", "p", "pgl:ex-inv",
         "ex-non", "nat-non", "p")
return(memodels)
}


## for Konza, with soil moisture
massivetable.pglmDF.smkonz <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.konzsma <- vcv.phylo(konzmatic$tree.sma, model="Brownian", corr=FALSE)
    vcv.konzsmab <- vcv.phylo(konzmatic$tree.smab, model="Brownian", corr=FALSE)
    vcv.konzsmanoout <- vcv.phylo(konzmatic$tree.smanoout, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
         data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp,  vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(meanFFD~native.exotic,
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    # merge files
    konzsmxsheds <- merge(konzmatic$smamod, konzmatic$smabmod, by="row.names",
        suffixes=c(".1d", ".20b"))
    konzsmxsheds$daysperdevxsheds <- (konzsmxsheds$daysperdev.1d+konzsmxsheds$daysperdev.20b)/2
    row.names(konzsmxsheds) <- konzsmxsheds$Row.names
    konzsmxsheds.o <- merge(konzmatic$smanooutmod, konzmatic$smabmod, by="row.names",
        suffixes=c(".1d", ".20b"))
    konzsmxsheds.o$daysperdevxsheds <- (konzsmxsheds.o$daysperdev.1d+konzsmxsheds.o$daysperdev.20b)/2
        row.names(konzsmxsheds.o) <- konzsmxsheds.o$Row.names
    konz.smxw.lm <-  lm(daysperdevxsheds~as.factor(native.exotic.1d), data=konzsmxsheds)
    konz.smxw.gee <-pglmEstLambda(daysperdevxsheds~as.factor(native.exotic.1d), 
        data=konzsmxsheds, vcv.konzsmab)
    konz.smxwo.lm <-  lm(daysperdevxsheds~as.factor(native.exotic.1d), data=konzsmxsheds.o)
    konz.smxwo.gee <-pglmEstLambda(daysperdevxsheds~as.factor(native.exotic.1d), 
        data=konzsmxsheds.o, vcv.konzsmanoout)
    ##
    smoutput <- c(konz.smxw.lm$df, konz.smxw.lm$coef[1], konz.smxw.lm$coef[1]+konz.smxw.lm$coef[2],
        summary(konz.smxw.lm)$coef[8], konz.smxw.gee$lambda, coef(konz.smxw.gee)[1], 
        coef( konz.smxw.gee)[1]+coef(konz.smxw.gee)[2], pglm.pval(konz.smxw.gee))
    smoutputo <- c(konz.smxwo.lm$df, konz.smxwo.lm$coef[1], konz.smxwo.lm$coef[1]+konz.smxwo.lm$coef[2],
        summary(konz.smxwo.lm)$coef[8], konz.smxwo.gee$lambda, coef(konz.smxwo.gee)[1], 
        coef( konz.smxwo.gee)[1]+coef(konz.smxwo.gee)[2], pglm.pval(konz.smxwo.gee))
    meanffd <- c(ffd.lm$df, ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1], coef(ffd.gee)[1]+coef(ffd.gee)[2],
        pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$df, gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
         summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1], coef(gdd.gee)[1]+coef(gdd.gee)[2],
         pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$df, gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
         summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
         coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$df, gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
         summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
         coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$df, gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
         summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
         coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$df, mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1], coef(mo3.gee)[1]+coef(mo3.gee)[2],
         pglm.pval(mo3.gee))
    matmod <- c(mat.lm$df, mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
          summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1], coef(mat.gee)[1]+coef(mat.gee)[2],
          pglm.pval(mat.gee))
    daterrows <- rbind(meanffd, gddmod, gddppmod, gddpp.pmod, gddpp.imod, smoutput, smoutputo,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (prec w/gdd)", "sensitivity (intxn)", "sensitivity (sm)",
         "sensitivity (sm, no outlier)","sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("DF",  "LM: exotic", "LM: native", "LM: p", "PGLM: lambda", 
        "PGLM: exotic", "PGLM: native", "PGLM: p")
return(memodels)
}


## for Konza, with soil moisture
massivetable.pglmDF.smkonznoout <- function(goofile){
    vcv.gdd <- vcv.phylo(goofile$tree.gdd, model="Brownian", corr=FALSE)
    vcv.gddpp <- vcv.phylo(goofile$tree.gddpp, model="Brownian", corr=FALSE)
    vcv.mat <- vcv.phylo(goofile$tree.mat, model="Brownian", corr=FALSE)
    vcv.mo3 <- vcv.phylo(goofile$tree.mo3, model="Brownian", corr=FALSE)
    vcv.konzsma <- vcv.phylo(konzmatic$tree.sma, model="Brownian", corr=FALSE)
    vcv.konzsmab <- vcv.phylo(konzmatic$tree.smab, model="Brownian", corr=FALSE)
    vcv.konzsmanoout <- vcv.phylo(konzmatic$tree.smanoout, model="Brownian", corr=FALSE)
    ffd.lm <- lm(meanFFD~as.factor(native.exotic), data=goofile$gdd)
    ffd.gee <- pglmEstLambda(meanFFD~as.factor(native.exotic), data=goofile$gdd,
        vcv.gdd)
    gdd.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gdd)
    gdd.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
         data=goofile$gdd, vcv.gdd)
    gddpp.lm <- lm(daysperHI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.gee <- pglmEstLambda(daysperHI~as.factor(native.exotic),
        data=goofile$gddpp,  vcv.gddpp)
    gddpp.p.lm <- lm(daysperMI~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.p.gee <- pglmEstLambda(daysperMI~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    gddpp.i.lm <- lm(daysperINTER~as.factor(native.exotic), data=goofile$gddpp)
    gddpp.i.gee <- pglmEstLambda(daysperINTER~as.factor(native.exotic),
        data=goofile$gddpp, vcv.gddpp)
    mo3.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$mo3)
    mo3.gee <- pglmEstLambda(daysperC~as.factor(native.exotic), data=goofile$mo3, 
        vcv.mo3)
    mat.lm <- lm(daysperC~as.factor(native.exotic), data=goofile$matmod)
    mat.gee <- pglmEstLambda(daysperC~as.factor(native.exotic),
        data=goofile$matmod, vcv.mat)
    # merge files
    konzsmxsheds <- merge(konzmatic$smamod, konzmatic$smabmod, by="row.names",
        suffixes=c(".1d", ".20b"))
    konzsmxsheds$daysperdevxsheds <- (konzsmxsheds$daysperdev.1d+konzsmxsheds$daysperdev.20b)/2
    row.names(konzsmxsheds) <- konzsmxsheds$Row.names
    konzsmxsheds.o <- merge(konzmatic$smanooutmod, konzmatic$smabmod, by="row.names",
        suffixes=c(".1d", ".20b"))
    konzsmxsheds.o$daysperdevxsheds <- (konzsmxsheds.o$daysperdev.1d+konzsmxsheds.o$daysperdev.20b)/2
        row.names(konzsmxsheds.o) <- konzsmxsheds.o$Row.names
    konz.smxw.lm <-  lm(daysperdevxsheds~as.factor(native.exotic.1d), data=konzsmxsheds)
    konz.smxw.gee <-pglmEstLambda(daysperdevxsheds~as.factor(native.exotic.1d), 
        data=konzsmxsheds, vcv.konzsmab)
    konz.smxwo.lm <-  lm(daysperdevxsheds~as.factor(native.exotic.1d), data=konzsmxsheds.o)
    konz.smxwo.gee <-pglmEstLambda(daysperdevxsheds~as.factor(native.exotic.1d), 
        data=konzsmxsheds.o, vcv.konzsmanoout)
    ##
    smoutput <- c(konz.smxw.lm$df, konz.smxw.lm$coef[1], konz.smxw.lm$coef[1]+konz.smxw.lm$coef[2],
        summary(konz.smxw.lm)$coef[8], konz.smxw.gee$lambda, coef(konz.smxw.gee)[1], 
        coef( konz.smxw.gee)[1]+coef(konz.smxw.gee)[2], pglm.pval(konz.smxw.gee))
    smoutputo <- c(konz.smxwo.lm$df, konz.smxwo.lm$coef[1], konz.smxwo.lm$coef[1]+konz.smxwo.lm$coef[2],
        summary(konz.smxwo.lm)$coef[8], konz.smxwo.gee$lambda, coef(konz.smxwo.gee)[1], 
        coef( konz.smxwo.gee)[1]+coef(konz.smxwo.gee)[2], pglm.pval(konz.smxwo.gee))
    meanffd <- c(ffd.lm$df, ffd.lm$coef[1], ffd.lm$coef[1]+ffd.lm$coef[2],
        summary(ffd.lm)$coef[8], ffd.gee$lambda, coef(ffd.gee)[1], coef(ffd.gee)[1]+coef(ffd.gee)[2],
        pglm.pval(ffd.gee))
    gddmod <- c(gdd.lm$df, gdd.lm$coef[1], gdd.lm$coef[1]+gdd.lm$coef[2],
         summary(gdd.lm)$coef[8], gdd.gee$lambda, coef(gdd.gee)[1], coef(gdd.gee)[1]+coef(gdd.gee)[2],
         pglm.pval(gdd.gee))
    gddppmod <- c(gddpp.lm$df, gddpp.lm$coef[1], gddpp.lm$coef[1]+gddpp.lm$coef[2],
         summary(gddpp.lm)$coef[8], gddpp.gee$lambda, coef(gddpp.gee)[1],
         coef(gddpp.gee)[1]+coef(gddpp.gee)[2], pglm.pval(gddpp.gee))
    gddpp.pmod <- c(gddpp.p.lm$df, gddpp.p.lm$coef[1], gddpp.p.lm$coef[1]+gddpp.p.lm$coef[2],
         summary(gddpp.p.lm)$coef[8], gddpp.p.gee$lambda, coef(gddpp.p.gee)[1],
         coef(gddpp.p.gee)[1]+coef(gddpp.p.gee)[2], pglm.pval(gddpp.p.gee))
    gddpp.imod <- c(gddpp.i.lm$df, gddpp.i.lm$coef[1], gddpp.i.lm$coef[1]+gddpp.i.lm$coef[2],
         summary(gddpp.i.lm)$coef[8], gddpp.i.gee$lambda, coef(gddpp.i.gee)[1],
         coef(gddpp.i.gee)[1]+coef(gddpp.i.gee)[2], pglm.pval(gddpp.i.gee))
    mo3mod <- c(mo3.lm$df, mo3.lm$coef[1], mo3.lm$coef[1]+mo3.lm$coef[2],
         summary(mo3.lm)$coef[8], mo3.gee$lambda, coef(mo3.gee)[1], coef(mo3.gee)[1]+coef(mo3.gee)[2],
         pglm.pval(mo3.gee))
    matmod <- c(mat.lm$df, mat.lm$coef[1], mat.lm$coef[1]+mat.lm$coef[2],
          summary(mat.lm)$coef[8], mat.gee$lambda, coef(mat.gee)[1], coef(mat.gee)[1]+coef(mat.gee)[2],
          pglm.pval(mat.gee))
    daterrows <- rbind(meanffd, gddmod, gddppmod, gddpp.imod, gddpp.pmod, smoutput,
          mo3mod, matmod)
    memodels <- data.frame(daterrows, row.names = c("mean FFD",
        "sensitivity (gdd)","sensitivity (gdd w/prec)" , 
        "sensitivity (intxn)", "sensitivity (prec w/gdd)", "sensitivity (soil moisture)",
         "sensitivity (3 month)", "sensitivity (MAT)"))
    names(memodels) <- c("DF",  "LM: exotic", "LM: native", "LM: p", "PGLM: lambda", 
        "PGLM: exotic", "PGLM: native", "PGLM: p")
return(memodels)
}
