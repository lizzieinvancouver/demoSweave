## Goes with analysesataglance.R ##
## wine vintage data ## 

getclmrows <- function(winecolname, climatecolname, breakyear){
    bordtbf80data <- subset(bbshapetemp, year<=breakyear)
    bordtbf80 <- clm(factor(bordtbf80data[[winecolname]],ordered=TRUE)~
        bordtbf80data[[climatecolname]])
    bordtaft80data <- subset(bbshapetemp, year>breakyear)
    bordtaft80 <- clm(factor(bordtaft80data[[winecolname]],
        ordered=TRUE)~bordtaft80data[[climatecolname]])
    bordprbf80data <- subset(bbshapeprec, year<=breakyear)
    bordprbf80 <- clm(factor(bordprbf80data[[winecolname]],ordered=TRUE)~
        bordprbf80data[[climatecolname]])
    bordpraft80data <- subset(bbshapeprec, year>breakyear)
    bordpraft80 <- clm(factor(bordpraft80data[[winecolname]],ordered=TRUE)~
        bordpraft80data[[climatecolname]])
    bordpdbf80data <- subset(bbshapepdsi, year<=breakyear)
    bordpdbf80 <- clm(factor(bordpdbf80data[[winecolname]],
        ordered=TRUE)~bordpdbf80data[[climatecolname]])
    bordpdaft80data <- subset(bbshapepdsi, year>breakyear)
    bordpdaft80 <- clm(factor(bordpdaft80data[[winecolname]], ordered=TRUE)~
        bordpdaft80data[[climatecolname]])
    bordghdbf80data <- subset(bbshapeghd, year<=breakyear)
    bordghdbf80 <- clm(factor(bordghdbf80data[[winecolname]], ordered=TRUE)~
        bordghdbf80data[[climatecolname]])
    bordghdaft80data <- subset(bbshapeghd, year>breakyear)
    bordghdaft80 <- clm(factor(bordghdaft80data[[winecolname]], ordered=TRUE)~
        bordghdaft80data[[climatecolname]])

    bordredrow <-cbind(tail(bordghdbf80$coefficients, n=1),
        tail(bordghdaft80$coefficients,n=1),
        tail(bordtbf80$coefficients, n=1),
        tail(bordtaft80$coefficients, n=1),
        tail(bordprbf80$coefficients, n=1),
        tail(bordpraft80$coefficients, n=1),
        tail(bordpdbf80$coefficients, n=1),
        tail(bordpdaft80$coefficients,n=1))
    return(bordredrow)
}


getclmrows.wpval <- function(winecolname, climatecolname, breakyear){
    bordtbf80data <- subset(bbshapetemp, year<=breakyear)
    bordtbf80 <- clm(factor(bordtbf80data[[winecolname]],
        ordered=TRUE)~bordtbf80data[[climatecolname]], ordered=TRUE)
    bordtaft80data <- subset(bbshapetemp, year>breakyear)
    bordtaft80 <- clm(factor(bordtaft80data[[winecolname]],
        ordered=TRUE)~bordtaft80data[[climatecolname]], ordered=TRUE)
    bordprbf80data <- subset(bbshapeprec, year<=breakyear)
    bordprbf80 <- clm(factor(bordprbf80data[[winecolname]],
        ordered=TRUE)~bordprbf80data[[climatecolname]], ordered=TRUE)
    bordpraft80data <- subset(bbshapeprec, year>breakyear)
    bordpraft80 <- clm(factor(bordpraft80data[[winecolname]],
        ordered=TRUE)~bordpraft80data[[climatecolname]], ordered=TRUE)
    bordpdbf80data <- subset(bbshapepdsi, year<=breakyear)
    bordpdbf80 <- clm(factor(bordpdbf80data[[winecolname]],
        ordered=TRUE)~bordpdbf80data[[climatecolname]], ordered=TRUE)
    bordpdaft80data <- subset(bbshapepdsi, year>breakyear)
    bordpdaft80 <- clm(factor(bordpdaft80data[[winecolname]],
        ordered=TRUE)~bordpdaft80data[[climatecolname]], ordered=TRUE)
    bordghdbf80data <- subset(bbshapeghd, year<=breakyear)
    bordghdbf80 <- clm(factor(bordghdbf80data[[winecolname]],
        ordered=TRUE)~bordghdbf80data[[climatecolname]], ordered=TRUE)
    bordghdaft80data <- subset(bbshapeghd, year>breakyear)
    bordghdaft80 <- clm(factor(bordghdaft80data[[winecolname]],
        ordered=TRUE)~bordghdaft80data[[climatecolname]], ordered=TRUE)
    
    bordredrow <-cbind(tail(coef(summary(bordghdbf80))[,4], n=1),
        tail(coef(summary(bordghdaft80))[,4],n=1),
        tail(coef(summary(bordtbf80))[,4], n=1),
        tail(coef(summary(bordtaft80))[,4], n=1),
        tail(coef(summary(bordprbf80))[,4], n=1),
        tail(coef(summary(bordpraft80))[,4], n=1),
        tail(coef(summary(bordpdbf80))[,4], n=1),
        tail(coef(summary(bordpdaft80))[,4],n=1))
    
     return(bordredrow)
}


getlmrows.wpval <- function(winecolname, climatecolname, breakyear){
    tbf80data <- subset(bbshapetemp, year<=breakyear)
    tbf80 <- lm(tbf80data[[winecolname]]~tbf80data[[climatecolname]])
    taft80data <- subset(bbshapetemp, year>breakyear)
    taft80 <- lm(taft80data[[winecolname]]~taft80data[[climatecolname]])
    prbf80data <- subset(bbshapeprec, year<=breakyear)
    prbf80 <- lm(prbf80data[[winecolname]]~prbf80data[[climatecolname]])
    praft80data <- subset(bbshapeprec, year>breakyear)
    praft80 <- lm(praft80data[[winecolname]]~praft80data[[climatecolname]])
    pdbf80data <- subset(bbshapepdsi, year<=breakyear)
    pdbf80 <- lm(pdbf80data[[winecolname]]~pdbf80data[[climatecolname]])
    pdaft80data <- subset(bbshapepdsi, year>breakyear)
    pdaft80 <- lm(pdaft80data[[winecolname]]~pdaft80data[[climatecolname]])
    ghdbf80data <- subset(bbshapeghd, year<=breakyear)
    ghdbf80 <- lm(ghdbf80data[[winecolname]]~ghdbf80data[[climatecolname]])
    ghdaft80data <- subset(bbshapeghd, year>breakyear)
    ghdaft80 <- lm(ghdaft80data[[winecolname]]~ghdaft80data[[climatecolname]])

    redrow <- cbind(coef(summary(tbf80))[2], coef(summary(taft80))[2],
        coef(summary(prbf80))[2],coef(summary(praft80))[2], 
        coef(summary(pdbf80))[2], coef(summary(pdaft80))[2],
        coef(summary(ghdbf80))[2], coef(summary(ghdaft80))[2])
    redrow.p <- cbind(coef(summary(tbf80))[8], coef(summary(taft80))[8],
        coef(summary(prbf80))[8],coef(summary(praft80))[8],
        coef(summary(pdbf80))[8], coef(summary(pdaft80))[8],
        coef(summary(ghdbf80))[8], coef(summary(ghdaft80))[8])

    return(rbind(redrow, redrow.p))
}
