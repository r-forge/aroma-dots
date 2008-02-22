#... passed to boxplot.stats
setMethodS3("boxplotStats","QualityAssessmentModel",
function(this, arrays=NULL,subset = NULL, verbose = FALSE, type=c("Nuse","Rle"),...) 
{
    type<-match.arg(type)
    plm<-getPlm(this)
    verbose <- Arguments$getVerbose(verbose)
    ces <- getChipEffects(plm)
    nbrOfUnits <- nbrOfUnits(getCdf(ces))
    if (!(is.null(subset))) {
        getFraction <- (length(subset) == 1) && (subset >= 0) && (subset < 1)
        if (!getFraction) {
            units <- Arguments$getIndices(subset, range = c(1, nbrOfUnits))
        }
        else {
            units <- seq(from=1,to=nbrOfUnits,by=1/subset)
        }
    }
    else {
        units <- 1:nbrOfUnits
    }
    verbose && enter(verbose, "Extracting ",type," values")
    calcFun<-switch(type,"Nuse"=extractNuse,"Rle"=extractRle)
    x<-calcFun(ces,units=units,verbose=verbose)
    verbose && exit(verbose)
  
    if(is.null(arrays)) arrays<-seq(ces)
    else arrays<-Arguments$getIndices(arrays, range = c(1,length(ces)))
    verbose && enter(verbose, "Calculating summaries for ", length(arrays), " arrays")
    boxplotStats <- list()
    for (kk in arrays) {
        boxplotStats[[kk]] <- boxplot.stats(x[,kk],...)
    }
    names(boxplotStats)<-getNames(ces)[arrays]
    rm(x)
    verbose && exit(verbose)

    return(boxplotStats)

})
