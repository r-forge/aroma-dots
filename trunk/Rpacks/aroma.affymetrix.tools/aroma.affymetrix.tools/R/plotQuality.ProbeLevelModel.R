setMethodS3("plotQuality","ProbeLevelModel",
function(this, arrays=NULL,subset = NULL, verbose = FALSE, plot.it=TRUE,
    type=c("Nuse","Rle"),main = type, xlim=NULL,ylim=NULL,outline=FALSE,las=2,...) 
{
    type<-match.arg(type)
    verbose <- Arguments$getVerbose(verbose)
    ces <- getChipEffects(this)
    cdfMono <- getCdf(ces)
    nbrOfUnits <- nbrOfUnits(cdfMono)
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
    
    mess<-switch(type,"Nuse"="standard errors","Rle"="chip effects")
    verbose && enter(verbose, "Extracting median ",mess)
    field<-switch(type,"Nuse"="stdvs","Rle"="intensities")
    verbose && enter(verbose, "Calculating median ",mess)
#will include the 'blank' indices; note will recalculate, so if done it once, will not save time to call 'subset'
#this is because of how getAverageLog works...
    if(length(units)<nbrOfUnits) indices<-unlist(getCellIndices(cdfMono, units=units)) 
    else indices<-"remaining" #save time and not recalculate if already calculated
    avg <- getAverageLog(ces, field = field, verbose = verbose,mean="median", sd="mad",indices=indices,force=FALSE) 
    verbose && exit(verbose)

    medianValue <- extractMatrix(avg, units = units, field="theta")
    medianValue <- log2(medianValue)
    verbose && exit(verbose)
    
    verbose && enter(verbose, "Extracting ", mess)
    
    field<-switch(type,"Nuse"="sdTheta","Rle"="theta")
    value<-extractMatrix(ces, units = units, field=field)
    value<-log2(value)
    verbose && exit(verbose)
    
    if(is.null(arrays)) arrays<-seq(ces)
    else arrays<-Arguments$getIndices(arrays, range = c(1,length(ces)))
    verbose && enter(verbose, "Calculating summaries for ", length(arrays), " arrays")
    boxplotStats <- list()
    for (kk in arrays) {
        x<-switch(type,"Nuse"=value[,kk]/medianValue,"Rle"=value[,kk]-medianValue)
        boxplotStats[[kk]] <- boxplot.stats(x)
        rm(x)
    }
    rm(avg)
    rm(medianValue)
    rm(value)
    verbose && exit(verbose)
    bxpStats <- list()
    bxpStats[["stats"]] <- do.call("cbind", lapply(boxplotStats, 
        function(x) {
            x[["stats"]]
        }))
    bxpStats[["conf"]] <- do.call("cbind", lapply(boxplotStats, 
        function(x) {
            x[["conf"]]
        }))
    bxpStats[["n"]] <- do.call("c", lapply(boxplotStats, function(x) {
        x[["n"]]
    }))
    bxpStats[["out"]] <- do.call("c", lapply(boxplotStats, function(x) {
        x[["out"]]
    }))
    igroup <- 0
    bxpStats[["group"]] <- do.call("c", lapply(boxplotStats, 
        function(x) {
            igroup <<- igroup + 1
            rep(igroup, length(x[["out"]]))
        }))
    bxpStats[["names"]]<-getNames(ces)[arrays]
    
    if(plot.it) {
        #fix the strange behavior of bxp if outline=FALSE
        addPlotValues<-list(...)
        addPlotChoices<-names(addPlotValues)
        if(!outline){
            if(!("horizontal" %in% addPlotChoices) || !addPlotValues["horizontal"]){ #standard plot
                if(is.null(ylim)) ylim<-range(as.vector(bxpStats[["stats"]]))
                
            }
            else{#horizontal plot
                if(is.null(xlim)) xlim<-range(as.vector(bxpStats[["stats"]]))
            }    
        }
        bxp(bxpStats, main = main, xlim=xlim,ylim=ylim,outline=outline,las=las,...)
    }
    invisible(bxpStats)
})
setMethodS3("plotRle","QualityAssessmentModel",
function(this,...) {
    plm<-getPlm(this)
    plotQuality(plm,type="Rle",...)
})
setMethodS3("plotNuse","QualityAssessmentModel",
function(this,...) {
    plm<-getPlm(this)
    plotQuality(plm,type="Nuse",...)
})
