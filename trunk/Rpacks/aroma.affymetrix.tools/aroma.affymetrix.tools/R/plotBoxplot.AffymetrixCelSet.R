setMethodS3("plotBoxplot","AffymetrixCelSet",function(this,subset=NULL,arrays=NULL, verbose = FALSE,
    field="intensities", trans=log2,types="pm",...){
    nbrOfArrays <- nbrOfArrays(this)
    
    #if (is.null(col)) {
#        col <- seq(length = nbrOfArrays)
#    }
#    else {
#        col <- rep(col, length.out = nbrOfArrays)
#    }
#    if (!is.null(lty)) 
#        lty <- rep(lty, length.out = nbrOfArrays)
#    if (!is.null(lwd)) 
#        lwd <- rep(lwd, length.out = nbrOfArrays)
    verbose <- Arguments$getVerbose(verbose)
    if(!is.null(arrays)){
        if(any(!is.numeric(arrays))) throw("Invalid array: must be numeric")
        else if(any(arrays>nbrOfArrays | arrays<1)) throw("Invalid array")
        arrays<-seq(this)[arrays]
    }
    else arrays<-seq(this)
    cdf <- getCdf(this)
    verbose && enter(verbose, "Identifying subset of probes")
    suppressWarnings({
        subset <- identifyCells(cdf, indices = subset, types = types, ..., verbose = less(verbose))
    })
    verbose && exit(verbose)
    boxplotStats<-list()
    for (cc in 1:length(arrays)) {
        df <- getFile(this, arrays[cc])
        boxplotStats[[cc]]<-plotBoxplot(df, subset = subset, plot=F,add = add, field=field, trans=trans,types=types, verbose = less(verbose))
    }
    bxpStats <- list()
    elementNames <-    c("stats","n","conf","out","group","names") #parts of the list
    for (name in elementNames[1:3]) {
        suppressWarnings(bxpStats[[name]] <- do.call("cbind", lapply(boxplotStats, function(x) {x[[name]]})))
    }
    bxpStats$out<-unlist(lapply(boxplotStats,function(x){x$out}))
    nOutPerGroup<-unlist(lapply(boxplotStats,function(x){length(x$out)}))
    bxpStats$group<-rep(1:length(arrays),times=nOutPerGroup)
    bxpStats$names<-getNames(this)[arrays]
    bxp(bxpStats,...)
    invisible(bxpStats)
})
