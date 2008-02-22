setMethodS3("plotBoxplot","AffymetrixCelFile",function(this, verbose = FALSE,statsArgs=NULL,...){
    boxplotStats<-do.call("boxplotStats", args=c(list(this=this),statsArgs))
    plotBoxplot(boxplotStats,...)
})
setMethodS3("plotBoxplot","AffymetrixCelFile",function(this,subset=1/2,verbose = FALSE,
    field="intensities", trans=log2,types="pm",...){
    verbose <- Arguments$getVerbose(verbose)
    cdf <- getCdf(this)
    verbose && enter(verbose, "Identifying subset of probes")
    suppressWarnings({
        subset <- identifyCells(cdf, indices = subset, types = types, verbose = less(verbose))
    })
    verbose && exit(verbose)
    verbose && enter(verbose, "Plotting the boxplot")
    verbose && cat(verbose, "Array: ", getName(this))
    suppressWarnings({
        verbose && enter(verbose, "Loading probe intensities")
        y <- getData(this, indices = subset, fields = "intensities")
        y <- y$intensities
        verbose && exit(verbose)
        y <- trans(y)
        verbose && cat(verbose, "Plotting")
        out<-boxplot(y,...)
    })
    verbose && exit(verbose)
    invisible(out)
})
