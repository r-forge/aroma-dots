#... passed to boxplot.stats
setMethodS3("boxplotStats","AffymetrixCelFile",
function(this, subset = NULL, verbose = FALSE,trans=log2,types="pm",...) 
{
    cdf <- getCdf(this)
    verbose <- Arguments$getVerbose(verbose)
    verbose && enter(verbose, "Identifying subset of probes")
    suppressWarnings({
        subset <- identifyCells(cdf, indices = subset, types = types, verbose = less(verbose))
    })
    verbose && exit(verbose)
    verbose && enter(verbose, "Calculating Boxplot Statistics")
    verbose && cat(verbose, "Array: ", getName(this))
    suppressWarnings({
        verbose && enter(verbose, "Loading probe intensities")
        y <- getData(this, indices = subset, fields = "intensities")
        y <- y$intensities
        verbose && exit(verbose)
        y <- trans(y)
        verbose && cat(verbose, "Plotting")
        out<-boxplot.stats(y,...)
    })
    verbose && exit(verbose)
    
    return(out)

})
