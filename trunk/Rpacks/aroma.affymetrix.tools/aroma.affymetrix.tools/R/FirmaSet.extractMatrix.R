#... passed to subset
setMethodS3( "extractMatrix","FirmaSet",
function (this, units = NULL, ..., field = c("intensities", "stdvs", "pixels"), 
    returnUgcMap = FALSE, verbose = FALSE) 
{
    cdf <- getCdf(this)
    if (is.null(units)) {
        nunits <- nbrOfUnits(cdf)
    }
    else {
        units <- Arguments$getIndices(units, range = c(1, nbrOfUnits(cdf)))
        nunits <- length(units)
    }
    field <- match.arg(field)
    verbose <- Arguments$getVerbose(verbose)
    if (verbose) {
        pushState(verbose)
        on.exit(popState(verbose))
    }
    verbose && enter(verbose, "Getting data for the array set")
    verbose && enter(verbose, "Getting unit-to-cell map")
    cf <- getFile(this, 1)
    ugcMap <- getCellMap(cf, units = units, verbose = less(verbose))
    ugcMap <- subset(ugcMap, ...)
    verbose && exit(verbose)
    if (nrow(ugcMap) == 0) 
        throw("Nothing to return.")
    arrayNames <- getNames(this)
    nbrOfArrays <- length(arrayNames)
    df <- matrix(NA, nrow = nrow(ugcMap), ncol = nbrOfArrays)
    colnames(df) <- arrayNames
    gc <- gc()
    verbose && print(verbose, gc)
    verbose && enter(verbose, "Retrieving sample intensities")
    for (aa in seq_len(nbrOfArrays)) {
        verbose && printf(verbose, "Array %d,\n", aa)
        cf <- getFile(this, aa)
        df[, aa] <- getDataFlat(cf, units = ugcMap, fields = field, 
            verbose = less(verbose))[, field]
        if (aa%%10 == 0) {
            gc <- gc()
            verbose && print(verbose, gc)
        }
    }
    verbose && exit(verbose)
    verbose && exit(verbose)
    if (returnUgcMap) 
        attr(df, "unitGroupCellMap") <- ugcMap
    df
})
