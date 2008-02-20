#... passed to subset
setMethodS3("extractMatrix","FirmaFile",
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
    ugcMap <- getCellMap(this, units = units, verbose = less(verbose))
    ugcMap <- subset(ugcMap, ...)
    verbose && exit(verbose)
    if (nrow(ugcMap) == 0) throw("Nothing to return.")
    df <- matrix(NA, nrow = nrow(ugcMap), ncol = 1)
    gc <- gc()
    verbose && print(verbose, gc)
    verbose && enter(verbose, "Retrieving sample values") #here
    df[, 1] <- getDataFlat(this, units = ugcMap, fields = field, verbose = less(verbose))[, field]
    verbose && exit(verbose)
    verbose && exit(verbose)
    if (returnUgcMap) attr(df, "unitGroupCellMap") <- ugcMap
    df
})
