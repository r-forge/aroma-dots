
fit.FirmaModel<-function (this, units = "remaining", ..., force = FALSE, verbose = FALSE) 
{
    ds <- getDataSet(this$.plm)
    cdf <- getCdf(ds)
    if (this$operateOn=="weights")
      ws <- calculateWeights(this, verbose = verbose)
    else
      ws <- calculateResiduals(this, verbose = verbose)
    nbrOfArrays <- length(ds)
    doRemaining <- FALSE
    if (is.null(units)) {
    }
    else if (is.numeric(units)) {
        units <- Arguments$getIndices(units, range = c(1, nbrOfUnits(cdf)))
    }
    else if (identical(units, "remaining")) {
        doRemaining <- TRUE
    }
    else {
        throw("Unknown mode of argument 'units': ", mode(units))
    }
    force <- Arguments$getLogical(force)
    verbose <- Arguments$getVerbose(verbose)
    if (verbose) {
        pushState(verbose)
        on.exit(popState(verbose))
    }
    verbose && enter(verbose, "Fitting model of class ", class(this)[1], ":")
    verbose && print(verbose, this)
    if (is.null(units)) {
        nbrOfUnits <- nbrOfUnits(cdf)
        units <- 1:nbrOfUnits
    }
    else if (doRemaining) {
        verbose && enter(verbose, "Identifying non-estimated units")
        units <- findUnitsTodo(this, verbose = less(verbose))
        nbrOfUnits <- length(units)
        verbose && exit(verbose)
    }
    else {
        units <- unique(units)
        nbrOfUnits <- length(units)
    }
    verbose && printf(verbose, "Getting FIRMA results for %d units.\n", nbrOfUnits)
    if (!doRemaining) {
        if (force) {
            verbose && printf(verbose, "All of these are forced to be fitted.\n")
        }
        else {
            units <- findUnitsTodo(this, units = units, verbose = less(verbose))
            nbrOfUnits <- length(units)
            verbose && printf(verbose, "Out of these, %d units need to be fitted.\n", nbrOfUnits)
        }
    }
    if (nbrOfUnits == 0) 
        return(NULL)
    fitUnit <- getFitUnitFunction(this)
    fs <- getFirmaScores(this, verbose = less(verbose))
    idxs <- 1:nbrOfUnits
    startTime <- processTime()
    timers <- list(total = 0, read = 0, fit = 0, writeFs = 0, gc = 0)
    verbose && cat(verbose, "Units: ")
    verbose && str(verbose, units)
    unitsFull <- units
    for (kk in seq(ws)) {
        units <- unitsFull
        tTotal <- processTime()
        ff <- getFile(fs, kk)
        #wf <- getFile(ws, kk)
        wf <- getFile(ws, kk)
        verbose && enter(verbose, "Array #", kk, ": ", getName(wf))
        if (!force) {
            units <- findUnitsTodo(ff, units = units, force = TRUE, 
                cache = FALSE, verbose = less(verbose))
            nbrOfUnits <- length(units)
            verbose && printf(verbose, "%d of the requested units need to be fitted.\n", 
                nbrOfUnits)
        }
        if (length(units) > 0) {
            tRead <- processTime()
            y <- readUnits(wf, units = units, stratifyBy = "pm", 
                ..., force = force, verbose = less(verbose))
            timers$read <- timers$read + (processTime() - tRead)
            tFit <- processTime()
            fit <- base::lapply(y, FUN = fitUnit)
            timers$fit <- timers$fit + (processTime() - tFit)
            y <- NULL
            verbose && str(verbose, fit[1])
            verbose && exit(verbose)
            verbose && enter(verbose, "Storing FIRMA results")
            tWriteFs <- processTime()
            intensities <- unlist(base::lapply(fit, function(unit) {
                base::lapply(unit, function(group) {
                  .subset2(group, "intensities")
                })
            }), use.names = FALSE)
            stdvs <- unlist(base::lapply(fit, function(unit) {
                base::lapply(unit, function(group) {
                  .subset2(group, "stdvs")
                })
            }), use.names = FALSE)
            pixels <- unlist(base::lapply(fit, function(unit) {
                base::lapply(unit, function(group) {
                  .subset2(group, "pixels")
                })
            }), use.names = FALSE)
            rm(fit)
            tGc <- processTime()
            gc <- gc()
            verbose && print(verbose, gc)
            timers$gc <- timers$gc + (processTime() - tGc)
            suppressWarnings({
                map <- getCellMap(ff, units = units, ..., verbose = less(verbose))
            })
            data <- data.frame(intensities = intensities, stdvs = stdvs, 
                pixels = pixels, cell = map[, "cell"])
            rm(intensities, stdvs, pixels, map)
            updateDataFlat(ff, units = units, data = data, verbose = less(verbose))
            rm(data)
            timers$writeFs <- timers$writeFs + (processTime() - 
                tWriteFs)
            verbose && exit(verbose)
            tGc <- processTime()
            gc <- gc()
            verbose && print(verbose, gc)
            timers$gc <- timers$gc + (processTime() - tGc)
            timers$total <- timers$total + (processTime() - tTotal)
        }
        else {
            verbose && exit(verbose)
        }
        clearCache(ff)
        clearCache(wf)
        rm(ff, wf)
    }
    totalTime <- processTime() - startTime
    if (verbose) {
        nunits <- length(units)
        t <- totalTime[3]
        printf(verbose, "Total time for all units across all %d arrays: %.2fs == %.2fmin\n", 
            nbrOfArrays, t, t/60)
        t <- totalTime[3]/nunits
        printf(verbose, "Total time per unit across all %d arrays: %.2fs/unit\n", 
            nbrOfArrays, t)
        t <- totalTime[3]/nunits/nbrOfArrays
        printf(verbose, "Total time per unit and array: %.3gms/unit & array\n", 
            1000 * t)
        t <- nbrOfUnits(cdf) * totalTime[3]/nunits/nbrOfArrays
        printf(verbose, "Total time for one array (%d units): %.2fmin = %.2fh\n", 
            nbrOfUnits(cdf), t/60, t/3600)
        t <- nbrOfUnits(cdf) * totalTime[3]/nunits
        printf(verbose, "Total time for complete data set: %.2fmin = %.2fh\n", 
            t/60, t/3600)
        timers$write <- timers$writeFs
        t <- base::lapply(timers, FUN = function(timer) unname(timer[3]))
        t <- unlist(t)
        t <- 100 * t/t["total"]
        printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%% (of which %.2f%% is for writing results), Explicit garbage collection: %.1f%%\n", 
            t["fit"], t["read"], t["write"], 100 * t["writeFs"]/t["write"], 
            t["gc"])
    }
    invisible(units)
}


calculateResiduals.FirmaModel<-function (this, ...) 
{
    calculateResiduals(this$.plm, ...)
}

setConstructorS3("FirmaModel",function (rmaPlm = NULL, summaryMethod =
"upperQuartile", operateOn="weights",tags = "*", ...) {
   operateOn<- match.arg(operateOn,c("residuals","weights"))
   if (!is.null(tags)) {
       tags <- Arguments$getCharacters(tags)
       tags <- trim(unlist(strsplit(tags, split = ",")))
       asteriskTag <- "FIRMA"
       tags[tags == "*"] <- asteriskTag
       tags <- paste(tags, collapse = ",")
       tags <- unlist(strsplit(tags, split = ","))
   }
   if (!is.character(summaryMethod)) {
       throw("Argument 'summaryMethod' must be a string.")
   }
   extend(UnitModel(..., tags = tags), "FirmaModel", .plm = rmaPlm,
       summaryMethod = summaryMethod, operateOn=operateOn,"cached:.fs" = NULL)
})



getFitUnitFunction.FirmaModel<-function (this,...) 
{
    fitfcn <- getFitFunction(this,, ...)
    fitUnit <- function(unit, ...) {
        base::lapply(unit, FUN = function(group) {
            if (length(group) > 0) {
                y <- .subset2(group, 1)
            }
            else {
                y <- NULL
            }
            y <- log(y, base = 2)
            fitfcn(y)
        })
    }
    fitUnit
}


getFitFunction.FirmaModel<-function (this,...) 
{
  if(this$operateOn=="weights") {
    if (this$summaryMethod == "upperQuartile") {
        fitfcn <- function(y) {
            J <- length(y)
            list(intensities = 2^(1 - quantile(y, probs = 0.75)), 
                stdvs = 1, pixels = 1)
        }
    }
    else if (this$summaryMethod == "median") {
        fitfcn <- function(y) {
            J <- length(y)
            1 - median(y)
            list(intensities = 2^(1 - median(y)), stdvs = 1, 
                pixels = 1)
        }
    }
    else if (this$summaryMethod == "max") {
        fitfcn <- function(y) {
            J <- length(y)
            list(intensities = 2^(1 - max(y)), stdvs = 1, pixels = 1)
        }
    }
    else {
        fitfcn <- function(y) {
            J <- length(y)
            list(intensities = 1, stdvs = 1, pixels = 1)
        }
    }
  } else {
    if (this$summaryMethod == "median") {
        fitfcn <- function(y) {
            J <- length(y)
            list(intensities = 2^median(y), stdvs = 1, pixels = 1)
        }
    }
    else if (this$summaryMethod == "mean") {
        fitfcn <- function(y) {
            J <- length(y)
            list(intensities = 2^mean(y), stdvs = 1, pixels = 1)
        }
    }
  }
    fitfcn
}

