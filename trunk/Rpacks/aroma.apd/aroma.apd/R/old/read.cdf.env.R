read.cdf.env <- function(fname, complementary.logic = TRUE, verbose = 0) {
    require(Biobase) || stop("Package not found: Biobase")

    fname <- file.path(dirname(fname), basename(fname))
    if (!file.exists(fname))
        stop(paste("file:", fname, "does not exist."))
    pmmm <- .Call("R_affx_get_pmmm_list",
                  as.character(fname),
                  complementary.logic,
                  verbose = verbose, 
                  PACKAGE="affxparser")
    if (is.null(pmmm) || length(pmmm) == 0)
        stop(paste("Error parsing:", fname))
    pmmm <- lapply(pmmm, function(x) {
        tmp <- t(x)
        colnames(tmp) <- c("pm", "mm")
        tmp})
    e <- new.env(hash = TRUE)
    multiassign(names(pmmm), value = pmmm, e)
    return(e)
}

############################################################################
# HISTORY:
# 2006-01-10
# o Added require(Biobase) after removing 'Depends: Biobase' in DESCRIPTION.
# o Added PACKAGE="affxparser" to all .Call() calls. /HB
# o Extracted to its own source file. /HB
############################################################################  
