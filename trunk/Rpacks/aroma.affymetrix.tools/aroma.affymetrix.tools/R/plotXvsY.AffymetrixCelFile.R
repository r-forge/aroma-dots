setMethodS3("plotXvsY","AffymetrixCelFile",function(this1,this2,indices=NULL, field="intensities",densityPlot=TRUE,loess=TRUE,...)
{
    require(geneplotter)
    nbrOfCells <- nbrOfCells(this1)
    if (is.null(indices)) {
    }
    else {
        indices <- Arguments$getIndices(indices, range = c(1, nbrOfCells))
    }
    if (nbrOfCells != nbrOfCells(this2)) {
        throw("the two CEL files have different number of cells: ", 
            nbrOfCells, " != ", nbrOfCells(this2))
    }
    y1 <- getData(this1, indices = indices, fields = field)[, 1]
    offset <- this1$offset
    if (is.null(offset)) 
        offset <- 0
    if (offset != 0) 
        cat("Offset: ", offset, "\n", sep = "")
    keep1 <- which(y1 != 0)
    y1 <- y1[keep1]
    y1 <- y1 + offset
    y1 <- log(y1, base = 2)
    if (length(y1) == 0) {
        y2 <- y1
    }
    else {
        offset <- this2$offset
        if (is.null(offset)) 
            offset <- 0
        if (offset != 0) 
            cat("Offset: ", offset, "\n", sep = "")
        if (is.null(indices)) {
            indices <- keep1
        }
        else {
            indices <- indices[keep1]
        }
        y2 <- getData(this2, indices = indices, fields = "intensities")[, 1]
        keep2 <- which(y2 != 0)
        if(length(keep2)!=length(y2)) warning("additional -Inf residuals in this2 than in this1")
        y1 <- y1[keep2]
        y2 <- y2[keep2]
        y2 <- y2 + offset
        y2 <- log(y2, base = 2)
    }
    fullData <- cbind(y1,y2)
    colnames(fullData) <- c(getName(this1), getName(this2))
    if(loess) scatter.smooth(fullData,...)
    else{ if(densityPlot) smoothScatter(fullData,...)
        else plot(fullData,...)}
    invisible(fullData)
}
) 
