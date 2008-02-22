##Convert list of output from boxplot.stat to bxp format...

setMethodS3("plotBoxplot","list",
function(boxplotStats,subset = NULL, verbose = FALSE,  
    ylim=NULL,outline=FALSE,las=2,...) 
{
    bxpStats<-convertBoxplotStats(boxplotStats)   
    #fix the strange behavior of bxp if outline=FALSE
    if(!outline && is.null(ylim)) ylim<-range(as.vector(bxpStats[["stats"]]))
    bxp(bxpStats, ylim=ylim,outline=outline,las=las,...)
    invisible(bxpStats)
})
