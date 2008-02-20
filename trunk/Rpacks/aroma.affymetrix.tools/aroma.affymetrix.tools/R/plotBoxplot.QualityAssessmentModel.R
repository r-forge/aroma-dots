##... passed to bxp
setMethodS3("plotBoxplot","QualityAssessmentModel",
function(this, subset = NULL, verbose = FALSE, arrays=NULL,
    type=c("Nuse","Rle"),main = paste("Boxplots of",type,"Values"),...) 
{
    type<-match.arg(type)
    boxplotStats<-boxplotStats(this,subset = subset, verbose = verbose, type=type,arrays=arrays)
    plotBoxplot(boxplotStats,main=main,...)
})
setMethodS3("plotBoxplot","list",
function(boxplotStats,subset = NULL, verbose = FALSE,  
    ylim=NULL,outline=FALSE,las=2,...) 
{
    bxpStats<-convertBoxplotStats(boxplotStats)   
    #fix the strange behavior of bxp if outline=FALSE
    if(!outline && is.null(ylim)) ylim<-range(as.vector(bxpStats[["stats"]]))
    bxp(bxpStats, ylim=ylim,outline=outline,las=las,...)
    invisible(boxplotStats)
})
