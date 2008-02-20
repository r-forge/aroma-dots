setMethodS3("plotBoxplot","AffymetrixCelSet",function(this, verbose = FALSE,statsArgs=NULL,...){
    boxplotStats<-do.call("boxplotStats", args=c(list(this=this),statsArgs))
    plotBoxplot(boxplotStats,...)
})
