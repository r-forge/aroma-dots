
setMethodS3("extractNuse","QualityAssessmentModel",
function(this,...) {
    plm<-getPlm(this)
    ces<-getChipEffectSet(plm)
    extractNuse(ces,...)
})
setMethodS3("extractRle","QualityAssessmentModel",
function(this,...) {
    plm<-getPlm(this)
    ces<-getChipEffectSet(plm)
    extractRle(ces,...)
})

setMethodS3("plotRle","QualityAssessmentModel",
function(this,...) {
    plotBoxplot(this,type="Rle",...)
})
setMethodS3("plotNuse","QualityAssessmentModel",
function(this,...) {
    plotBoxplot(this,type="Nuse",...)
})
setMethodS3("plotBoxplot","QualityAssessmentModel",
function(this, subset = NULL, verbose = FALSE, arrays=NULL,
    type=c("Nuse","Rle"),main = paste("Boxplots of",type,"Values"),...) 
{
    type<-match.arg(type)
    boxplotStats<-boxplotStats(this,subset = subset, verbose = verbose, type=type,arrays=arrays)
    plotBoxplot(boxplotStats,main=main,...)
})
