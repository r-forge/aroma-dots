
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

#setMethodS3("plotRle","QualityAssessmentModel",
#function(this,...) {
#    plotBoxplot(this,type="Rle",...)
#})
#setMethodS3("plotNuse","QualityAssessmentModel",
#function(this,...) {
#    plotBoxplot(this,type="Nuse",...)
#})
