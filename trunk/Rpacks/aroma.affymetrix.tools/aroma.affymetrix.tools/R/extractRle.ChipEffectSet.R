setMethodS3("extractRle","ChipEffectSet",
function(this,verbose=FALSE,...) {
    theta <- extractMatrix(this,field="theta",...)
    avg<-getAverageLog(this,field="intensities", mean="median",verbose=verbose)
    thetaR <- extractMatrix(avg,field="theta",...)
    M <- sweep(log2(theta),1,FUN="-",STATS=log2(thetaR))
    return(M)
})
setMethodS3("extractNuse","ChipEffectSet",
function(this,verbose=FALSE,...) {
    theta <- extractMatrix(this, field="sdTheta",...)
    avg<-getAverageLog(this,field="stdvs", mean="median",verbose=verbose)
    thetaR <- extractMatrix(avg,field="theta",...)
    M <- sweep(log2(theta),1,FUN="/",STATS=log2(thetaR))
    return(M)
})
