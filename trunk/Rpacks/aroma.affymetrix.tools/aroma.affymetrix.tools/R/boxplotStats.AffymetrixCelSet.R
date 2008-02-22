setMethodS3("boxplotStats","AffymetrixCelSet",function(this, verbose = FALSE,...){
    nbrOfArrays <- nbrOfArrays(this)
    verbose <- Arguments$getVerbose(verbose)
    arrays<-seq(this)
    cdf <- getCdf(this)
    verbose && enter(verbose, "Calculating summaries for ", length(arrays), " arrays")
    boxplotStats<-list()
    for (cc in 1:length(arrays)) {
        df <- getFile(this, arrays[cc])
        boxplotStats[[cc]]<-boxplotStats(this=df,...)
    }
    names(boxplotStats)<-getNames(this)[arrays]
    verbose && exit(verbose)
    return(boxplotStats)   
})
