
setMethodS3("importFromBpMap", "AromaCellMatchScoreFile", function(this,srcPathname,rows=NULL, ...,verbose=TRUE) {
  if( is.null(rows) )
    stop("Must provide the chip dimensions: 'rows' argument is NULL")
  srcPathname <- Arguments$getReadablePathname(srcPathname)
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }
  verbose && enter(verbose,"Updating sequences from ",srcPathname)
  bp<-readBpmap(srcPathname,readMatchScore=TRUE)
  for(i in 1:length(bp)) {
    verbose && enter(verbose,"Updating ",bp[[i]]$seqInfo$fullname[1])
	ms<-round(bp[[i]]$matchscore*1e6)
	w<-which(ms >= 1 & ms <= 10)
	if(length(w)>0)
      updateMatchScores(this,cells=bp[[i]]$pmy[w]*rows+bp[[i]]$pmx[w]+1,scores=ms[w])
	rm(ms,w);
    verbose && exit(verbose)
  }
  verbose && exit(verbose)
  invisible(this);
})