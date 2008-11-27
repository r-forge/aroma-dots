setMethodS3("importFromBpMap", "AromaCellSequenceFile", function(this,srcPathname,rows=NULL, ...,verbose=TRUE) {
  if( is.null(rows) )
    stop("Must provide the chip dimensions: 'rows' argument is NULL")
  srcPathname <- Arguments$getReadablePathname(srcPathname)
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }
  verbose && enter(verbose,"Updating sequences from ",srcPathname)
  bp<-readBpmap(srcPathname)
  for(i in 1:length(bp)) {
    verbose && enter(verbose,"Updating ",bp[[i]]$seqInfo$fullname[1])
    updateSequences(this,cells=bp[[i]]$pmy*rows+bp[[i]]$pmx+1,seqs=bp[[i]]$probeseq)
    verbose && exit(verbose)
  }
  verbose && exit(verbose)
  invisible(this);
})