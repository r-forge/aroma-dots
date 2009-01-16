bpmapCluster2Cdf <- function(filename, cdfName, nProbes=30, gapDist=3000, rows=NULL,  groupName="Hs", cols=NULL, field="fullname", verbose=10, stringRemove="Hs:Sun.Nov.19.17:25:02.2006;") {
  require("affxparser") || throw("Package not loaded: affxparser");
  require("R.utils") || throw("Package not loaded: R.utils");

  # Argument 'groupName':
  groupName <- Arguments$getCharacter(groupName);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  BpmapUnit2df <- function(u) {
    o <- order(u[["startpos"]]);
    mmx <- u[["mmx"]];
    mmy <- u[["mmy"]];
    mmx <- if (all(mmx == 0) | is.null(mmx)) 0 else mmx;
    mmy <- if (all(mmy == 0) | is.null(mmy)) 0 else mmy;
    data.frame(seqname=u$seqInfo[[field]],groupname=u$seqInfo$groupname, u[c("pmx","pmy")], mmx=mmx, mmy=mmy, u[c("probeseq","strand","startpos","matchscore")], stringsAsFactors=FALSE)[o,];
  } # BpmapUnit2df()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'filename':
  filename <- Arguments$getReadablePathname(filename);

  # Argument 'cdfName':
  cdfName <- Arguments$getCharacter(cdfName);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Generating CDF from BPMAP");

  verbose && enter(verbose, "Reading BPMAP file");
  verbose && cat(verbose, "Source pathname: ", filename);
  bpmaplist <- readBpmap(filename, readMatchScore=TRUE);
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting X/Y locations");
  bpmapdflist <- lapply(bpmaplist, FUN=BpmapUnit2df);
  verbose && exit(verbose);

  rm(bpmaplist); gc()  # save some space

  if (is.null(rows)) {
    rows <- max(sapply(bpmapdflist, FUN=function(u) max(c(u$mmx,u$pmx))));
    verbose && cat(verbose, "NB: 'rows' of CDF are being set as ", rows, 
        ". If this is not correct, stop now and specify 'rows' argument.");
  }

  if (is.null(cols)) {
    cols <- max(sapply(bpmapdflist,FUN=function(u) max(c(u$mmy,u$pmy))));
    verbose && cat(verbose, "NB: 'cols' of CDF are being set as ", cols,
        ". If this is not correct, stop now and specify 'col' argument.");
  }
    
  verbose && enter(verbose, "Creating structure for ", length(bpmapdflist), " units");
  e <- vector("list", 1);
  startps <- l <- ll <- vector("list", 50000);
  count <- 0;
  naValue <- as.character(NA);
  nm <- rep(naValue, length(l));
  for (ii in seq(along=bpmapdflist)) {
    name <- names(bpmapdflist)[ii];

    # Access ones
    bpmapdf <- bpmapdflist[[ii]];
    np <- nrow(bpmapdf);

    ch <- gsub(stringRemove, "", bpmapdf$seqname[1]);

    #if (all(bpmapdf$mmx == 0) & all(bpmapdf$startpos > 0)) {
    if (all(bpmapdf$startpos > 0) & bpmapdf$groupname[1]==groupName) {
      sp <- bpmapdf$startpos;
      d <- diff(sp);
      w <- whichVector(d > gapDist);
      ends <- c(w, np);
      starts <- c(1, w+1);
      k <- whichVector((ends-starts) > nProbes);
      verbose && cat(verbose, length(k), " ROIs for ", name, ".");

      # Access ones
      pmx <- bpmapdf$pmx;
      pmy <- bpmapdf$pmy;
	  mmx <- bpmapdf$mmx;
      mmy <- bpmapdf$mmy;


      for (jj in seq(along=k)) {
        w <- starts[k[jj]]:ends[k[jj]];
        np <- length(w);

		if (all(bpmapdf$mmx==0)) {
		  # PM only
          e[[1]] <- list(x=pmx[w], y=pmy[w], pbase=rep("A", np), tbase=rep("T", np), 
		                 atom=0:(np-1), indexpos=0:(np-1), groupdirection="sense", natoms=np, ncellsperatom=1); 
		} else {
		  # PM+MM
          e[[1]] <- list(x=c(pmx[w],mmx[w]), y=c(pmy[w],mmy[w]), pbase=rep("A", np*2), tbase=rep(c("T","A"), each=np), 
		                 atom=rep(0:(np-1),2), indexpos=rep(0:(np-1),2), groupdirection="sense", natoms=np, ncellsperatom=2); 
		}
        names(e) <- paste(ch, "FROM", sp[starts[k[jj]]], "TO", sp[ends[k[jj]]], sep="");
        na <- sum(unlist(sapply(e,FUN=function(u) u$natoms)));
        nc <- sum(unlist(sapply(e,FUN=function(u) u$natoms*u$ncellsperatom)));
        count <- count + 1;

        l[[count]] <- list(unittype=1, unitdirection=1, groups=e, natoms=na, ncells=nc, ncellsperatom=nc/na, unitnumber=ii);
        startps[[count]] <- sp[w];
        nm[count] <- names(e);
        #if (verbose) { if (count %% 250 == 0) cat(verbose, ".") }
      } # for (jj ...)

      #if (verbose) cat("\n");
    } else {
      # keep all probes
      verbose && cat(verbose, "Skipping all ", np, " probes for ", name, ".");
      next;
    }
  } # for (ii ...)
  verbose && exit(verbose);


  l <- l[1:count];
  names(l) <- nm[1:count];

  verbose && enter(verbose, "Writing PPS file");
  startps <- startps[1:count];
  names(startps) <- names(l);
  saveObject(startps, file=sprintf("%s.pps", cdfName));
  rm(startps);
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing CDF file");
  hdr <- list(probesets=length(l), qcprobesets=0, reference="", chiptype=cdfName, filename=sprintf("%s.cdf", cdfName), nqcunits=0, nunits=length(l), rows=rows, cols=cols, refseq="", nrows=rows, ncols=cols);
  verbose && cat(verbose, "Output pathname: ", hdr$filename);
  verbose && str(verbose, hdr);
  writeCdf(hdr$filename, cdfheader=hdr, cdf=l, cdfqc=NULL, overwrite=TRUE, verbose=verbose);
  verbose && exit(verbose);

  res <- list(cdfList=l, cdfHeader=hdr);

  verbose && exit(verbose);

  invisible(res);
} # BpmapCluster2Cdf()

############################################################################
# HISTORY: 
# 2009-01-14 [MR]
# o fixed 1 small bug ('hdr' out of place)
# o changed the verbose output commands
# 2008-11-28 [HB]
# o Tidying up code.
# o Now using R.utils::Verbose statements.
# 2008-11-xx [MR]
# o Created.
############################################################################





